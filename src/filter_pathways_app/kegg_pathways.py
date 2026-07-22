"""Enrich UniProt proteins with KEGG pathway annotations.

UniProt's ``xref_kegg`` cross-reference maps a protein to KEGG *gene* IDs
(e.g. ``hsa:351``).  KEGG genes in turn participate in KEGG *pathways*.  This
module bridges that gap so the app can show, for each protein, the KEGG
pathways it belongs to:

1. Fetch the ``xref_kegg`` cross-reference for the UniProt IDs and parse out the
   KEGG gene IDs (:func:`fetch_kegg_gene_map`).
2. Resolve which pathways each gene participates in, using
   :func:`acore.io.kegg.link_kegg_batch` (:func:`fetch_gene_pathways`).
3. Fill in the human-readable pathway names via the KEGG ``list`` endpoint
   (:func:`fetch_pathway_names`).

:func:`query_kegg_pathways` ties these together and returns rows in the same
long format (``identifier``, ``source``, ``annotation``) used elsewhere in the
app, so KEGG pathways can be filtered, pivoted and exported like any other
annotation source.
"""

from __future__ import annotations

import re as _re
from collections import defaultdict

import pandas as pd
import requests
from acore.io.kegg import link_kegg_batch
from acore.io.uniprot import fetch_annotations

# Source label used for the KEGG pathway rows in the long-format DataFrame.
KEGG_PATHWAY_SOURCE = "KEGG Pathway"

KEGG_API_BASE_URL = "https://rest.kegg.jp"

# A KEGG gene ID is an organism code (2-4 letters) followed by ":" and a locus,
# e.g. "hsa:351" or "yli:2912002".
_KEGG_GENE_PATTERN = _re.compile(r"[A-Za-z]{2,4}:[^\s;]+")

# The organism prefix of a KEGG pathway ID, e.g. "hsa" in "hsa00010".
_PATHWAY_ORG_PATTERN = _re.compile(r"^([A-Za-z]+)\d+$")


# ---------------------------------------------------------------------------
# Pure parsing helpers (network-free, doctested)
# ---------------------------------------------------------------------------
def parse_kegg_gene_ids(value: str) -> list[str]:
    """Parse the KEGG gene IDs out of a UniProt ``xref_kegg`` field value.

    Parameters
    ----------
    value : str
        Raw ``xref_kegg`` value, e.g. ``"hsa:351;"`` or
        ``"yli:2912002;yli:2910294;"``.

    Returns
    -------
    list[str]
        The KEGG gene IDs, in order of appearance and de-duplicated.

    Examples
    --------
    >>> parse_kegg_gene_ids("hsa:351;")
    ['hsa:351']
    >>> parse_kegg_gene_ids("yli:2912002;yli:2910294;")
    ['yli:2912002', 'yli:2910294']
    >>> parse_kegg_gene_ids("")
    []
    """
    if not isinstance(value, str):
        return []
    return list(dict.fromkeys(_KEGG_GENE_PATTERN.findall(value)))


def parse_gene_pathway_links(raw: str) -> dict[str, list[str]]:
    """Parse the tab-delimited output of the KEGG ``link/pathway`` endpoint.

    Parameters
    ----------
    raw : str
        Raw response text, one ``gene<TAB>path:pathway`` pair per line, e.g.::

            hsa:351\\tpath:hsa00010
            hsa:351\\tpath:hsa04930

    Returns
    -------
    dict[str, list[str]]
        Mapping of KEGG gene ID to the pathway IDs it belongs to.  The
        ``path:`` prefix is stripped from the pathway IDs.

    Examples
    --------
    >>> parse_gene_pathway_links("hsa:351\\tpath:hsa00010\\nhsa:351\\tpath:hsa04930\\n")
    {'hsa:351': ['hsa00010', 'hsa04930']}
    >>> parse_gene_pathway_links("")
    {}
    """
    gene_to_pathways: defaultdict[str, list[str]] = defaultdict(list)
    for line in raw.strip().splitlines():
        line = line.strip()
        if not line or "\t" not in line:
            continue
        gene_id, pathway_id = line.split("\t", maxsplit=1)
        pathway_id = pathway_id.strip().removeprefix("path:")
        if pathway_id not in gene_to_pathways[gene_id]:
            gene_to_pathways[gene_id].append(pathway_id)
    return dict(gene_to_pathways)


def parse_pathway_names(raw: str) -> dict[str, str]:
    """Parse the tab-delimited output of the KEGG ``list/pathway`` endpoint.

    Parameters
    ----------
    raw : str
        Raw response text, one ``pathway<TAB>name - organism`` pair per line,
        e.g. ``"hsa00010\\tGlycolysis / Gluconeogenesis - Homo sapiens (human)"``.

    Returns
    -------
    dict[str, str]
        Mapping of pathway ID to its (organism-stripped) name.

    Examples
    --------
    >>> parse_pathway_names(
    ...     "hsa00010\\tGlycolysis / Gluconeogenesis - Homo sapiens (human)\\n"
    ... )
    {'hsa00010': 'Glycolysis / Gluconeogenesis'}
    """
    names: dict[str, str] = {}
    for line in raw.strip().splitlines():
        line = line.strip()
        if not line or "\t" not in line:
            continue
        pathway_id, name = line.split("\t", maxsplit=1)
        names[pathway_id.strip()] = _strip_organism(name.strip())
    return names


def _strip_organism(name: str) -> str:
    """Strip the trailing ``" - <organism>"`` KEGG appends to pathway names.

    Examples
    --------
    >>> _strip_organism("Glycolysis / Gluconeogenesis - Homo sapiens (human)")
    'Glycolysis / Gluconeogenesis'
    >>> _strip_organism("Metabolic pathways")
    'Metabolic pathways'
    """
    return name.rsplit(" - ", maxsplit=1)[0]


def pathway_organism(pathway_id: str) -> str | None:
    """Return the organism code of a KEGG pathway ID, or ``None``.

    Examples
    --------
    >>> pathway_organism("hsa00010")
    'hsa'
    >>> pathway_organism("yli00592")
    'yli'
    >>> pathway_organism("garbage") is None
    True
    """
    match = _PATHWAY_ORG_PATTERN.match(pathway_id.strip())
    return match.group(1) if match else None


def format_pathway_annotation(pathway_id: str, name: str) -> str:
    """Format a KEGG pathway as an annotation string ``"name [pathway_id]"``.

    Mirrors the ``"term [GO:...]"`` style used for Gene Ontology annotations so
    KEGG pathways read consistently in the results table.

    Examples
    --------
    >>> format_pathway_annotation("hsa00010", "Glycolysis / Gluconeogenesis")
    'Glycolysis / Gluconeogenesis [hsa00010]'
    >>> format_pathway_annotation("hsa99999", "")
    'hsa99999'
    """
    name = name.strip()
    return f"{name} [{pathway_id}]" if name else pathway_id


# ---------------------------------------------------------------------------
# Network-backed helpers
# ---------------------------------------------------------------------------
def fetch_kegg_gene_map(uniprot_ids: list[str] | pd.Index) -> dict[str, list[str]]:
    """Fetch the ``xref_kegg`` cross-reference and map proteins to KEGG genes.

    Parameters
    ----------
    uniprot_ids : list[str] | pd.Index
        UniProt accession IDs.

    Returns
    -------
    dict[str, list[str]]
        Mapping of UniProt accession to its KEGG gene IDs.  Proteins without a
        KEGG cross-reference are omitted.
    """
    raw = fetch_annotations(uniprot_ids, fields="accession,xref_kegg")
    # Columns are ["From", "Entry", "KEGG"]; the KEGG gene IDs are the last one.
    gene_map: dict[str, list[str]] = {}
    for _, row in raw.iterrows():
        genes = parse_kegg_gene_ids(row.iloc[-1])
        if genes:
            gene_map[row["From"]] = genes
    return gene_map


def fetch_gene_pathways(gene_ids: list[str]) -> dict[str, list[str]]:
    """Map KEGG gene IDs to the pathways they participate in.

    Uses :func:`acore.io.kegg.link_kegg_batch` (KEGG ``link/pathway`` endpoint),
    which batches requests to respect the KEGG API limits.

    Parameters
    ----------
    gene_ids : list[str]
        KEGG gene IDs, e.g. ``["hsa:351", "yli:2912002"]``.

    Returns
    -------
    dict[str, list[str]]
        Mapping of KEGG gene ID to its pathway IDs (``path:`` prefix stripped).
    """
    unique_ids = list(dict.fromkeys(gene_ids))
    if not unique_ids:
        return {}
    raw = link_kegg_batch("pathway", unique_ids)
    return parse_gene_pathway_links(raw)


def fetch_pathway_names(pathway_ids: list[str]) -> dict[str, str]:
    """Fetch human-readable names for KEGG pathway IDs.

    Pathway names are organism-specific, so this issues one KEGG ``list/pathway``
    request per organism (rather than one per pathway) and builds a lookup.

    Parameters
    ----------
    pathway_ids : list[str]
        Pathway IDs such as ``["hsa00010", "hsa04930"]``.

    Returns
    -------
    dict[str, str]
        Mapping of pathway ID to name.  Pathways whose name could not be
        resolved are omitted.
    """
    organisms = {
        org for org in (pathway_organism(pid) for pid in pathway_ids) if org is not None
    }
    names: dict[str, str] = {}
    for org in sorted(organisms):
        response = requests.get(f"{KEGG_API_BASE_URL}/list/pathway/{org}", timeout=30)
        response.raise_for_status()
        names.update(parse_pathway_names(response.text))
    return names


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def query_kegg_pathways(uniprot_ids: list[str] | pd.Index) -> pd.DataFrame:
    """Fetch the KEGG pathways associated with the given UniProt proteins.

    Combines the three steps (UniProt ``xref_kegg`` → KEGG genes → KEGG pathways
    → pathway names) into a single long-format DataFrame that matches the shape
    returned by :func:`filter_pathways_app.filter_pathways.query_uniprot`.

    Parameters
    ----------
    uniprot_ids : list[str] | pd.Index
        UniProt accession IDs.

    Returns
    -------
    pd.DataFrame
        Long-format DataFrame with columns ``identifier``, ``source`` (always
        ``"KEGG Pathway"``) and ``annotation`` (``"name [pathway_id]"``).
        Empty (with those columns) when no pathways are found.
    """
    columns = ["identifier", "source", "annotation"]
    gene_map = fetch_kegg_gene_map(uniprot_ids)
    if not gene_map:
        return pd.DataFrame(columns=columns)

    all_genes = sorted({gene for genes in gene_map.values() for gene in genes})
    gene_pathways = fetch_gene_pathways(all_genes)

    all_pathways = sorted({pid for pids in gene_pathways.values() for pid in pids})
    if not all_pathways:
        return pd.DataFrame(columns=columns)

    pathway_names = fetch_pathway_names(all_pathways)

    rows: list[dict[str, str]] = []
    for identifier, genes in gene_map.items():
        # De-duplicate pathways across the protein's (possibly multiple) genes.
        seen: set[str] = set()
        for gene in genes:
            for pathway_id in gene_pathways.get(gene, []):
                if pathway_id in seen:
                    continue
                seen.add(pathway_id)
                rows.append(
                    {
                        "identifier": identifier,
                        "source": KEGG_PATHWAY_SOURCE,
                        "annotation": format_pathway_annotation(
                            pathway_id, pathway_names.get(pathway_id, "")
                        ),
                    }
                )
    return pd.DataFrame(rows, columns=columns)
