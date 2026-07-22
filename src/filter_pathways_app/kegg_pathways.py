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

The same KEGG gene IDs also link to KEGG Orthology (KO) terms and Enzyme
Commission (EC) numbers via the KEGG ``link`` endpoint.  KO terms enable
cross-species lookup; EC numbers are what KEGG pathway maps are annotated with.
:func:`query_kegg_terms` fetches these and returns them in the same long format.
"""

from __future__ import annotations

import re as _re
from collections import defaultdict

import pandas as pd
import requests
from acore.io.kegg import fetch_kegg_ko_descriptions, link_kegg_batch
from acore.io.uniprot import fetch_annotations

# Source label used for the KEGG pathway rows in the long-format DataFrame.
KEGG_PATHWAY_SOURCE = "KEGG Pathway"

# Source labels for the KEGG orthology and enzyme rows.  KO terms enable
# cross-species lookup; EC numbers are what KEGG pathway maps are annotated with.
KEGG_KO_SOURCE = "KEGG Orthology (KO)"
KEGG_EC_SOURCE = "KEGG EC number"

# Source label of the raw UniProt ``xref_kegg`` rows in the long-format
# DataFrame (matches the "KEGG" display name in ``UNIPROT_FIELDS``).
KEGG_XREF_SOURCE = "KEGG"

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


def parse_gene_links(raw: str, prefix: str = "") -> dict[str, list[str]]:
    """Parse the tab-delimited output of a KEGG ``link/<target>`` endpoint.

    All ``link`` endpoints share the ``gene<TAB>target`` shape, differing only
    in the target prefix (``path:``, ``ko:``, ``ec:``).

    Parameters
    ----------
    raw : str
        Raw response text, one ``gene<TAB>prefix:term`` pair per line.
    prefix : str, optional
        Prefix to strip from each target term (e.g. ``"path:"``, ``"ko:"``,
        ``"ec:"``).  Default strips nothing.

    Returns
    -------
    dict[str, list[str]]
        Mapping of KEGG gene ID to its (de-duplicated, order-preserving) target
        terms.

    Examples
    --------
    >>> parse_gene_links("hsa:351\\tko:K04520\\nhsa:351\\tko:K04520\\n", "ko:")
    {'hsa:351': ['K04520']}
    >>> parse_gene_links("yli:2912002\\tec:2.3.1.16\\n", "ec:")
    {'yli:2912002': ['2.3.1.16']}
    >>> parse_gene_links("")
    {}
    """
    gene_to_terms: defaultdict[str, list[str]] = defaultdict(list)
    for line in raw.strip().splitlines():
        line = line.strip()
        if not line or "\t" not in line:
            continue
        gene_id, term = line.split("\t", maxsplit=1)
        term = term.strip().removeprefix(prefix)
        if term not in gene_to_terms[gene_id]:
            gene_to_terms[gene_id].append(term)
    return dict(gene_to_terms)


def parse_gene_pathway_links(raw: str) -> dict[str, list[str]]:
    """Parse the tab-delimited output of the KEGG ``link/pathway`` endpoint.

    Thin wrapper over :func:`parse_gene_links` that strips the ``path:`` prefix.

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
    return parse_gene_links(raw, "path:")


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


def format_term_annotation(term_id: str, name: str = "") -> str:
    """Format a KEGG term as an annotation string ``"name [term_id]"``.

    Mirrors the ``"term [GO:...]"`` style used for Gene Ontology annotations so
    KEGG pathways, KO terms and EC numbers read consistently in the results
    table.  Falls back to the bare identifier when no name is available (as for
    EC numbers, which carry no separate name here).

    Examples
    --------
    >>> format_term_annotation("K04520", "amyloid beta precursor protein")
    'amyloid beta precursor protein [K04520]'
    >>> format_term_annotation("2.3.1.16")
    '2.3.1.16'
    """
    name = name.strip()
    return f"{name} [{term_id}]" if name else term_id


def format_pathway_annotation(pathway_id: str, name: str) -> str:
    """Format a KEGG pathway as an annotation string ``"name [pathway_id]"``.

    Thin wrapper over :func:`format_term_annotation`.

    Examples
    --------
    >>> format_pathway_annotation("hsa00010", "Glycolysis / Gluconeogenesis")
    'Glycolysis / Gluconeogenesis [hsa00010]'
    >>> format_pathway_annotation("hsa99999", "")
    'hsa99999'
    """
    return format_term_annotation(pathway_id, name)


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


def gene_map_from_annotations(annotations: pd.DataFrame) -> dict[str, list[str]]:
    """Build the UniProt→KEGG-gene map from an already-fetched annotations DF.

    When the ``KEGG`` return field was requested, :func:`query_uniprot` already
    contains rows whose ``annotation`` is the raw ``xref_kegg`` value (e.g.
    ``"yli:2912002;"``).  Parsing those rows yields the same mapping as
    :func:`fetch_kegg_gene_map` without a second UniProt request.

    Parameters
    ----------
    annotations : pd.DataFrame
        Long-format DataFrame with ``identifier``, ``source`` and
        ``annotation`` columns, as returned by
        :func:`filter_pathways_app.filter_pathways.query_uniprot`.

    Returns
    -------
    dict[str, list[str]]
        Mapping of UniProt accession to its KEGG gene IDs.  Empty when the
        DataFrame carries no ``KEGG`` cross-reference rows.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     "identifier": ["Q05493", "Q05493", "B5FVB1"],
    ...     "source": ["KEGG", "Entry", "KEGG"],
    ...     "annotation": ["yli:2912002;", "Q05493", "yli:7009397;"],
    ... })
    >>> gene_map_from_annotations(df)
    {'Q05493': ['yli:2912002'], 'B5FVB1': ['yli:7009397']}
    """
    if annotations is None or annotations.empty:
        return {}
    kegg_rows = annotations[annotations["source"] == KEGG_XREF_SOURCE]
    gene_map: dict[str, list[str]] = {}
    for identifier, group in kegg_rows.groupby("identifier", sort=False):
        genes = list(
            dict.fromkeys(
                gene
                for value in group["annotation"]
                for gene in parse_kegg_gene_ids(value)
            )
        )
        if genes:
            gene_map[identifier] = genes
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


def fetch_gene_ko(gene_ids: list[str]) -> dict[str, list[str]]:
    """Map KEGG gene IDs to their KEGG Orthology (KO) terms.

    KO terms group orthologous genes across organisms, so they are the handle
    for cross-species lookup.  Uses the KEGG ``link/ko`` endpoint.

    Parameters
    ----------
    gene_ids : list[str]
        KEGG gene IDs, e.g. ``["hsa:351", "yli:2912002"]``.

    Returns
    -------
    dict[str, list[str]]
        Mapping of KEGG gene ID to its KO terms (``ko:`` prefix stripped, e.g.
        ``"K04520"``).
    """
    unique_ids = list(dict.fromkeys(gene_ids))
    if not unique_ids:
        return {}
    raw = link_kegg_batch("ko", unique_ids)
    return parse_gene_links(raw, "ko:")


def fetch_gene_ec(gene_ids: list[str]) -> dict[str, list[str]]:
    """Map KEGG gene IDs to their Enzyme Commission (EC) numbers.

    EC numbers are what KEGG annotates its pathway/reaction maps with.  Uses the
    KEGG ``link/enzyme`` endpoint.  Non-enzyme genes simply yield no EC numbers.

    Parameters
    ----------
    gene_ids : list[str]
        KEGG gene IDs, e.g. ``["hsa:351", "yli:2912002"]``.

    Returns
    -------
    dict[str, list[str]]
        Mapping of KEGG gene ID to its EC numbers (``ec:`` prefix stripped, e.g.
        ``"2.3.1.16"``).  Genes without an EC number are omitted.
    """
    unique_ids = list(dict.fromkeys(gene_ids))
    if not unique_ids:
        return {}
    raw = link_kegg_batch("enzyme", unique_ids)
    return parse_gene_links(raw, "ec:")


def fetch_ko_descriptions(ko_terms: list[str]) -> dict[str, str]:
    """Fetch human-readable names for KEGG KO terms.

    Wraps :func:`acore.io.kegg.fetch_kegg_ko_descriptions` and keys the result
    by the bare KO identifier (``ko:`` prefix stripped) to match the terms
    produced by :func:`fetch_gene_ko`.

    Parameters
    ----------
    ko_terms : list[str]
        KO identifiers such as ``["K04520", "K07513"]``.

    Returns
    -------
    dict[str, str]
        Mapping of KO term to its common description.  Terms whose description
        could not be resolved are omitted.
    """
    unique_terms = list(dict.fromkeys(ko_terms))
    if not unique_terms:
        return {}
    described = fetch_kegg_ko_descriptions(unique_terms)
    names: dict[str, str] = {}
    for _, row in described.iterrows():
        term = str(row["ko_term"]).removeprefix("ko:")
        description = str(row["common_description"]).strip()
        if description:
            names[term] = description
    return names


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
def _build_term_rows(
    gene_map: dict[str, list[str]],
    gene_terms: dict[str, list[str]],
    source: str,
    names: dict[str, str] | None = None,
) -> list[dict[str, str]]:
    """Build long-format rows mapping each protein to its (deduplicated) terms.

    Shared by the pathway/KO/EC orchestration: for every protein, walk its KEGG
    genes, collect the terms each gene links to, drop duplicates across genes,
    and format them into ``identifier``/``source``/``annotation`` rows.
    """
    names = names or {}
    rows: list[dict[str, str]] = []
    for identifier, genes in gene_map.items():
        seen: set[str] = set()
        for gene in genes:
            for term in gene_terms.get(gene, []):
                if term in seen:
                    continue
                seen.add(term)
                rows.append(
                    {
                        "identifier": identifier,
                        "source": source,
                        "annotation": format_term_annotation(
                            term, names.get(term, "")
                        ),
                    }
                )
    return rows


def query_kegg_pathways(
    uniprot_ids: list[str] | pd.Index,
    gene_map: dict[str, list[str]] | None = None,
) -> pd.DataFrame:
    """Fetch the KEGG pathways associated with the given UniProt proteins.

    Combines the three steps (UniProt ``xref_kegg`` → KEGG genes → KEGG pathways
    → pathway names) into a single long-format DataFrame that matches the shape
    returned by :func:`filter_pathways_app.filter_pathways.query_uniprot`.

    Parameters
    ----------
    uniprot_ids : list[str] | pd.Index
        UniProt accession IDs.
    gene_map : dict[str, list[str]] | None, optional
        Pre-resolved UniProt→KEGG-gene mapping (e.g. from
        :func:`gene_map_from_annotations` when the ``KEGG`` field was already
        fetched).  When provided, the extra UniProt ``xref_kegg`` request is
        skipped; when ``None`` (the default) it is fetched via
        :func:`fetch_kegg_gene_map`.

    Returns
    -------
    pd.DataFrame
        Long-format DataFrame with columns ``identifier``, ``source`` (always
        ``"KEGG Pathway"``) and ``annotation`` (``"name [pathway_id]"``).
        Empty (with those columns) when no pathways are found.
    """
    columns = ["identifier", "source", "annotation"]
    if gene_map is None:
        gene_map = fetch_kegg_gene_map(uniprot_ids)
    if not gene_map:
        return pd.DataFrame(columns=columns)

    all_genes = sorted({gene for genes in gene_map.values() for gene in genes})
    gene_pathways = fetch_gene_pathways(all_genes)

    all_pathways = sorted({pid for pids in gene_pathways.values() for pid in pids})
    if not all_pathways:
        return pd.DataFrame(columns=columns)

    pathway_names = fetch_pathway_names(all_pathways)
    rows = _build_term_rows(
        gene_map, gene_pathways, KEGG_PATHWAY_SOURCE, pathway_names
    )
    return pd.DataFrame(rows, columns=columns)


def query_kegg_terms(
    uniprot_ids: list[str] | pd.Index,
    gene_map: dict[str, list[str]] | None = None,
    include_ko: bool = True,
    include_ec: bool = True,
) -> pd.DataFrame:
    """Fetch the KEGG Orthology (KO) and Enzyme (EC) terms for UniProt proteins.

    Resolves each protein's KEGG gene(s) to KO terms (for cross-species lookup)
    and EC numbers (used to annotate KEGG pathway maps), returning them in the
    same long format as :func:`query_kegg_pathways`.

    Parameters
    ----------
    uniprot_ids : list[str] | pd.Index
        UniProt accession IDs.
    gene_map : dict[str, list[str]] | None, optional
        Pre-resolved UniProt→KEGG-gene mapping (see :func:`query_kegg_pathways`).
        When ``None`` it is fetched via :func:`fetch_kegg_gene_map`.
    include_ko, include_ec : bool, optional
        Whether to include KO terms / EC numbers.  Both default to ``True``.

    Returns
    -------
    pd.DataFrame
        Long-format DataFrame with columns ``identifier``, ``source``
        (``"KEGG Orthology (KO)"`` or ``"KEGG EC number"``) and ``annotation``.
        Empty (with those columns) when nothing is found or both flags are off.
    """
    columns = ["identifier", "source", "annotation"]
    if gene_map is None:
        gene_map = fetch_kegg_gene_map(uniprot_ids)
    if not gene_map or not (include_ko or include_ec):
        return pd.DataFrame(columns=columns)

    all_genes = sorted({gene for genes in gene_map.values() for gene in genes})
    rows: list[dict[str, str]] = []

    if include_ko:
        gene_ko = fetch_gene_ko(all_genes)
        all_ko = sorted({ko for kos in gene_ko.values() for ko in kos})
        ko_names = fetch_ko_descriptions(all_ko)
        rows += _build_term_rows(gene_map, gene_ko, KEGG_KO_SOURCE, ko_names)

    if include_ec:
        gene_ec = fetch_gene_ec(all_genes)
        rows += _build_term_rows(gene_map, gene_ec, KEGG_EC_SOURCE)

    return pd.DataFrame(rows, columns=columns)
