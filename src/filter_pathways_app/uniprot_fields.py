"""UniProt return field definitions.

A curated subset of fields available from the UniProt API.
See: https://www.uniprot.org/help/return_fields
"""

# Each entry: field_name -> (api_field_id, description)
UNIPROT_FIELDS: dict[str, tuple[str, str]] = {
    # --- Names & Taxonomy ---
    "Gene names": (
        "gene_names",
        "All gene names (primary, synonyms, ORF, ordered locus)",
    ),
    "Organism": ("organism_name", "Scientific name of the source organism"),
    "Protein names": (
        "protein_name",
        "Full protein name including recommended/alternative names",
    ),
    # --- Function ---
    "Function [CC]": ("cc_function", "Functional annotation free-text comment"),
    "Pathway [CC]": ("cc_pathway", "Metabolic pathway annotations"),
    "Involvement in disease": ("cc_disease", "Disease associations"),
    "Catalytic activity": (
        "cc_catalytic_activity",
        "Catalytic reactions catalysed by enzyme",
    ),
    # --- Gene Ontology ---
    "GO (biological process)": ("go_p", "Gene Ontology biological process terms"),
    "GO (molecular function)": ("go_f", "Gene Ontology molecular function terms"),
    "GO (cellular component)": ("go_c", "Gene Ontology cellular component terms"),
    # --- Pathways & Interactions ---
    "Reactome": ("xref_reactome", "Reactome pathway cross-references"),
    "KEGG": ("xref_kegg", "KEGG pathway cross-references"),
    "BioCyc": ("xref_biocyc", "BioCyc pathway cross-references"),
    "WikiPathways": ("xref_wikipathways", "WikiPathways cross-references"),
    # --- Subcellular location ---
    "Subcellular location [CC]": (
        "cc_subcellular_location",
        "Subcellular localisation annotation",
    ),
    # --- Post-translational modification ---
    "Post-translational modification": ("cc_ptm", "PTM annotations"),
    # --- Sequence ---
    "Sequence length": ("length", "Length of the canonical protein sequence"),
    "Mass": ("mass", "Molecular mass of the canonical sequence"),
}

# Ordered list of display names (for UI ordering)
FIELD_DISPLAY_NAMES: list[str] = list(UNIPROT_FIELDS.keys())


def fields_to_api_string(display_names: list[str]) -> str:
    """Convert a list of display names to an API fields string.

    Parameters
    ----------
    display_names : list[str]
        Display names of the fields to include (keys from UNIPROT_FIELDS).

    Returns
    -------
    str
        Comma-separated API field identifiers suitable for the UniProt API.
        The ``accession`` field is always prepended.

    Examples
    --------
    >>> fields_to_api_string(["GO (biological process)"])
    'accession,go_p'
    >>> fields_to_api_string([])
    'accession'
    """
    api_ids = ["accession"]
    for name in display_names:
        api_id, _ = UNIPROT_FIELDS[name]
        api_ids.append(api_id)
    return ",".join(api_ids)
