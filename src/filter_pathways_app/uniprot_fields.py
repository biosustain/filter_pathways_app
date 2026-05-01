"""UniProtKB return field definitions.

These display labels mirror the current UniProtKB return-fields help page.
See: https://www.uniprot.org/help/return_fields
"""

# Each entry: field_name -> (api_field_id, description)
UNIPROT_FIELDS: dict[str, tuple[str, str]] = {
    # --- Names & Taxonomy ---
    "Entry": ("accession", "Primary UniProtKB accession"),
    "Entry Name": ("id", "UniProtKB entry name"),
    "Gene Names": (
        "gene_names",
        "All gene names (primary, synonyms, ORF, ordered locus)",
    ),
    "Gene Names (primary)": ("gene_primary", "Primary gene name"),
    "Gene Names (synonym)": ("gene_synonym", "Gene name synonyms"),
    "Gene Names (ordered locus)": (
        "gene_oln",
        "Ordered locus gene names",
    ),
    "Gene Names (ORF)": ("gene_orf", "ORF gene names"),
    "Organism": ("organism_name", "Scientific name of the source organism"),
    "Organism ID": ("organism_id", "NCBI taxonomy identifier"),
    "Protein names": (
        "protein_name",
        "Full protein name including recommended/alternative names",
    ),
    "Proteomes": ("xref_proteomes", "Proteome cross-references"),
    "Taxonomic lineage": ("lineage", "Taxonomic lineage names"),
    "Taxonomic lineage (IDs)": (
        "lineage_ids",
        "Taxonomic lineage identifiers",
    ),
    "Virus hosts": ("virus_hosts", "Virus host organisms"),
    # --- Function ---
    "Absorption": ("absorption", "Absorption comment"),
    "Active site": ("ft_act_site", "Annotated active sites"),
    "Activity regulation": (
        "cc_activity_regulation",
        "Activity regulation comment",
    ),
    "Binding site": ("ft_binding", "Annotated binding sites"),
    "Function [CC]": ("cc_function", "Functional annotation free-text comment"),
    "Catalytic activity": (
        "cc_catalytic_activity",
        "Catalytic reactions catalysed by enzyme",
    ),
    "Cofactor": ("cc_cofactor", "Cofactor annotations"),
    "DNA binding": ("ft_dna_bind", "Annotated DNA-binding regions"),
    "EC number": ("ec", "Enzyme Commission numbers"),
    "Kinetics": ("kinetics", "Enzyme kinetics comment"),
    "Pathway": ("cc_pathway", "Metabolic pathway annotations"),
    "pH dependence": ("ph_dependence", "pH dependence comment"),
    "Redox potential": ("redox_potential", "Redox potential comment"),
    "Rhea ID": ("rhea", "Rhea reaction identifiers"),
    "Site": ("ft_site", "Annotated functional sites"),
    "Temperature dependence": (
        "temp_dependence",
        "Temperature dependence comment",
    ),
    # --- Miscellaneous ---
    "Annotation": ("annotation_score", "Annotation score"),
    "Caution": ("cc_caution", "Curatorial cautions"),
    "Comment Count": ("comment_count", "Number of comment blocks"),
    "Features": ("feature_count", "Number of sequence features"),
    "Keyword ID": ("keywordid", "UniProt keyword identifiers"),
    "Keywords": ("keyword", "UniProt keywords"),
    "Miscellaneous [CC]": (
        "cc_miscellaneous",
        "Miscellaneous comment",
    ),
    "Protein existence": (
        "protein_existence",
        "Protein existence evidence level",
    ),
    "Reviewed": ("reviewed", "Reviewed status"),
    "Tools": ("tools", "External tool references"),
    "UniParc": ("uniparc_id", "UniParc identifier"),
    # --- Interaction ---
    "Interacts with": ("cc_interaction", "Protein interaction comment"),
    "Subunit structure [CC]": ("cc_subunit", "Subunit structure comment"),
    # --- Expression ---
    "Developmental stage": (
        "cc_developmental_stage",
        "Developmental stage comment",
    ),
    "Induction": ("cc_induction", "Induction comment"),
    "Tissue specificity": (
        "cc_tissue_specificity",
        "Tissue specificity comment",
    ),
    # --- Gene Ontology ---
    "Gene Ontology (biological process)": (
        "go_p",
        "Gene Ontology biological process terms",
    ),
    "Gene Ontology (cellular component)": (
        "go_c",
        "Gene Ontology cellular component terms",
    ),
    "Gene Ontology (GO)": ("go", "All Gene Ontology terms"),
    "Gene Ontology (molecular function)": (
        "go_f",
        "Gene Ontology molecular function terms",
    ),
    "Gene Ontology IDs": ("go_id", "Gene Ontology identifiers"),
    # --- Pathology & Biotech ---
    "Allergenic properties": (
        "cc_allergen",
        "Allergenic properties comment",
    ),
    "Biotechnological use": (
        "cc_biotechnology",
        "Biotechnological use comment",
    ),
    "Disruption phenotype": (
        "cc_disruption_phenotype",
        "Disruption phenotype comment",
    ),
    "Involvement in disease": ("cc_disease", "Disease associations"),
    "Mutagenesis": ("ft_mutagen", "Annotated mutagenesis sites"),
    "Pharmaceutical use": (
        "cc_pharmaceutical",
        "Pharmaceutical use comment",
    ),
    "Toxic dose": ("cc_toxic_dose", "Toxic dose comment"),
    # --- Cross-references ---
    "Reactome": ("xref_reactome", "Reactome pathway cross-references"),
    "KEGG": ("xref_kegg", "KEGG pathway cross-references"),
    "BioCyc": ("xref_biocyc", "BioCyc pathway cross-references"),
    "WikiPathways": ("xref_wikipathways", "WikiPathways cross-references"),
    # --- Subcellular location ---
    "Intramembrane": ("ft_intramem", "Annotated intramembrane regions"),
    "Subcellular location [CC]": (
        "cc_subcellular_location",
        "Subcellular localisation annotation",
    ),
    "Topological domain": ("ft_topo_dom", "Annotated topological domains"),
    "Transmembrane": ("ft_transmem", "Annotated transmembrane regions"),
    # --- PTM / Processing ---
    "Chain": ("ft_chain", "Annotated chain features"),
    "Cross-link": ("ft_crosslnk", "Annotated cross-links"),
    "Disulfide bond": ("ft_disulfid", "Annotated disulfide bonds"),
    "Glycosylation": ("ft_carbohyd", "Annotated glycosylation sites"),
    "Initiator methionine": (
        "ft_init_met",
        "Annotated initiator methionines",
    ),
    "Lipidation": ("ft_lipid", "Annotated lipidation sites"),
    "Modified residue": ("ft_mod_res", "Annotated modified residues"),
    "Peptide": ("ft_peptide", "Annotated peptide features"),
    "Post-translational modification": ("cc_ptm", "PTM annotations"),
    "Propeptide": ("ft_propep", "Annotated propeptides"),
    "Signal peptide": ("ft_signal", "Annotated signal peptides"),
    "Transit peptide": ("ft_transit", "Annotated transit peptides"),
    # --- Structure ---
    "3D": ("structure_3d", "3D structure availability"),
    "Beta strand": ("ft_strand", "Annotated beta strands"),
    "Helix": ("ft_helix", "Annotated helices"),
    "Turn": ("ft_turn", "Annotated turns"),
    # --- Publications ---
    "PubMed ID": ("lit_pubmed_id", "PubMed identifiers"),
    # --- Date of ---
    "Date of creation": ("date_created", "Entry creation date"),
    "Date of last modification": (
        "date_modified",
        "Last entry modification date",
    ),
    "Date of last sequence modification": (
        "date_sequence_modified",
        "Last sequence modification date",
    ),
    "Entry version": ("version", "Entry version number"),
    # --- Family & Domains ---
    "Coiled coil": ("ft_coiled", "Annotated coiled-coil regions"),
    "Compositional bias": (
        "ft_compbias",
        "Annotated compositional bias regions",
    ),
    "Domain[CC]": ("cc_domain", "Domain comment"),
    "Domain[FT]": ("ft_domain", "Annotated domains"),
    "Motif": ("ft_motif", "Annotated motifs"),
    "Protein families": ("protein_families", "Protein family assignments"),
    "Region": ("ft_region", "Annotated regions"),
    "Repeat": ("ft_repeat", "Annotated repeats"),
    "Zinc finger": ("ft_zn_fing", "Annotated zinc finger regions"),
    # --- Sequences ---
    "Alternative products": (
        "cc_alternative_products",
        "Alternative products comment",
    ),
    "Alternative sequence": ("ft_var_seq", "Annotated alternative sequences"),
    "Erroneous gene model prediction": (
        "cc_sc_epred",
        "Erroneous gene model prediction comment",
    ),
    "Fragment": ("fragment", "Fragment status"),
    "Gene encoded in": ("encoded_in", "Location of the encoded gene"),
    "Length": ("length", "Length of the canonical protein sequence"),
    "Mass": ("mass", "Molecular mass of the canonical sequence"),
    "Mass spectrometry": (
        "cc_mass_spectrometry",
        "Mass spectrometry comment",
    ),
    "Natural variant": ("ft_variant", "Annotated natural variants"),
    "Non-adjacent residues": (
        "ft_non_cons",
        "Annotated non-adjacent residues",
    ),
    "Non-standard residue": (
        "ft_non_std",
        "Annotated non-standard residues",
    ),
    "Non-terminal residue": (
        "ft_non_ter",
        "Annotated non-terminal residues",
    ),
    "Polymorphism": ("cc_polymorphism", "Polymorphism comment"),
    "RNA editing": ("cc_rna_editing", "RNA editing comment"),
    "Sequence": ("sequence", "Canonical protein sequence"),
    "Sequence caution": (
        "cc_sequence_caution",
        "Sequence caution comment",
    ),
    "Sequence conflict": ("ft_conflict", "Annotated sequence conflicts"),
    "Sequence uncertainty": (
        "ft_unsure",
        "Annotated sequence uncertainties",
    ),
    "Sequence version": ("sequence_version", "Sequence version number"),
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
    >>> fields_to_api_string(["Gene Ontology (biological process)"])
    'accession,go_p'
    >>> fields_to_api_string(["Entry", "Pathway"])
    'accession,cc_pathway'
    >>> fields_to_api_string([])
    'accession'
    """
    api_ids = ["accession"]
    seen_api_ids = set(api_ids)
    for name in display_names:
        api_id, _ = UNIPROT_FIELDS[name]
        if api_id not in seen_api_ids:
            api_ids.append(api_id)
            seen_api_ids.add(api_id)
    return ",".join(api_ids)
