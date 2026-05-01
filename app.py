"""Streamlit app for querying UniProt and filtering pathway/annotation data.

Usage
-----
Run from the repository root::

    streamlit run app.py
"""

import streamlit as st

from filter_pathways_app.filter_pathways import (
    export_to_csv,
    filter_annotations,
    query_uniprot,
)
from filter_pathways_app.uniprot_fields import (
    FIELD_DISPLAY_NAMES,
    UNIPROT_FIELDS,
    fields_to_api_string,
)

# region Page config
# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="UniProt Pathway Filter",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 UniProt Pathway & Annotation Filter")
st.markdown(
    "Query UniProt for protein annotations and filter them by pathway keywords. "
    "Powered by [acore](https://github.com/biosustain/acore) and the "
    "[UniProt REST API](https://www.uniprot.org/help/return_fields)."
)
# endregion

# region Sidebar inputs
# ---------------------------------------------------------------------------
# Sidebar – inputs
# ---------------------------------------------------------------------------
with st.sidebar:
    st.header("⚙️ Query settings")

    # UniProt IDs
    ids_text = st.text_area(
        "UniProt accession IDs (one per line)",
        value="P05067\nP12345",
        height=150,
        help="Enter UniProt accession IDs, one per line (e.g. P05067).",
    )

    # Field selection
    st.subheader("Annotation fields")
    st.markdown(
        "Choose which UniProt "
        "[return fields](https://www.uniprot.org/help/return_fields) "
        "to fetch."
    )
    selected_fields = st.multiselect(
        "Fields",
        options=FIELD_DISPLAY_NAMES,
        default=[
            "Gene Ontology (biological process)",
            "Gene Ontology (molecular function)",
            "Gene Ontology (cellular component)",
            "Pathway",
            "Reactome",
        ],
        format_func=lambda name: f"{name}  —  `{UNIPROT_FIELDS[name][0]}`",
        help="Select one or more annotation fields to retrieve from UniProt.",
    )

    fetch_btn = st.button("🔍 Fetch annotations", type="primary", width="stretch")
# endregion

# region Main area filtering
# ---------------------------------------------------------------------------
# Main area – filtering
# ---------------------------------------------------------------------------
col1, col2 = st.columns([3, 2])

with col2:
    st.subheader("🔎 Filter annotations")
    keywords_text = st.text_area(
        "Keywords (one per line)",
        placeholder="e.g.\napoptosis\nmitochondria\ncell cycle",
        height=150,
        help=(
            "Enter keywords to filter annotations. "
            "A row is kept if it contains any of the keywords."
        ),
    )
    case_sensitive = st.checkbox("Case-sensitive matching", value=False)
    filter_btn = st.button("Apply filter", width="stretch")
# endregion

# region Session state helpers
# ---------------------------------------------------------------------------
# Session state helpers
# ---------------------------------------------------------------------------
if "annotations" not in st.session_state:
    st.session_state["annotations"] = None
if "filtered" not in st.session_state:
    st.session_state["filtered"] = None
# endregion

# region Fetch logic
# ---------------------------------------------------------------------------
# Fetch logic
# ---------------------------------------------------------------------------
if fetch_btn:
    uniprot_ids = [line.strip() for line in ids_text.splitlines() if line.strip()]
    if not uniprot_ids:
        st.error("Please enter at least one UniProt accession ID.")
    elif not selected_fields:
        st.error("Please select at least one annotation field.")
    else:
        fields_str = fields_to_api_string(selected_fields)
        with st.spinner(f"Fetching annotations for {len(uniprot_ids)} protein(s)…"):
            try:
                df = query_uniprot(uniprot_ids, fields=fields_str)
                st.session_state["annotations"] = df
                st.session_state["filtered"] = None
                st.success(
                    f"Fetched {len(df):,} annotation rows for "
                    f"{df['identifier'].nunique()} protein(s)."
                )
            except Exception as exc:
                st.error(f"Error fetching data from UniProt: {exc}")
# endregion

# region Filter logic
# ---------------------------------------------------------------------------
# Filter logic
# ---------------------------------------------------------------------------
if filter_btn:
    if st.session_state["annotations"] is None:
        st.warning("Fetch annotations first before applying a filter.")
    else:
        keywords = [kw.strip() for kw in keywords_text.splitlines() if kw.strip()]
        filtered = filter_annotations(
            st.session_state["annotations"],
            keywords=keywords,
            case_sensitive=case_sensitive,
        )
        st.session_state["filtered"] = filtered
        if filtered.empty:
            st.warning("No annotations matched the given keywords.")
        else:
            st.success(f"Filter matched {len(filtered):,} rows.")
# endregion

# region Display results
# ---------------------------------------------------------------------------
# Display results
# ---------------------------------------------------------------------------
with col1:
    annotations = st.session_state["annotations"]
    filtered = st.session_state["filtered"]

    if annotations is not None:
        display_df = filtered if filtered is not None else annotations
        label = "Filtered annotations" if filtered is not None else "All annotations"

        st.subheader(f"📋 {label}")
        st.caption(
            f"{len(display_df):,} rows · "
            f"{display_df['identifier'].nunique() if not display_df.empty else 0}"
            " proteins · "
            f"{display_df['source'].nunique() if not display_df.empty else 0}"
            " annotation sources"
        )

        # Show summary pivot (proteins × sources)
        if not display_df.empty:
            pivot = (
                display_df.groupby(["identifier", "source"])["annotation"]
                .apply(lambda s: " | ".join(s.astype(str)))
                .unstack("source")
                .reset_index()
            )
            st.dataframe(pivot, width="stretch", height=400)
        else:
            st.info("No data to display.")

        # Export section
        st.subheader("💾 Export")
        export_df = display_df if not display_df.empty else annotations

        csv_bytes = export_to_csv(export_df)
        st.download_button(
            label="⬇️ Download as CSV",
            data=csv_bytes,
            file_name="uniprot_annotations.csv",
            mime="text/csv",
            width="stretch",
        )

        # Also offer JSON
        json_bytes = export_df.to_json(orient="records", indent=2).encode("utf-8")
        st.download_button(
            label="⬇️ Download as JSON",
            data=json_bytes,
            file_name="uniprot_annotations.json",
            mime="application/json",
            width="stretch",
        )
    else:
        st.info("Enter UniProt IDs and click **Fetch annotations** to start.")
# endregion
