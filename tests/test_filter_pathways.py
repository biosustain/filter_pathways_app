"""Tests for the filter_pathways and uniprot_fields modules."""

import pandas as pd
import pytest

from filter_pathways_app.filter_pathways import export_to_csv, filter_annotations
from filter_pathways_app.uniprot_fields import (
    FIELD_DISPLAY_NAMES,
    UNIPROT_FIELDS,
    fields_to_api_string,
)

# ---------------------------------------------------------------------------
# uniprot_fields tests
# ---------------------------------------------------------------------------


def test_fields_to_api_string_empty():
    result = fields_to_api_string([])
    assert result == "accession"


def test_fields_to_api_string_single():
    result = fields_to_api_string(["GO (biological process)"])
    assert result == "accession,go_p"


def test_fields_to_api_string_multiple():
    result = fields_to_api_string(["GO (biological process)", "Reactome"])
    assert result == "accession,go_p,xref_reactome"


def test_all_field_display_names_in_dict():
    for name in FIELD_DISPLAY_NAMES:
        assert name in UNIPROT_FIELDS, f"{name!r} missing from UNIPROT_FIELDS"


def test_uniprot_fields_structure():
    for name, (api_id, description) in UNIPROT_FIELDS.items():
        assert isinstance(api_id, str) and api_id, f"Invalid api_id for {name!r}"
        assert (
            isinstance(description, str) and description
        ), f"Invalid description for {name!r}"


def test_fields_to_api_string_unknown_field():
    with pytest.raises(KeyError):
        fields_to_api_string(["nonexistent_field"])


# ---------------------------------------------------------------------------
# filter_annotations tests
# ---------------------------------------------------------------------------


@pytest.fixture()
def sample_annotations() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "identifier": ["P1", "P1", "P1", "P2", "P2"],
            "source": ["go_p", "go_p", "go_c", "go_p", "cc_pathway"],
            "annotation": [
                "cell apoptosis [GO:0042981]",
                "cell division [GO:0051301]",
                "mitochondria membrane [GO:0031966]",
                "DNA repair [GO:0006281]",
                "Glycolysis / Gluconeogenesis",
            ],
        }
    )


def test_filter_annotations_single_keyword(sample_annotations):
    result = filter_annotations(sample_annotations, keywords=["apoptosis"])
    assert len(result) == 1
    assert result.iloc[0]["annotation"] == "cell apoptosis [GO:0042981]"


def test_filter_annotations_multiple_keywords(sample_annotations):
    result = filter_annotations(
        sample_annotations, keywords=["apoptosis", "mitochondr"]
    )
    assert len(result) == 2


def test_filter_annotations_case_insensitive(sample_annotations):
    result = filter_annotations(
        sample_annotations, keywords=["APOPTOSIS"], case_sensitive=False
    )
    assert len(result) == 1


def test_filter_annotations_case_sensitive_no_match(sample_annotations):
    result = filter_annotations(
        sample_annotations, keywords=["APOPTOSIS"], case_sensitive=True
    )
    assert len(result) == 0


def test_filter_annotations_empty_keywords(sample_annotations):
    result = filter_annotations(sample_annotations, keywords=[])
    assert len(result) == 0
    assert list(result.columns) == list(sample_annotations.columns)


def test_filter_annotations_whitespace_only_keywords(sample_annotations):
    result = filter_annotations(sample_annotations, keywords=["  ", "\t"])
    assert len(result) == 0


def test_filter_annotations_no_match(sample_annotations):
    result = filter_annotations(sample_annotations, keywords=["spaghetti"])
    assert len(result) == 0


def test_filter_annotations_all_proteins(sample_annotations):
    result = filter_annotations(
        sample_annotations, keywords=["cell", "DNA", "Glycolysis"]
    )
    assert set(result["identifier"]) == {"P1", "P2"}


# ---------------------------------------------------------------------------
# export_to_csv tests
# ---------------------------------------------------------------------------


def test_export_to_csv_returns_bytes():
    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    data = export_to_csv(df)
    assert isinstance(data, bytes)


def test_export_to_csv_bom():
    df = pd.DataFrame({"a": [1]})
    data = export_to_csv(df)
    assert data.startswith(b"\xef\xbb\xbf"), "CSV should start with UTF-8 BOM"


def test_export_to_csv_roundtrip():
    df = pd.DataFrame({"identifier": ["P1"], "annotation": ["apoptosis"]})
    data = export_to_csv(df)
    import io

    result = pd.read_csv(io.BytesIO(data), encoding="utf-8-sig")
    assert list(result.columns) == ["identifier", "annotation"]
    assert result.iloc[0]["identifier"] == "P1"
