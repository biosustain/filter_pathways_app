"""Tests for the kegg_pathways module (network-free parsing/orchestration)."""

import pandas as pd

import filter_pathways_app.kegg_pathways as kp
from filter_pathways_app.kegg_pathways import (
    KEGG_PATHWAY_SOURCE,
    format_pathway_annotation,
    gene_map_from_annotations,
    parse_gene_pathway_links,
    parse_kegg_gene_ids,
    parse_pathway_names,
    pathway_organism,
    query_kegg_pathways,
)

# ---------------------------------------------------------------------------
# parse_kegg_gene_ids
# ---------------------------------------------------------------------------


def test_parse_kegg_gene_ids_single():
    assert parse_kegg_gene_ids("hsa:351;") == ["hsa:351"]


def test_parse_kegg_gene_ids_multiple():
    assert parse_kegg_gene_ids("yli:2912002;yli:2910294;") == [
        "yli:2912002",
        "yli:2910294",
    ]


def test_parse_kegg_gene_ids_deduplicates():
    assert parse_kegg_gene_ids("hsa:351;hsa:351;") == ["hsa:351"]


def test_parse_kegg_gene_ids_empty_and_nonstring():
    assert parse_kegg_gene_ids("") == []
    assert parse_kegg_gene_ids(pd.NA) == []
    assert parse_kegg_gene_ids(None) == []


# ---------------------------------------------------------------------------
# parse_gene_pathway_links
# ---------------------------------------------------------------------------


def test_parse_gene_pathway_links():
    raw = "hsa:351\tpath:hsa00010\nhsa:351\tpath:hsa04930\nyli:1\tpath:yli00010\n"
    result = parse_gene_pathway_links(raw)
    assert result == {
        "hsa:351": ["hsa00010", "hsa04930"],
        "yli:1": ["yli00010"],
    }


def test_parse_gene_pathway_links_empty():
    assert parse_gene_pathway_links("") == {}


def test_parse_gene_pathway_links_ignores_blank_lines():
    raw = "\nhsa:351\tpath:hsa00010\n\n"
    assert parse_gene_pathway_links(raw) == {"hsa:351": ["hsa00010"]}


# ---------------------------------------------------------------------------
# parse_pathway_names / helpers
# ---------------------------------------------------------------------------


def test_parse_pathway_names_strips_organism():
    raw = (
        "hsa00010\tGlycolysis / Gluconeogenesis - Homo sapiens (human)\n"
        "hsa04930\tType II diabetes mellitus - Homo sapiens (human)\n"
    )
    assert parse_pathway_names(raw) == {
        "hsa00010": "Glycolysis / Gluconeogenesis",
        "hsa04930": "Type II diabetes mellitus",
    }


def test_pathway_organism():
    assert pathway_organism("hsa00010") == "hsa"
    assert pathway_organism("yli00592") == "yli"
    assert pathway_organism("not-a-pathway") is None


def test_format_pathway_annotation():
    assert (
        format_pathway_annotation("hsa00010", "Glycolysis / Gluconeogenesis")
        == "Glycolysis / Gluconeogenesis [hsa00010]"
    )
    # Falls back to the bare ID when no name is available.
    assert format_pathway_annotation("hsa99999", "") == "hsa99999"


# ---------------------------------------------------------------------------
# gene_map_from_annotations (reuse the KEGG xref already in the annotations DF)
# ---------------------------------------------------------------------------


def _annotations_df():
    return pd.DataFrame(
        {
            "identifier": ["Q05493", "Q05493", "B5FVB1", "P87200"],
            "source": ["KEGG", "Entry", "KEGG", "Pathway"],
            "annotation": ["yli:2912002;", "Q05493", "yli:7009397;", "some pathway"],
        }
    )


def test_gene_map_from_annotations_parses_kegg_rows():
    assert gene_map_from_annotations(_annotations_df()) == {
        "Q05493": ["yli:2912002"],
        "B5FVB1": ["yli:7009397"],
    }


def test_gene_map_from_annotations_matches_fetch_kegg_gene_map(monkeypatch):
    # The DF-derived map must equal what a fresh UniProt fetch would produce.
    df = _annotations_df()

    def fake_fetch_annotations(ids, fields):
        return pd.DataFrame(
            {
                "From": ["Q05493", "B5FVB1"],
                "Entry": ["Q05493", "B5FVB1"],
                "KEGG": ["yli:2912002;", "yli:7009397;"],
            }
        )

    monkeypatch.setattr(kp, "fetch_annotations", fake_fetch_annotations)
    assert gene_map_from_annotations(df) == kp.fetch_kegg_gene_map(
        ["Q05493", "B5FVB1"]
    )


def test_gene_map_from_annotations_no_kegg_rows():
    df = pd.DataFrame(
        {
            "identifier": ["P1"],
            "source": ["Entry"],
            "annotation": ["P1"],
        }
    )
    assert gene_map_from_annotations(df) == {}
    assert gene_map_from_annotations(df.iloc[0:0]) == {}


# ---------------------------------------------------------------------------
# query_kegg_pathways orchestration (network calls monkeypatched)
# ---------------------------------------------------------------------------


def test_query_kegg_pathways_builds_long_format(monkeypatch):
    monkeypatch.setattr(
        kp,
        "fetch_kegg_gene_map",
        lambda ids: {"P1": ["hsa:351"], "P2": ["yli:1"]},
    )
    monkeypatch.setattr(
        kp,
        "fetch_gene_pathways",
        lambda genes: {
            "hsa:351": ["hsa00010", "hsa04930"],
            "yli:1": ["yli00010"],
        },
    )
    monkeypatch.setattr(
        kp,
        "fetch_pathway_names",
        lambda pids: {
            "hsa00010": "Glycolysis / Gluconeogenesis",
            "hsa04930": "Type II diabetes mellitus",
            "yli00010": "Glycolysis / Gluconeogenesis",
        },
    )

    result = query_kegg_pathways(["P1", "P2"])

    assert list(result.columns) == ["identifier", "source", "annotation"]
    assert set(result["source"]) == {KEGG_PATHWAY_SOURCE}
    assert set(result["identifier"]) == {"P1", "P2"}
    assert (
        "Glycolysis / Gluconeogenesis [hsa00010]"
        in result[result["identifier"] == "P1"]["annotation"].tolist()
    )
    assert len(result) == 3


def test_query_kegg_pathways_uses_supplied_gene_map(monkeypatch):
    # A provided gene_map must be used as-is, without calling fetch_kegg_gene_map.
    def _boom(ids):
        raise AssertionError("fetch_kegg_gene_map should not be called")

    monkeypatch.setattr(kp, "fetch_kegg_gene_map", _boom)
    monkeypatch.setattr(
        kp, "fetch_gene_pathways", lambda genes: {"hsa:351": ["hsa00010"]}
    )
    monkeypatch.setattr(
        kp, "fetch_pathway_names", lambda pids: {"hsa00010": "Glycolysis"}
    )

    result = query_kegg_pathways(["P1"], gene_map={"P1": ["hsa:351"]})

    assert len(result) == 1
    assert result.iloc[0]["annotation"] == "Glycolysis [hsa00010]"


def test_query_kegg_pathways_no_kegg_xref(monkeypatch):
    monkeypatch.setattr(kp, "fetch_kegg_gene_map", lambda ids: {})
    result = query_kegg_pathways(["P1"])
    assert result.empty
    assert list(result.columns) == ["identifier", "source", "annotation"]


def test_query_kegg_pathways_no_pathways(monkeypatch):
    monkeypatch.setattr(kp, "fetch_kegg_gene_map", lambda ids: {"P1": ["hsa:351"]})
    monkeypatch.setattr(kp, "fetch_gene_pathways", lambda genes: {})
    result = query_kegg_pathways(["P1"])
    assert result.empty


def test_query_kegg_pathways_deduplicates_across_genes(monkeypatch):
    # A protein with two KEGG genes that share a pathway -> single row.
    monkeypatch.setattr(
        kp, "fetch_kegg_gene_map", lambda ids: {"P1": ["hsa:1", "hsa:2"]}
    )
    monkeypatch.setattr(
        kp,
        "fetch_gene_pathways",
        lambda genes: {"hsa:1": ["hsa00010"], "hsa:2": ["hsa00010"]},
    )
    monkeypatch.setattr(
        kp, "fetch_pathway_names", lambda pids: {"hsa00010": "Glycolysis"}
    )
    result = query_kegg_pathways(["P1"])
    assert len(result) == 1
    assert result.iloc[0]["annotation"] == "Glycolysis [hsa00010]"
