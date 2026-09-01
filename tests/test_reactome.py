"""Tests for local Reactome pathway evidence integration."""

from pathlib import Path

import pytest

from src.evidence_integration.models import EvidenceType
from src.evidence_integration.reactome import (
    collect_pathway_evidence,
    load_reactome_mapping,
)


def write_mapping(tmp_path: Path, lines: list[str]) -> Path:
    """Create a temporary Reactome mapping file for testing."""
    path = tmp_path / "NCBI2Reactome.txt"
    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
    return path


def test_load_reactome_mapping_retains_human_rows(tmp_path):
    lines = [
        (
            "7157\tR-HSA-123\thttps://reactome.org/123\t"
            "DNA repair\tTAS\tHomo sapiens"
        ),
        (
            "22059\tR-MMU-456\thttps://reactome.org/456\t"
            "Mouse pathway\tTAS\tMus musculus"
        ),
    ]

    path = write_mapping(tmp_path, lines)

    mapping = load_reactome_mapping(path)

    assert "7157" in mapping
    assert "22059" not in mapping

    assert mapping["7157"][0]["pathway_id"] == "R-HSA-123"
    assert mapping["7157"][0]["pathway_name"] == "DNA repair"
    assert mapping["7157"][0]["species"] == "Homo sapiens"


def test_multiple_pathways_for_same_gene_are_preserved(tmp_path):
    lines = [
        (
            "7157\tR-HSA-456\thttps://reactome.org/456\t"
            "Cell cycle\tTAS\tHomo sapiens"
        ),
        (
            "7157\tR-HSA-123\thttps://reactome.org/123\t"
            "DNA repair\tTAS\tHomo sapiens"
        ),
    ]

    path = write_mapping(tmp_path, lines)

    mapping = load_reactome_mapping(path)

    assert len(mapping["7157"]) == 2

    assert [
        item["pathway_id"]
        for item in mapping["7157"]
    ] == [
        "R-HSA-123",
        "R-HSA-456",
    ]


def test_missing_mapping_file_is_rejected(tmp_path):
    missing_path = tmp_path / "missing.txt"

    with pytest.raises(FileNotFoundError):
        load_reactome_mapping(missing_path)


def test_directory_instead_of_file_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="not a file"):
        load_reactome_mapping(tmp_path)


def test_invalid_column_count_is_rejected(tmp_path):
    path = write_mapping(
        tmp_path,
        [
            "7157\tR-HSA-123\ttoo-few-columns",
        ],
    )

    with pytest.raises(
        ValueError,
        match="expected 6 columns",
    ):
        load_reactome_mapping(path)


def test_collect_pathway_evidence_constructs_records():
    mapping = {
        "7157": [
            {
                "pathway_id": "R-HSA-123",
                "pathway_url": "https://reactome.org/123",
                "pathway_name": "DNA repair",
                "evidence_code": "TAS",
                "species": "Homo sapiens",
            },
            {
                "pathway_id": "R-HSA-456",
                "pathway_url": "https://reactome.org/456",
                "pathway_name": "Cell cycle",
                "evidence_code": "TAS",
                "species": "Homo sapiens",
            },
        ]
    }

    records = collect_pathway_evidence(
        "TP53",
        "7157",
        reactome_mapping=mapping,
        retrieved_at="2026-08-26T12:00:00+05:30",
        source_version="97",
    )

    assert len(records) == 2

    assert all(
        record.gene_id == "TP53"
        for record in records
    )

    assert all(
        record.evidence_type is EvidenceType.PATHWAY
        for record in records
    )

    assert all(
        record.source == "Reactome"
        for record in records
    )

    assert all(
        record.source_version == "97"
        for record in records
    )


def test_gene_without_pathways_returns_empty_list():
    records = collect_pathway_evidence(
        "UNKNOWN",
        "999999",
        reactome_mapping={},
        retrieved_at="2026-08-26T12:00:00+05:30",
        source_version="97",
    )

    assert records == []


def test_pathway_records_are_sorted_deterministically():
    mapping = {
        "7157": [
            {
                "pathway_id": "R-HSA-999",
                "pathway_url": "https://reactome.org/999",
                "pathway_name": "Later pathway",
                "evidence_code": "TAS",
                "species": "Homo sapiens",
            },
            {
                "pathway_id": "R-HSA-111",
                "pathway_url": "https://reactome.org/111",
                "pathway_name": "Earlier pathway",
                "evidence_code": "TAS",
                "species": "Homo sapiens",
            },
        ]
    }

    records = collect_pathway_evidence(
        "TP53",
        "7157",
        reactome_mapping=mapping,
        retrieved_at="2026-08-26T12:00:00+05:30",
        source_version="97",
    )

    assert [
        record.evidence_id
        for record in records
    ] == [
        "R-HSA-111",
        "R-HSA-999",
    ]


@pytest.mark.parametrize(
    "gene_symbol",
    [
        "",
        " ",
        "\t",
    ],
)
def test_invalid_gene_symbol_is_rejected(gene_symbol):
    with pytest.raises(
        ValueError,
        match="gene_symbol",
    ):
        collect_pathway_evidence(
            gene_symbol,
            "7157",
            reactome_mapping={},
            retrieved_at="2026-08-26T12:00:00+05:30",
            source_version="97",
        )


@pytest.mark.parametrize(
    "ncbi_gene_id",
    [
        "",
        " ",
        "\t",
    ],
)
def test_invalid_ncbi_gene_id_is_rejected(ncbi_gene_id):
    with pytest.raises(
        ValueError,
        match="ncbi_gene_id",
    ):
        collect_pathway_evidence(
            "TP53",
            ncbi_gene_id,
            reactome_mapping={},
            retrieved_at="2026-08-26T12:00:00+05:30",
            source_version="97",
        )


@pytest.mark.parametrize(
    "source_version",
    [
        "",
        " ",
        "\t",
    ],
)
def test_invalid_source_version_is_rejected(source_version):
    with pytest.raises(
        ValueError,
        match="source_version",
    ):
        collect_pathway_evidence(
            "TP53",
            "7157",
            reactome_mapping={},
            retrieved_at="2026-08-26T12:00:00+05:30",
            source_version=source_version,
        )