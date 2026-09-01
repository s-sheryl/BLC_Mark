import gzip

import pytest

from src.evidence_integration.models import EvidenceType
from src.evidence_integration.ncbi_local import (
    NCBIGeneRecord,
    collect_local_functional_evidence,
    load_human_gene_info,
)


def write_gene_info(
    tmp_path,
    lines,
):
    path = tmp_path / "Homo_sapiens.gene_info.gz"

    header = (
        "#tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\t"
        "dbXrefs\tchromosome\tmap_location\tdescription\n"
    )

    with gzip.open(
        path,
        "wt",
        encoding="utf-8",
    ) as handle:
        handle.write(header)

        for line in lines:
            handle.write(
                line + "\n"
            )

    return path


def test_load_human_gene_info(tmp_path):
    path = write_gene_info(
        tmp_path,
        [
            (
                "9606\t7157\tTP53\t-\tP53|BCC7\t"
                "MIM:191170|Ensembl:ENSG00000141510\t"
                "17\t17p13.1\ttumor protein p53"
            )
        ],
    )

    index = load_human_gene_info(
        path
    )

    assert "TP53" in index

    record = index["TP53"]

    assert record.gene_id == "7157"
    assert record.ensembl_id == "ENSG00000141510"
    assert record.description == "tumor protein p53"
    assert record.aliases == (
        "BCC7",
        "P53",
    )


def test_gene_without_ensembl_mapping_is_allowed(tmp_path):
    path = write_gene_info(
        tmp_path,
        [
            (
                "9606\t123\tGENE1\t-\t-\t"
                "MIM:123\t1\t1p1\tExample gene"
            )
        ],
    )

    index = load_human_gene_info(
        path
    )

    assert index["GENE1"].ensembl_id is None


def test_missing_file_is_rejected(tmp_path):
    with pytest.raises(
        FileNotFoundError
    ):
        load_human_gene_info(
            tmp_path / "missing.gz"
        )


def test_local_functional_evidence():
    index = {
        "TP53": NCBIGeneRecord(
            gene_id="7157",
            symbol="TP53",
            ensembl_id="ENSG00000141510",
            description="tumor protein p53",
            aliases=("P53",),
        )
    }

    evidence = (
        collect_local_functional_evidence(
            "TP53",
            ncbi_index=index,
            retrieved_at=(
                "2026-08-26T15:00:00+05:30"
            ),
            source_version="2026-08-05",
        )
    )

    assert evidence is not None
    assert evidence.evidence_type is EvidenceType.FUNCTIONAL
    assert evidence.source == "NCBI Gene"
    assert evidence.source_version == "2026-08-05"
    assert evidence.evidence_id == "NCBI_GENE:7157"


def test_unknown_gene_returns_none():
    assert (
        collect_local_functional_evidence(
            "UNKNOWN",
            ncbi_index={},
            retrieved_at=(
                "2026-08-26T15:00:00+05:30"
            ),
            source_version="2026-08-05",
        )
        is None
    )


def test_gene_without_description_returns_none():
    index = {
        "GENE1": NCBIGeneRecord(
            gene_id="123",
            symbol="GENE1",
            ensembl_id=None,
            description="",
            aliases=(),
        )
    }

    assert (
        collect_local_functional_evidence(
            "GENE1",
            ncbi_index=index,
            retrieved_at=(
                "2026-08-26T15:00:00+05:30"
            ),
            source_version="2026-08-05",
        )
        is None
    )

def test_duplicate_symbols_are_treated_as_ambiguous(tmp_path):
    path = write_gene_info(
        tmp_path,
        [
            (
                "9606\t111\tRNR1\t-\t-\t"
                "Ensembl:ENSG00000111111\t"
                "1\t1p1\tFirst RNR1 record"
            ),
            (
                "9606\t222\tRNR1\t-\t-\t"
                "Ensembl:ENSG00000222222\t"
                "2\t2p2\tSecond RNR1 record"
            ),
        ],
    )

    index = load_human_gene_info(
        path
    )

    assert "RNR1" not in index