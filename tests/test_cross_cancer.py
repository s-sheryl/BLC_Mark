import pytest

from src.evidence_integration.cross_cancer import (
    build_cross_cancer_index,
    collect_cross_cancer_evidence,
)
from src.evidence_integration.models import EvidenceType


def test_build_cross_cancer_index():
    candidates = {
        "TCGA-BRCA": {"TP53", "BRCA1"},
        "TCGA-LUAD": {"TP53", "EGFR"},
        "TCGA-COAD": {"TP53"},
    }

    index = build_cross_cancer_index(candidates)

    assert index["TP53"] == (
        "TCGA-BRCA",
        "TCGA-COAD",
        "TCGA-LUAD",
    )

    assert index["BRCA1"] == ("TCGA-BRCA",)
    assert index["EGFR"] == ("TCGA-LUAD",)


def test_cross_cancer_record_created_for_shared_gene():
    index = {
        "TP53": (
            "TCGA-BRCA",
            "TCGA-COAD",
            "TCGA-LUAD",
        )
    }

    records = collect_cross_cancer_evidence(
        "TP53",
        "TCGA-BRCA",
        cross_cancer_index=index,
        retrieved_at="2026-08-26T12:00:00+05:30",
    )

    assert len(records) == 1

    record = records[0]

    assert record.evidence_type is EvidenceType.CROSS_CANCER
    assert record.source == "BLC Mark Phase 3"
    assert record.cancer_cohort == "TCGA-BRCA"
    assert "3 cohorts" in record.description


def test_single_cohort_gene_has_no_cross_cancer_evidence():
    index = {
        "BRCA1": ("TCGA-BRCA",),
    }

    records = collect_cross_cancer_evidence(
        "BRCA1",
        "TCGA-BRCA",
        cross_cancer_index=index,
        retrieved_at="2026-08-26T12:00:00+05:30",
    )

    assert records == []


def test_gene_not_present_in_requested_cohort_returns_empty():
    index = {
        "TP53": (
            "TCGA-BRCA",
            "TCGA-LUAD",
        ),
    }

    records = collect_cross_cancer_evidence(
        "TP53",
        "TCGA-COAD",
        cross_cancer_index=index,
        retrieved_at="2026-08-26T12:00:00+05:30",
    )

    assert records == []


def test_index_is_deterministic():
    candidates = {
        "TCGA-LUAD": {"TP53"},
        "TCGA-BRCA": {"TP53"},
        "TCGA-COAD": {"TP53"},
    }

    index = build_cross_cancer_index(candidates)

    assert list(index) == ["TP53"]

    assert index["TP53"] == (
        "TCGA-BRCA",
        "TCGA-COAD",
        "TCGA-LUAD",
    )


def test_invalid_candidate_mapping_type_is_rejected():
    with pytest.raises(
        TypeError,
        match="mapping",
    ):
        build_cross_cancer_index([])


def test_candidate_collection_must_be_set():
    with pytest.raises(
        TypeError,
        match="set",
    ):
        build_cross_cancer_index(
            {
                "TCGA-BRCA": ["TP53"],
            }
        )


@pytest.mark.parametrize(
    "gene_symbol",
    ["", " ", "\t"],
)
def test_blank_gene_symbol_is_rejected(gene_symbol):
    with pytest.raises(
        ValueError,
        match="gene_symbol",
    ):
        collect_cross_cancer_evidence(
            gene_symbol,
            "TCGA-BRCA",
            cross_cancer_index={},
            retrieved_at="2026-08-26T12:00:00+05:30",
        )


@pytest.mark.parametrize(
    "cancer_cohort",
    ["", " ", "\t"],
)
def test_blank_cancer_cohort_is_rejected(cancer_cohort):
    with pytest.raises(
        ValueError,
        match="cancer_cohort",
    ):
        collect_cross_cancer_evidence(
            "TP53",
            cancer_cohort,
            cross_cancer_index={},
            retrieved_at="2026-08-26T12:00:00+05:30",
        )