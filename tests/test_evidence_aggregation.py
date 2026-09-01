import pytest

from src.evidence_integration.aggregation import (
    GeneEvidenceProfile,
    aggregate_gene_evidence,
)
from src.evidence_integration.models import (
    EvidenceRecord,
    EvidenceType,
)


def make_record(
    evidence_type: EvidenceType,
    evidence_id: str,
    description: str,
) -> EvidenceRecord:
    return EvidenceRecord(
        gene_id="TP53",
        evidence_type=evidence_type,
        source="ExampleSource",
        evidence_id=evidence_id,
        description=description,
    )


def test_gene_evidence_profile_constructs():
    record = make_record(
        EvidenceType.FUNCTIONAL,
        "FUNC:1",
        "Functional evidence.",
    )

    profile = GeneEvidenceProfile(
        gene_id="TP53",
        cancer_cohort="TCGA-BRCA",
        evidence_records=(record,),
    )

    assert profile.gene_id == "TP53"
    assert profile.cancer_cohort == "TCGA-BRCA"
    assert profile.evidence_records == (record,)


def test_aggregate_combines_multiple_evidence_types():
    records = [
        make_record(
            EvidenceType.FUNCTIONAL,
            "NCBI:7157",
            "Functional evidence.",
        ),
        make_record(
            EvidenceType.PATHWAY,
            "R-HSA-123",
            "Pathway evidence.",
        ),
        make_record(
            EvidenceType.CANCER_ASSOCIATION,
            "MONDO_123",
            "Cancer association evidence.",
        ),
    ]

    profile = aggregate_gene_evidence(
        "TP53",
        "TCGA-BRCA",
        records,
    )

    assert len(profile.evidence_records) == 3

    assert {
        record.evidence_type
        for record in profile.evidence_records
    } == {
        EvidenceType.FUNCTIONAL,
        EvidenceType.PATHWAY,
        EvidenceType.CANCER_ASSOCIATION,
    }


def test_duplicate_evidence_records_are_removed():
    record = make_record(
        EvidenceType.PATHWAY,
        "R-HSA-123",
        "DNA repair",
    )

    profile = aggregate_gene_evidence(
        "TP53",
        "TCGA-BRCA",
        [record, record],
    )

    assert len(profile.evidence_records) == 1


def test_output_order_is_deterministic():
    records = [
        make_record(
            EvidenceType.PATHWAY,
            "R-HSA-999",
            "Later",
        ),
        make_record(
            EvidenceType.FUNCTIONAL,
            "NCBI:7157",
            "Functional",
        ),
        make_record(
            EvidenceType.PATHWAY,
            "R-HSA-111",
            "Earlier",
        ),
    ]

    profile = aggregate_gene_evidence(
        "TP53",
        "TCGA-BRCA",
        records,
    )

    assert [
        record.evidence_id
        for record in profile.evidence_records
    ] == [
        "NCBI:7157",
        "R-HSA-111",
        "R-HSA-999",
    ]


def test_empty_evidence_is_allowed():
    profile = aggregate_gene_evidence(
        "TP53",
        "TCGA-BRCA",
        [],
    )

    assert profile.evidence_records == ()


def test_mismatched_gene_record_is_rejected():
    record = EvidenceRecord(
        gene_id="BRCA1",
        evidence_type=EvidenceType.FUNCTIONAL,
        source="ExampleSource",
        evidence_id="EXAMPLE:1",
        description="Evidence",
    )

    with pytest.raises(
        ValueError,
        match="another gene",
    ):
        aggregate_gene_evidence(
            "TP53",
            "TCGA-BRCA",
            [record],
        )


def test_profile_rejects_mismatched_record_gene():
    record = EvidenceRecord(
        gene_id="BRCA1",
        evidence_type=EvidenceType.FUNCTIONAL,
        source="ExampleSource",
        description="Evidence",
    )

    with pytest.raises(
        ValueError,
        match="does not match",
    ):
        GeneEvidenceProfile(
            gene_id="TP53",
            cancer_cohort="TCGA-BRCA",
            evidence_records=(record,),
        )


def test_profile_requires_tuple_records():
    with pytest.raises(
        TypeError,
        match="tuple",
    ):
        GeneEvidenceProfile(
            gene_id="TP53",
            cancer_cohort="TCGA-BRCA",
            evidence_records=[],
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
        GeneEvidenceProfile(
            gene_id="TP53",
            cancer_cohort=cancer_cohort,
            evidence_records=(),
        )


def test_profile_is_immutable():
    profile = GeneEvidenceProfile(
        gene_id="TP53",
        cancer_cohort="TCGA-BRCA",
        evidence_records=(),
    )

    with pytest.raises(Exception):
        profile.gene_id = "BRCA1"