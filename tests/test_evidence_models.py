import pytest

from src.evidence_integration.models import EvidenceRecord, EvidenceType


def test_evidence_type_values_are_stable():
    assert EvidenceType.FUNCTIONAL.value == "functional"
    assert EvidenceType.PATHWAY.value == "pathway"
    assert EvidenceType.CANCER_ASSOCIATION.value == "cancer_association"
    assert EvidenceType.CLINICAL.value == "clinical"
    assert EvidenceType.CROSS_CANCER.value == "cross_cancer"


def test_valid_evidence_record_constructs():
    record = EvidenceRecord(
        gene_id="TP53",
        evidence_type=EvidenceType.CANCER_ASSOCIATION,
        source="ExampleSource",
        description="Gene has documented cancer-associated evidence.",
        source_version="1.0",
        evidence_id="EXAMPLE:123",
        cancer_cohort="TCGA-BRCA",
        retrieved_at="2026-08-26T12:00:00+05:30",
        source_url="https://example.org/123",
    )

    assert record.gene_id == "TP53"
    assert record.evidence_type is EvidenceType.CANCER_ASSOCIATION
    assert record.source == "ExampleSource"
    assert record.evidence_id == "EXAMPLE:123"


def test_whitespace_is_removed_from_required_strings():
    record = EvidenceRecord(
        gene_id=" TP53 ",
        evidence_type=EvidenceType.FUNCTIONAL,
        source=" ExampleSource ",
        description=" Example description. ",
    )

    assert record.gene_id == "TP53"
    assert record.source == "ExampleSource"
    assert record.description == "Example description."


@pytest.mark.parametrize(
    "gene_id",
    ["", " ", "\t"],
)
def test_blank_gene_id_is_rejected(gene_id):
    with pytest.raises(ValueError, match="gene_id"):
        EvidenceRecord(
            gene_id=gene_id,
            evidence_type=EvidenceType.FUNCTIONAL,
            source="ExampleSource",
            description="Evidence",
        )


def test_blank_source_is_rejected():
    with pytest.raises(ValueError, match="source"):
        EvidenceRecord(
            gene_id="TP53",
            evidence_type=EvidenceType.FUNCTIONAL,
            source=" ",
            description="Evidence",
        )


def test_blank_description_is_rejected():
    with pytest.raises(ValueError, match="description"):
        EvidenceRecord(
            gene_id="TP53",
            evidence_type=EvidenceType.FUNCTIONAL,
            source="ExampleSource",
            description=" ",
        )


def test_invalid_evidence_type_is_rejected():
    with pytest.raises(TypeError, match="EvidenceType"):
        EvidenceRecord(
            gene_id="TP53",
            evidence_type="functional",
            source="ExampleSource",
            description="Evidence",
        )


def test_evidence_record_is_immutable():
    record = EvidenceRecord(
        gene_id="TP53",
        evidence_type=EvidenceType.FUNCTIONAL,
        source="ExampleSource",
        description="Evidence",
    )

    with pytest.raises(Exception):
        record.gene_id = "BRCA1"