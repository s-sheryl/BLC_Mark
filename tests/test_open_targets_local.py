import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from src.evidence_integration.models import EvidenceType
from src.evidence_integration.open_targets_local import (
    BLC_MARK_DISEASE_IDS,
    OpenTargetsAssociation,
    collect_local_cancer_association_evidence,
    load_open_targets_associations,
)


def write_associations(tmp_path):
    path = tmp_path / "associations.parquet"

    table = pa.table(
        {
            "diseaseId": [
                "MONDO_0007254",
                "MONDO_0005061",
            ],
            "targetId": [
                "ENSG00000141510",
                "ENSG00000141510",
            ],
            "associationScore": [
                0.75,
                0.72,
            ],
            "evidenceCount": [
                10,
                8,
            ],
        }
    )

    pq.write_table(
        table,
        path,
    )

    return path


def test_cohort_disease_mapping():
    assert (
        BLC_MARK_DISEASE_IDS["TCGA-BRCA"]
        == "MONDO_0007254"
    )

    assert (
        BLC_MARK_DISEASE_IDS["TCGA-LUAD"]
        == "MONDO_0005061"
    )

    assert (
        BLC_MARK_DISEASE_IDS["TCGA-COAD"]
        == "MONDO_0002271"
    )


def test_load_open_targets_associations(tmp_path):
    path = write_associations(
        tmp_path
    )

    index = load_open_targets_associations(
        path
    )

    association = index[
        (
            "MONDO_0007254",
            "ENSG00000141510",
        )
    ]

    assert isinstance(
        association,
        OpenTargetsAssociation,
    )

    assert association.association_score == 0.75
    assert association.evidence_count == 10


def test_missing_file_is_rejected(tmp_path):
    with pytest.raises(
        FileNotFoundError
    ):
        load_open_targets_associations(
            tmp_path / "missing.parquet"
        )


def test_collect_brca_association():
    index = {
        (
            "MONDO_0007254",
            "ENSG00000141510",
        ): OpenTargetsAssociation(
            disease_id="MONDO_0007254",
            target_id="ENSG00000141510",
            association_score=0.75,
            evidence_count=10,
        )
    }

    records = (
        collect_local_cancer_association_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-BRCA",
            association_index=index,
            retrieved_at=(
                "2026-08-26T16:00:00+05:30"
            ),
            source_version="26.06",
        )
    )

    assert len(records) == 1

    record = records[0]

    assert (
        record.evidence_type
        is EvidenceType.CANCER_ASSOCIATION
    )

    assert record.source == "Open Targets Platform"
    assert record.cancer_cohort == "TCGA-BRCA"
    assert record.evidence_id == "MONDO_0007254"
    assert "0.750000" in record.description
    assert "evidence count=10" in record.description


def test_missing_association_returns_empty():
    records = (
        collect_local_cancer_association_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-BRCA",
            association_index={},
            retrieved_at=(
                "2026-08-26T16:00:00+05:30"
            ),
            source_version="26.06",
        )
    )

    assert records == []


def test_unsupported_cohort_is_rejected():
    with pytest.raises(
        ValueError,
        match="Unsupported",
    ):
        collect_local_cancer_association_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-OV",
            association_index={},
            retrieved_at=(
                "2026-08-26T16:00:00+05:30"
            ),
            source_version="26.06",
        )