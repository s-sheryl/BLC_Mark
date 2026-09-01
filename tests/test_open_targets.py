from unittest.mock import patch

import pytest

from src.evidence_integration.models import EvidenceType
from src.evidence_integration.open_targets import (
    collect_cancer_association_evidence,
    fetch_target_disease_associations,
)


def test_fetch_target_returns_target_record():
    payload = {
        "data": {
            "target": {
                "id": "ENSG00000141510",
                "approvedSymbol": "TP53",
                "associatedDiseases": {
                    "count": 1,
                    "rows": [],
                },
            }
        }
    }

    with patch(
        "src.evidence_integration.open_targets._post_graphql",
        return_value=payload,
    ):
        target = fetch_target_disease_associations(
            "ENSG00000141510"
        )

    assert target is not None
    assert target["approvedSymbol"] == "TP53"


def test_missing_target_returns_none():
    payload = {
        "data": {
            "target": None,
        }
    }

    with patch(
        "src.evidence_integration.open_targets._post_graphql",
        return_value=payload,
    ):
        target = fetch_target_disease_associations(
            "ENSG00000000000"
        )

    assert target is None


def test_collect_cancer_association_constructs_relevant_records():
    target = {
        "id": "ENSG00000141510",
        "approvedSymbol": "TP53",
        "associatedDiseases": {
            "count": 3,
            "rows": [
                {
                    "score": 0.75,
                    "disease": {
                        "id": "MONDO_0005575",
                        "name": "colorectal cancer",
                    },
                },
                {
                    "score": 0.72,
                    "disease": {
                        "id": "MONDO_0005061",
                        "name": "lung adenocarcinoma",
                    },
                },
                {
                    "score": 0.90,
                    "disease": {
                        "id": "MONDO_9999999",
                        "name": "unrelated metabolic disorder",
                    },
                },
            ],
        },
    }

    with patch(
        "src.evidence_integration.open_targets."
        "fetch_target_disease_associations",
        return_value=target,
    ):
        records = collect_cancer_association_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-COAD",
            retrieved_at="2026-08-26T12:00:00+05:30",
            source_version="26.06",
        )

    assert len(records) == 1

    record = records[0]

    assert record.evidence_type is EvidenceType.CANCER_ASSOCIATION
    assert record.source == "Open Targets Platform"
    assert record.source_version == "26.06"
    assert record.evidence_id == "MONDO_0005575"
    assert record.cancer_cohort == "TCGA-COAD"
    assert "colorectal cancer" in record.description


def test_lung_association_is_retained_for_luad():
    target = {
        "id": "ENSG00000141510",
        "approvedSymbol": "TP53",
        "associatedDiseases": {
            "count": 2,
            "rows": [
                {
                    "score": 0.72,
                    "disease": {
                        "id": "MONDO_0005061",
                        "name": "lung adenocarcinoma",
                    },
                },
                {
                    "score": 0.75,
                    "disease": {
                        "id": "MONDO_0005575",
                        "name": "colorectal cancer",
                    },
                },
            ],
        },
    }

    with patch(
        "src.evidence_integration.open_targets."
        "fetch_target_disease_associations",
        return_value=target,
    ):
        records = collect_cancer_association_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-LUAD",
            retrieved_at="2026-08-26T12:00:00+05:30",
            source_version="26.06",
        )

    assert len(records) == 1
    assert records[0].evidence_id == "MONDO_0005061"


def test_breast_association_is_retained_for_brca():
    target = {
        "id": "ENSG00000141510",
        "approvedSymbol": "TP53",
        "associatedDiseases": {
            "count": 2,
            "rows": [
                {
                    "score": 0.70,
                    "disease": {
                        "id": "MONDO_0016419",
                        "name": "hereditary breast carcinoma",
                    },
                },
                {
                    "score": 0.72,
                    "disease": {
                        "id": "MONDO_0005061",
                        "name": "lung adenocarcinoma",
                    },
                },
            ],
        },
    }

    with patch(
        "src.evidence_integration.open_targets."
        "fetch_target_disease_associations",
        return_value=target,
    ):
        records = collect_cancer_association_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-BRCA",
            retrieved_at="2026-08-26T12:00:00+05:30",
            source_version="26.06",
        )

    assert len(records) == 1
    assert records[0].evidence_id == "MONDO_0016419"


def test_symbol_mismatch_returns_no_evidence():
    target = {
        "id": "ENSG00000141510",
        "approvedSymbol": "OTHER",
        "associatedDiseases": {
            "count": 0,
            "rows": [],
        },
    }

    with patch(
        "src.evidence_integration.open_targets."
        "fetch_target_disease_associations",
        return_value=target,
    ):
        records = collect_cancer_association_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-BRCA",
            retrieved_at="2026-08-26T12:00:00+05:30",
            source_version="26.06",
        )

    assert records == []


def test_missing_disease_fields_are_skipped():
    target = {
        "id": "ENSG00000141510",
        "approvedSymbol": "TP53",
        "associatedDiseases": {
            "count": 1,
            "rows": [
                {
                    "score": 0.9,
                    "disease": {
                        "id": "",
                        "name": "",
                    },
                },
            ],
        },
    }

    with patch(
        "src.evidence_integration.open_targets."
        "fetch_target_disease_associations",
        return_value=target,
    ):
        records = collect_cancer_association_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-BRCA",
            retrieved_at="2026-08-26T12:00:00+05:30",
            source_version="26.06",
        )

    assert records == []


def test_non_cancer_disease_is_excluded():
    target = {
        "id": "ENSG00000141510",
        "approvedSymbol": "TP53",
        "associatedDiseases": {
            "count": 1,
            "rows": [
                {
                    "score": 0.95,
                    "disease": {
                        "id": "MONDO_123",
                        "name": "cardiomyopathy",
                    },
                },
            ],
        },
    }

    with patch(
        "src.evidence_integration.open_targets."
        "fetch_target_disease_associations",
        return_value=target,
    ):
        records = collect_cancer_association_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-BRCA",
            retrieved_at="2026-08-26T12:00:00+05:30",
            source_version="26.06",
        )

    assert records == []


@pytest.mark.parametrize(
    "ensembl_id",
    ["", " ", "\t"],
)
def test_invalid_ensembl_id_is_rejected(ensembl_id):
    with pytest.raises(ValueError, match="ensembl_id"):
        fetch_target_disease_associations(
            ensembl_id
        )


def test_invalid_page_size_is_rejected():
    with pytest.raises(ValueError, match="page_size"):
        fetch_target_disease_associations(
            "ENSG00000141510",
            page_size=0,
        )


def test_unsupported_cancer_cohort_is_rejected():
    with pytest.raises(
        ValueError,
        match="Unsupported BLC Mark cancer cohort",
    ):
        collect_cancer_association_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-OV",
            retrieved_at="2026-08-26T12:00:00+05:30",
            source_version="26.06",
        )


def test_records_are_deterministically_sorted():
    target = {
        "id": "ENSG00000141510",
        "approvedSymbol": "TP53",
        "associatedDiseases": {
            "count": 2,
            "rows": [
                {
                    "score": 0.9,
                    "disease": {
                        "id": "MONDO_999",
                        "name": "colon cancer subtype",
                    },
                },
                {
                    "score": 0.8,
                    "disease": {
                        "id": "MONDO_111",
                        "name": "colorectal cancer",
                    },
                },
            ],
        },
    }

    with patch(
        "src.evidence_integration.open_targets."
        "fetch_target_disease_associations",
        return_value=target,
    ):
        records = collect_cancer_association_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-COAD",
            retrieved_at="2026-08-26T12:00:00+05:30",
            source_version="26.06",
        )

    assert [
        record.evidence_id
        for record in records
    ] == [
        "MONDO_111",
        "MONDO_999",
    ]