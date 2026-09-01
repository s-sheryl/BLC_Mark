import pytest

from src.biomarker_prioritization.models import (
    ComponentScores,
    EvidenceAvailability,
    PrioritizationInput,
    PrioritizedBiomarker,
)
from src.biomarker_prioritization.qc import (
    build_phase5_qc_report,
)


def _biomarker(
    gene_id: str,
    *,
    final_score: float | None,
    availability: EvidenceAvailability,
    clinical_category: str | None = None,
    cohort_count: int = 1,
) -> PrioritizedBiomarker:
    raw_input = PrioritizationInput(
        gene_id=gene_id,
        cancer_cohort="TCGA-BRCA",
        effect_size=1.0,
        effect_size_label="log2_fold_change",
        adjusted_p_value=0.01,
        cancer_association_score=(
            None
            if availability
            == EvidenceAvailability.UNAVAILABLE
            else 0.0
        ),
        clinical_category=clinical_category,
        clinical_direction=None,
        cross_cancer_cohort_count=(
            cohort_count
        ),
        cancer_association_availability=(
            availability
        ),
        clinical_availability=(
            availability
        ),
        cross_cancer_availability=(
            availability
        ),
    )

    scores = ComponentScores(
        de_score=0.5,
        cancer_association_score=(
            None
            if availability
            == EvidenceAvailability.UNAVAILABLE
            else 0.0
        ),
        clinical_score=(
            None
            if availability
            == EvidenceAvailability.UNAVAILABLE
            else 0.0
        ),
        cross_cancer_score=(
            None
            if availability
            == EvidenceAvailability.UNAVAILABLE
            else 0.0
        ),
    )

    return PrioritizedBiomarker(
        gene_id=gene_id,
        cancer_cohort="TCGA-BRCA",
        raw_input=raw_input,
        component_scores=scores,
        final_score=final_score,
        rank=(
            1
            if final_score is not None
            else None
        ),
    )


def test_qc_report_counts_scored_and_unavailable():
    biomarkers = [
        _biomarker(
            "TP53",
            final_score=0.8,
            availability=(
                EvidenceAvailability.AVAILABLE
            ),
            clinical_category=(
                "validated prognostic - unfavorable"
            ),
            cohort_count=3,
        ),
        _biomarker(
            "?|90288",
            final_score=None,
            availability=(
                EvidenceAvailability.UNAVAILABLE
            ),
            cohort_count=1,
        ),
    ]

    report = build_phase5_qc_report(
        biomarkers,
        candidate_count=2,
    )

    assert report[
        "candidate_count"
    ] == 2

    assert report[
        "prioritized_biomarker_count"
    ] == 2

    assert report[
        "scored_candidate_count"
    ] == 1

    assert report[
        "unavailable_score_count"
    ] == 1


def test_qc_report_tracks_score_range():
    biomarkers = [
        _biomarker(
            "A",
            final_score=0.2,
            availability=(
                EvidenceAvailability.AVAILABLE
            ),
        ),
        _biomarker(
            "B",
            final_score=0.9,
            availability=(
                EvidenceAvailability.AVAILABLE
            ),
        ),
    ]

    report = build_phase5_qc_report(
        biomarkers,
        candidate_count=2,
    )

    assert report[
        "minimum_final_score"
    ] == pytest.approx(0.2)

    assert report[
        "maximum_final_score"
    ] == pytest.approx(0.9)


def test_qc_report_tracks_availability():
    biomarkers = [
        _biomarker(
            "A",
            final_score=0.2,
            availability=(
                EvidenceAvailability.NO_SUPPORT
            ),
        ),
        _biomarker(
            "B",
            final_score=None,
            availability=(
                EvidenceAvailability.UNAVAILABLE
            ),
        ),
    ]

    report = build_phase5_qc_report(
        biomarkers,
        candidate_count=2,
    )

    assert report[
        "evidence_availability_counts"
    ]["no_support"] == 3

    assert report[
        "evidence_availability_counts"
    ]["unavailable"] == 3


def test_qc_report_rejects_count_mismatch():
    biomarkers = [
        _biomarker(
            "TP53",
            final_score=0.8,
            availability=(
                EvidenceAvailability.AVAILABLE
            ),
        )
    ]

    with pytest.raises(ValueError):
        build_phase5_qc_report(
            biomarkers,
            candidate_count=2,
        )


def test_qc_report_rejects_negative_candidate_count():
    with pytest.raises(ValueError):
        build_phase5_qc_report(
            [],
            candidate_count=-1,
        )