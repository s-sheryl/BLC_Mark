import pytest

from src.biomarker_prioritization.models import (
    ComponentScores,
    EvidenceAvailability,
    PrioritizationInput,
    PrioritizedBiomarker,
)


def _valid_input() -> PrioritizationInput:
    return PrioritizationInput(
        gene_id="TP53",
        cancer_cohort="TCGA-BRCA",
        effect_size=2.5,
        effect_size_label="log2_fold_change",
        adjusted_p_value=0.001,
        cancer_association_score=0.8,
        clinical_category="validated prognostic",
        clinical_direction="unfavorable",
        cross_cancer_cohort_count=3,
        cancer_association_availability=EvidenceAvailability.AVAILABLE,
        clinical_availability=EvidenceAvailability.AVAILABLE,
        cross_cancer_availability=EvidenceAvailability.AVAILABLE,
        functional_description="tumor protein p53",
        pathway_count=12,
    )


def test_prioritization_input_accepts_valid_data():
    item = _valid_input()

    assert item.gene_id == "TP53"
    assert item.cancer_cohort == "TCGA-BRCA"
    assert item.cross_cancer_cohort_count == 3


def test_prioritization_input_rejects_empty_gene_id():
    with pytest.raises(ValueError):
        PrioritizationInput(
            gene_id="",
            cancer_cohort="TCGA-BRCA",
            effect_size=1.0,
            effect_size_label="log2_fold_change",
            adjusted_p_value=0.01,
            cancer_association_score=None,
            clinical_category=None,
            clinical_direction=None,
            cross_cancer_cohort_count=1,
            cancer_association_availability=EvidenceAvailability.NO_SUPPORT,
            clinical_availability=EvidenceAvailability.NO_SUPPORT,
            cross_cancer_availability=EvidenceAvailability.NO_SUPPORT,
        )


def test_prioritization_input_rejects_bad_adjusted_p_value():
    data = _valid_input().__dict__.copy()
    data["adjusted_p_value"] = 1.5

    with pytest.raises(ValueError):
        PrioritizationInput(**data)


def test_prioritization_input_rejects_invalid_cancer_score():
    data = _valid_input().__dict__.copy()
    data["cancer_association_score"] = 1.2

    with pytest.raises(ValueError):
        PrioritizationInput(**data)


def test_prioritization_input_rejects_cross_cancer_count_below_one():
    data = _valid_input().__dict__.copy()
    data["cross_cancer_cohort_count"] = 0

    with pytest.raises(ValueError):
        PrioritizationInput(**data)


def test_component_scores_accept_none_for_unavailable_component():
    scores = ComponentScores(
        de_score=0.9,
        cancer_association_score=None,
        clinical_score=None,
        cross_cancer_score=0.5,
    )

    assert scores.cancer_association_score is None


def test_component_scores_reject_values_above_one():
    with pytest.raises(ValueError):
        ComponentScores(
            de_score=1.1,
            cancer_association_score=0.5,
            clinical_score=0.5,
            cross_cancer_score=0.5,
        )


def test_prioritized_biomarker_accepts_valid_result():
    raw_input = _valid_input()

    scores = ComponentScores(
        de_score=0.9,
        cancer_association_score=0.8,
        clinical_score=1.0,
        cross_cancer_score=1.0,
    )

    result = PrioritizedBiomarker(
        gene_id="TP53",
        cancer_cohort="TCGA-BRCA",
        raw_input=raw_input,
        component_scores=scores,
        final_score=0.925,
        rank=1,
    )

    assert result.final_score == pytest.approx(0.925)
    assert result.rank == 1


def test_prioritized_biomarker_rejects_invalid_final_score():
    raw_input = _valid_input()

    scores = ComponentScores(
        de_score=0.9,
        cancer_association_score=0.8,
        clinical_score=1.0,
        cross_cancer_score=1.0,
    )

    with pytest.raises(ValueError):
        PrioritizedBiomarker(
            gene_id="TP53",
            cancer_cohort="TCGA-BRCA",
            raw_input=raw_input,
            component_scores=scores,
            final_score=1.5,
            rank=1,
        )


def test_prioritized_biomarker_rejects_rank_zero():
    raw_input = _valid_input()

    scores = ComponentScores(
        de_score=0.9,
        cancer_association_score=0.8,
        clinical_score=1.0,
        cross_cancer_score=1.0,
    )

    with pytest.raises(ValueError):
        PrioritizedBiomarker(
            gene_id="TP53",
            cancer_cohort="TCGA-BRCA",
            raw_input=raw_input,
            component_scores=scores,
            final_score=0.9,
            rank=0,
        )