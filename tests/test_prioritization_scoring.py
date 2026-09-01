import pytest

from src.biomarker_prioritization.models import ComponentScores
from src.biomarker_prioritization.scoring import (
    calculate_de_percentile_scores,
    calculate_final_score,
    score_cancer_association,
    score_clinical_category,
    score_cross_cancer_recurrence,
)


def test_de_percentiles_preserve_original_order():
    scores = calculate_de_percentile_scores(
        [2.0, -1.0, 4.0]
    )

    assert scores == pytest.approx(
        (0.5, 0.0, 1.0)
    )


def test_de_percentiles_use_absolute_effect_size():
    scores = calculate_de_percentile_scores(
        [-3.0, 1.0, 2.0]
    )

    assert scores == pytest.approx(
        (1.0, 0.0, 0.5)
    )


def test_de_percentiles_average_ties():
    scores = calculate_de_percentile_scores(
        [1.0, -1.0, 3.0]
    )

    assert scores == pytest.approx(
        (0.25, 0.25, 1.0)
    )


def test_single_de_candidate_receives_one():
    assert calculate_de_percentile_scores(
        [2.0]
    ) == pytest.approx((1.0,))


def test_de_percentiles_reject_empty_input():
    with pytest.raises(ValueError):
        calculate_de_percentile_scores([])


def test_cancer_association_preserves_valid_score():
    assert score_cancer_association(0.42) == pytest.approx(
        0.42
    )


def test_cancer_association_preserves_unavailable():
    assert score_cancer_association(None) is None


def test_cancer_association_rejects_out_of_range():
    with pytest.raises(ValueError):
        score_cancer_association(1.1)


@pytest.mark.parametrize(
    ("category", "expected"),
    [
        ("unprognostic - favorable", 0.0),
        ("unprognostic - unfavorable", 0.0),
        ("potential prognostic - favorable", 0.5),
        ("potential prognostic - unfavorable", 0.5),
        ("validated prognostic - favorable", 1.0),
        ("validated prognostic - unfavorable", 1.0),
    ],
)
def test_clinical_category_mapping(category, expected):
    assert score_clinical_category(category) == pytest.approx(
        expected
    )


def test_clinical_category_preserves_unavailable():
    assert score_clinical_category(None) is None


def test_clinical_category_rejects_unknown_category():
    with pytest.raises(ValueError):
        score_clinical_category("unknown")


@pytest.mark.parametrize(
    ("count", "expected"),
    [
        (1, 0.0),
        (2, 0.5),
        (3, 1.0),
    ],
)
def test_cross_cancer_mapping(count, expected):
    assert score_cross_cancer_recurrence(count) == pytest.approx(
        expected
    )


def test_cross_cancer_rejects_out_of_scope_count():
    with pytest.raises(ValueError):
        score_cross_cancer_recurrence(4)


def test_final_score_uses_equal_v1_weights():
    scores = ComponentScores(
        de_score=1.0,
        cancer_association_score=0.8,
        clinical_score=0.5,
        cross_cancer_score=1.0,
    )

    result = calculate_final_score(scores)

    assert result == pytest.approx(0.825)


def test_final_score_is_none_when_component_unavailable():
    scores = ComponentScores(
        de_score=1.0,
        cancer_association_score=None,
        clinical_score=0.5,
        cross_cancer_score=1.0,
    )

    assert calculate_final_score(scores) is None


def test_final_score_all_zero_is_zero():
    scores = ComponentScores(
        de_score=0.0,
        cancer_association_score=0.0,
        clinical_score=0.0,
        cross_cancer_score=0.0,
    )

    assert calculate_final_score(scores) == pytest.approx(
        0.0
    )