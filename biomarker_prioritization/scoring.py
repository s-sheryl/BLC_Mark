"""Scoring functions for BLC Mark Phase 5 biomarker prioritization."""

from collections.abc import Sequence
from math import isfinite

from src.biomarker_prioritization.models import (
    CANCER_ASSOCIATION_WEIGHT,
    CLINICAL_WEIGHT,
    CROSS_CANCER_WEIGHT,
    DE_WEIGHT,
    ComponentScores,
)


SCORING_VERSION = "1.0"


def calculate_de_percentile_scores(
    effect_sizes: Sequence[float],
) -> tuple[float, ...]:
    """
    Convert absolute effect-size magnitudes to within-cohort percentile ranks.

    Percentiles use average ranks for ties and are scaled to the closed
    interval [0, 1].

    For a cohort containing a single candidate, the candidate receives 1.0.
    """

    if not effect_sizes:
        raise ValueError("effect_sizes must contain at least one value.")

    magnitudes: list[float] = []

    for value in effect_sizes:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TypeError("Every effect size must be numeric.")

        numeric_value = float(value)

        if not isfinite(numeric_value):
            raise ValueError("Effect sizes must be finite.")

        magnitudes.append(abs(numeric_value))

    if len(magnitudes) == 1:
        return (1.0,)

    indexed = sorted(
        enumerate(magnitudes),
        key=lambda item: (item[1], item[0]),
    )

    ranks = [0.0] * len(magnitudes)

    position = 0

    while position < len(indexed):
        end = position

        while (
            end + 1 < len(indexed)
            and indexed[end + 1][1] == indexed[position][1]
        ):
            end += 1

        average_rank = ((position + 1) + (end + 1)) / 2.0
        percentile = (average_rank - 1.0) / (len(indexed) - 1.0)

        for tied_position in range(position, end + 1):
            original_index = indexed[tied_position][0]
            ranks[original_index] = percentile

        position = end + 1

    return tuple(ranks)


def score_cancer_association(
    association_score: float | None,
) -> float | None:
    """
    Return the Open Targets association score unchanged.

    None represents unavailable evidence, not biological negative evidence.
    """

    if association_score is None:
        return None

    if isinstance(association_score, bool) or not isinstance(
        association_score,
        (int, float),
    ):
        raise TypeError(
            "association_score must be numeric or None."
        )

    value = float(association_score)

    if not isfinite(value):
        raise ValueError("association_score must be finite.")

    if not (0.0 <= value <= 1.0):
        raise ValueError(
            "association_score must be between 0 and 1."
        )

    return value


def score_clinical_category(
    clinical_category: str | None,
) -> float | None:
    """
    Map HPA prognostic evidence strength to the Version 1 score.

    Favorable versus unfavorable direction does not change evidence strength.
    """

    if clinical_category is None:
        return None

    if not isinstance(clinical_category, str):
        raise TypeError(
            "clinical_category must be a string or None."
        )

    normalized = clinical_category.strip().lower()

    if not normalized:
        raise ValueError(
            "clinical_category must not be an empty string."
        )

    if normalized.startswith("unprognostic"):
        return 0.0

    if normalized.startswith("potential prognostic"):
        return 0.5

    if normalized.startswith("validated prognostic"):
        return 1.0

    raise ValueError(
        f"Unsupported clinical category: {clinical_category!r}."
    )


def score_cross_cancer_recurrence(
    cohort_count: int,
) -> float:
    """
    Map significant-candidate recurrence across the three V1 cohorts.

    1 cohort -> 0.0
    2 cohorts -> 0.5
    3 cohorts -> 1.0
    """

    if isinstance(cohort_count, bool) or not isinstance(
        cohort_count,
        int,
    ):
        raise TypeError("cohort_count must be an int.")

    mapping = {
        1: 0.0,
        2: 0.5,
        3: 1.0,
    }

    if cohort_count not in mapping:
        raise ValueError(
            "cohort_count must be 1, 2, or 3 in Version 1."
        )

    return mapping[cohort_count]


def calculate_final_score(
    component_scores: ComponentScores,
) -> float | None:
    """
    Calculate the fixed-weight Version 1 prioritization score.

    If any component is unavailable (None), no complete final score is
    produced. The candidate is retained and the unavailable component remains
    explicit rather than being silently treated as zero.
    """

    if not isinstance(component_scores, ComponentScores):
        raise TypeError(
            "component_scores must be a ComponentScores."
        )

    values = (
        component_scores.de_score,
        component_scores.cancer_association_score,
        component_scores.clinical_score,
        component_scores.cross_cancer_score,
    )

    if any(value is None for value in values):
        return None

    return (
        DE_WEIGHT * component_scores.de_score
        + CANCER_ASSOCIATION_WEIGHT
        * component_scores.cancer_association_score
        + CLINICAL_WEIGHT * component_scores.clinical_score
        + CROSS_CANCER_WEIGHT
        * component_scores.cross_cancer_score
    )