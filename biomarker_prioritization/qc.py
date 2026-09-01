"""Quality control for BLC Mark Phase 5 biomarker prioritization."""

from collections import Counter

from src.biomarker_prioritization.models import (
    EvidenceAvailability,
    PrioritizedBiomarker,
)


def build_phase5_qc_report(
    biomarkers: list[PrioritizedBiomarker]
    | tuple[PrioritizedBiomarker, ...],
    *,
    candidate_count: int,
) -> dict:
    """Build deterministic Phase 5 prioritization QC statistics."""

    if (
        isinstance(candidate_count, bool)
        or not isinstance(
            candidate_count,
            int,
        )
    ):
        raise TypeError(
            "candidate_count must be an int."
        )

    if candidate_count < 0:
        raise ValueError(
            "candidate_count cannot be negative."
        )

    if len(biomarkers) != candidate_count:
        raise ValueError(
            "Biomarker count does not match candidate_count."
        )

    scored_candidate_count = 0
    unavailable_score_count = 0

    availability_counts: Counter[str] = (
        Counter()
    )

    clinical_category_counts: Counter[
        str
    ] = Counter()

    cross_cancer_counts: Counter[
        str
    ] = Counter()

    final_scores: list[float] = []

    for biomarker in biomarkers:
        if not isinstance(
            biomarker,
            PrioritizedBiomarker,
        ):
            raise TypeError(
                "Every biomarker must be a PrioritizedBiomarker."
            )

        if biomarker.final_score is None:
            unavailable_score_count += 1
        else:
            scored_candidate_count += 1
            final_scores.append(
                biomarker.final_score
            )

        raw_input = biomarker.raw_input

        for availability in (
            raw_input.cancer_association_availability,
            raw_input.clinical_availability,
            raw_input.cross_cancer_availability,
        ):
            if not isinstance(
                availability,
                EvidenceAvailability,
            ):
                raise TypeError(
                    "Evidence availability values must be "
                    "EvidenceAvailability instances."
                )

            availability_counts[
                availability.value
            ] += 1

        if (
            raw_input.clinical_category
            is not None
        ):
            clinical_category_counts[
                raw_input.clinical_category
            ] += 1

        cross_cancer_counts[
            str(
                raw_input.cross_cancer_cohort_count
            )
        ] += 1

    if final_scores:
        minimum_final_score = min(
            final_scores
        )
        maximum_final_score = max(
            final_scores
        )
    else:
        minimum_final_score = None
        maximum_final_score = None

    return {
        "candidate_count": candidate_count,
        "prioritized_biomarker_count": (
            len(biomarkers)
        ),
        "scored_candidate_count": (
            scored_candidate_count
        ),
        "unavailable_score_count": (
            unavailable_score_count
        ),
        "minimum_final_score": (
            minimum_final_score
        ),
        "maximum_final_score": (
            maximum_final_score
        ),
        "evidence_availability_counts": dict(
            sorted(
                availability_counts.items()
            )
        ),
        "clinical_category_counts": dict(
            sorted(
                clinical_category_counts.items()
            )
        ),
        "cross_cancer_cohort_counts": dict(
            sorted(
                cross_cancer_counts.items()
            )
        ),
    }