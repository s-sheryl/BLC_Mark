"""Deterministic ranking for BLC Mark Phase 5 biomarker prioritization."""

from dataclasses import replace
from collections.abc import Sequence

from src.biomarker_prioritization.models import PrioritizedBiomarker


RANKING_VERSION = "1.0"


def rank_biomarkers(
    biomarkers: Sequence[PrioritizedBiomarker],
) -> tuple[PrioritizedBiomarker, ...]:
    """
    Rank Phase 5 biomarker candidates deterministically.

    Rules
    -----
    1. Candidates with a numeric final_score are ranked before candidates
       whose final_score is unavailable.
    2. Numeric final_score is sorted from highest to lowest.
    3. Exact score ties are ordered by gene_id ascending.
    4. The lexical tie-break is for reproducibility only and does not imply
       additional biological evidence.
    5. Candidates with final_score=None are retained with rank=None.
    """

    validated: list[PrioritizedBiomarker] = []

    for biomarker in biomarkers:
        if not isinstance(biomarker, PrioritizedBiomarker):
            raise TypeError(
                "Every biomarker must be a PrioritizedBiomarker."
            )

        validated.append(biomarker)

    scored = [
        biomarker
        for biomarker in validated
        if biomarker.final_score is not None
    ]

    unavailable = [
        biomarker
        for biomarker in validated
        if biomarker.final_score is None
    ]

    scored_sorted = sorted(
        scored,
        key=lambda biomarker: (
            -biomarker.final_score,
            biomarker.gene_id,
        ),
    )

    unavailable_sorted = sorted(
        unavailable,
        key=lambda biomarker: biomarker.gene_id,
    )

    ranked: list[PrioritizedBiomarker] = []

    for rank, biomarker in enumerate(
        scored_sorted,
        start=1,
    ):
        ranked.append(
            replace(
                biomarker,
                rank=rank,
            )
        )

    for biomarker in unavailable_sorted:
        ranked.append(
            replace(
                biomarker,
                rank=None,
            )
        )

    return tuple(ranked)