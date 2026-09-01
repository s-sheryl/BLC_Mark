"""Result writing for BLC Mark Phase 5 biomarker prioritization."""

import csv
import json
from pathlib import Path

from src.biomarker_prioritization.models import (
    PrioritizedBiomarker,
)


PRIORITIZATION_OUTPUT_COLUMNS = (
    "gene_id",
    "cancer_cohort",
    "effect_size",
    "effect_size_label",
    "adjusted_p_value",
    "de_score",
    "cancer_association_raw_score",
    "cancer_association_score",
    "cancer_association_availability",
    "clinical_category",
    "clinical_direction",
    "clinical_score",
    "clinical_availability",
    "cross_cancer_cohort_count",
    "cross_cancer_score",
    "cross_cancer_availability",
    "functional_description",
    "pathway_count",
    "final_score",
    "rank",
    "scoring_version",
)


def write_prioritization_results(
    biomarkers: list[PrioritizedBiomarker]
    | tuple[PrioritizedBiomarker, ...],
    output_csv: str | Path,
) -> Path:
    """Write explainable Phase 5 biomarker rankings to deterministic CSV."""

    output_csv = Path(output_csv)

    output_csv.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    rows: list[dict] = []

    for biomarker in biomarkers:
        if not isinstance(
            biomarker,
            PrioritizedBiomarker,
        ):
            raise TypeError(
                "Every biomarker must be a PrioritizedBiomarker."
            )

        raw_input = biomarker.raw_input
        scores = biomarker.component_scores

        rows.append(
            {
                "gene_id": biomarker.gene_id,
                "cancer_cohort": biomarker.cancer_cohort,
                "effect_size": raw_input.effect_size,
                "effect_size_label": (
                    raw_input.effect_size_label
                ),
                "adjusted_p_value": (
                    raw_input.adjusted_p_value
                ),
                "de_score": scores.de_score,
                "cancer_association_raw_score": (
                    raw_input.cancer_association_score
                ),
                "cancer_association_score": (
                    scores.cancer_association_score
                ),
                "cancer_association_availability": (
                    raw_input
                    .cancer_association_availability
                    .value
                ),
                "clinical_category": (
                    raw_input.clinical_category
                ),
                "clinical_direction": (
                    raw_input.clinical_direction
                ),
                "clinical_score": (
                    scores.clinical_score
                ),
                "clinical_availability": (
                    raw_input.clinical_availability.value
                ),
                "cross_cancer_cohort_count": (
                    raw_input.cross_cancer_cohort_count
                ),
                "cross_cancer_score": (
                    scores.cross_cancer_score
                ),
                "cross_cancer_availability": (
                    raw_input
                    .cross_cancer_availability
                    .value
                ),
                "functional_description": (
                    raw_input.functional_description
                ),
                "pathway_count": (
                    raw_input.pathway_count
                ),
                "final_score": (
                    biomarker.final_score
                ),
                "rank": biomarker.rank,
                "scoring_version": (
                    biomarker.scoring_version
                ),
            }
        )

    rows.sort(
        key=lambda row: (
            row["rank"] is None,
            (
                row["rank"]
                if row["rank"] is not None
                else float("inf")
            ),
            str(row["gene_id"]),
        )
    )

    with output_csv.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                PRIORITIZATION_OUTPUT_COLUMNS
            ),
        )

        writer.writeheader()
        writer.writerows(rows)

    return output_csv


def write_json(
    payload: dict,
    output_path: str | Path,
) -> Path:
    """Write deterministic UTF-8 Phase 5 JSON output."""

    output_path = Path(
        output_path
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            payload,
            handle,
            indent=2,
            sort_keys=True,
        )

        handle.write("\n")

    return output_path