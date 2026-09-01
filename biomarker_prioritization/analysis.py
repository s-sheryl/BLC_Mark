"""End-to-end orchestration for BLC Mark Phase 5 biomarker prioritization."""

import re
from collections import defaultdict

import pandas as pd

from src.biomarker_prioritization.models import (
    ComponentScores,
    EvidenceAvailability,
    PrioritizationConfiguration,
    PrioritizationInput,
    PrioritizedBiomarker,
)
from src.biomarker_prioritization.ranking import rank_biomarkers
from src.biomarker_prioritization.scoring import (
    calculate_de_percentile_scores,
    calculate_final_score,
    score_cancer_association,
    score_clinical_category,
    score_cross_cancer_recurrence,
)
from src.biomarker_prioritization.validation import (
    validate_prioritization_inputs,
)
from src.evidence_integration.identifiers import normalize_gene_symbol


ANALYSIS_VERSION = "1.0"


_OPEN_TARGETS_SCORE_PATTERN = re.compile(
    r"association score=([0-9]*\.?[0-9]+)"
)

_CROSS_CANCER_PATTERN = re.compile(
    r"candidate in (\d+) cohorts"
)


def _records_by_gene(
    evidence: pd.DataFrame,
) -> dict[str, list[dict]]:
    """Group Phase 4 long-format evidence rows by gene ID."""

    grouped: dict[str, list[dict]] = defaultdict(list)

    for record in evidence.to_dict(
        orient="records"
    ):
        grouped[str(record["gene_id"])].append(
            record
        )

    return dict(grouped)


def _records_of_type(
    records: list[dict],
    evidence_type: str,
) -> list[dict]:
    """Return evidence records of one Phase 4 evidence type."""

    return [
        record
        for record in records
        if record["evidence_type"] == evidence_type
    ]


def _single_scoring_record(
    records: list[dict],
    evidence_type: str,
) -> dict | None:
    """
    Return the single record for a scoring evidence type.

    Version 1 expects at most one cancer-association, clinical,
    or cross-cancer record per candidate gene.
    """

    matching = _records_of_type(
        records,
        evidence_type,
    )

    if len(matching) > 1:
        raise ValueError(
            f"Gene has multiple {evidence_type!r} "
            "scoring records; Version 1 expects at most one."
        )

    if not matching:
        return None

    return matching[0]


def _parse_cancer_association_score(
    record: dict,
) -> float:
    """Extract the Open Targets association score from Phase 4 evidence."""

    description = str(record["description"])

    match = _OPEN_TARGETS_SCORE_PATTERN.search(
        description
    )

    if match is None:
        raise ValueError(
            "Cancer-association evidence does not contain "
            "an Open Targets association score."
        )

    return score_cancer_association(
        float(match.group(1))
    )


def _parse_clinical_record(
    record: dict,
) -> tuple[str, str | None]:
    """
    Extract HPA prognostic category and direction.

    Example:
        'potential prognostic - unfavorable; reported value=...'
    """

    description = str(record["description"])

    category = description.split(
        ";",
        maxsplit=1,
    )[0].strip()

    if not category:
        raise ValueError(
            "Clinical evidence contains an empty prognostic category."
        )

    direction: str | None = None

    normalized = category.lower()

    if normalized.endswith(" - favorable"):
        direction = "favorable"
    elif normalized.endswith(" - unfavorable"):
        direction = "unfavorable"

    # Validate the category against the frozen scoring rules.
    score_clinical_category(category)

    return category, direction


def _parse_cross_cancer_count(
    record: dict,
) -> int:
    """Extract significant-cohort recurrence count from Phase 4 evidence."""

    description = str(record["description"])

    match = _CROSS_CANCER_PATTERN.search(
        description
    )

    if match is None:
        raise ValueError(
            "Cross-cancer evidence does not contain a cohort count."
        )

    count = int(match.group(1))

    # Validate against the frozen V1 1/2/3-cohort scoring contract.
    score_cross_cancer_recurrence(count)

    return count


def _build_prioritization_input(
    phase3_row: pd.Series,
    records: list[dict],
    cancer_cohort: str,
) -> PrioritizationInput:
    """Build one raw Phase 5 candidate input."""

    gene_id = str(
        phase3_row["gene_id"]
    )

    identifier = normalize_gene_symbol(
        gene_id
    )

    # Phase 4 explicitly does not query unresolved legacy ?|<id>
    # identifiers. Their Phase 3 DE signal remains valid, but Phase 4
    # scoring evidence is unavailable rather than negative.
    if not identifier.resolvable:
        return PrioritizationInput(
            gene_id=gene_id,
            cancer_cohort=cancer_cohort,
            effect_size=float(
                phase3_row["effect_size"]
            ),
            effect_size_label=str(
                phase3_row["effect_size_label"]
            ),
            adjusted_p_value=float(
                phase3_row["adjusted_p_value"]
            ),
            cancer_association_score=None,
            clinical_category=None,
            clinical_direction=None,
            cross_cancer_cohort_count=1,
            cancer_association_availability=(
                EvidenceAvailability.UNAVAILABLE
            ),
            clinical_availability=(
                EvidenceAvailability.UNAVAILABLE
            ),
            cross_cancer_availability=(
                EvidenceAvailability.UNAVAILABLE
            ),
            functional_description=None,
            pathway_count=0,
        )

    cancer_record = _single_scoring_record(
        records,
        "cancer_association",
    )

    if cancer_record is None:
        cancer_score = 0.0
        cancer_availability = (
            EvidenceAvailability.NO_SUPPORT
        )
    else:
        cancer_score = (
            _parse_cancer_association_score(
                cancer_record
            )
        )
        cancer_availability = (
            EvidenceAvailability.AVAILABLE
        )

    clinical_record = _single_scoring_record(
        records,
        "clinical",
    )

    if clinical_record is None:
        clinical_category = None
        clinical_direction = None
        clinical_availability = (
            EvidenceAvailability.NO_SUPPORT
        )
    else:
        (
            clinical_category,
            clinical_direction,
        ) = _parse_clinical_record(
            clinical_record
        )
        clinical_availability = (
            EvidenceAvailability.AVAILABLE
        )

    cross_record = _single_scoring_record(
        records,
        "cross_cancer",
    )

    if cross_record is None:
        # Absence of a cross-cancer record means the gene is significant
        # only in the current cohort.
        cross_cancer_count = 1
        cross_availability = (
            EvidenceAvailability.NO_SUPPORT
        )
    else:
        cross_cancer_count = (
            _parse_cross_cancer_count(
                cross_record
            )
        )
        cross_availability = (
            EvidenceAvailability.AVAILABLE
        )

    functional_records = _records_of_type(
        records,
        "functional",
    )

    if len(functional_records) > 1:
        raise ValueError(
            "Gene has multiple functional records; "
            "Version 1 expects at most one."
        )

    functional_description = (
        str(functional_records[0]["description"])
        if functional_records
        else None
    )

    pathway_count = len(
        _records_of_type(
            records,
            "pathway",
        )
    )

    return PrioritizationInput(
        gene_id=gene_id,
        cancer_cohort=cancer_cohort,
        effect_size=float(
            phase3_row["effect_size"]
        ),
        effect_size_label=str(
            phase3_row["effect_size_label"]
        ),
        adjusted_p_value=float(
            phase3_row["adjusted_p_value"]
        ),
        cancer_association_score=cancer_score,
        clinical_category=clinical_category,
        clinical_direction=clinical_direction,
        cross_cancer_cohort_count=cross_cancer_count,
        cancer_association_availability=(
            cancer_availability
        ),
        clinical_availability=(
            clinical_availability
        ),
        cross_cancer_availability=(
            cross_availability
        ),
        functional_description=(
            functional_description
        ),
        pathway_count=pathway_count,
    )


def run_prioritization_analysis(
    configuration: PrioritizationConfiguration,
) -> tuple[PrioritizedBiomarker, ...]:
    """
    Run Phase 5 scoring and ranking for one cancer cohort.

    Phase 3 significance determines eligibility.
    Phase 4 supplies evidence.
    This function performs no file writing.
    """

    if not isinstance(
        configuration,
        PrioritizationConfiguration,
    ):
        raise TypeError(
            "configuration must be a PrioritizationConfiguration."
        )

    (
        candidates,
        evidence,
        _metadata,
    ) = validate_prioritization_inputs(
        configuration
    )

    grouped_evidence = _records_by_gene(
        evidence
    )

    raw_inputs: list[
        PrioritizationInput
    ] = []

    for _, row in candidates.iterrows():
        gene_id = str(
            row["gene_id"]
        )

        raw_inputs.append(
            _build_prioritization_input(
                row,
                grouped_evidence.get(
                    gene_id,
                    [],
                ),
                configuration.cancer_cohort,
            )
        )

    de_scores = calculate_de_percentile_scores(
        [
            candidate.effect_size
            for candidate in raw_inputs
        ]
    )

    unranked: list[
        PrioritizedBiomarker
    ] = []

    for raw_input, de_score in zip(
        raw_inputs,
        de_scores,
        strict=True,
    ):
        if (
            raw_input.cancer_association_availability
            == EvidenceAvailability.UNAVAILABLE
        ):
            cancer_score = None
        else:
            cancer_score = (
                raw_input.cancer_association_score
            )

        if (
            raw_input.clinical_availability
            == EvidenceAvailability.UNAVAILABLE
        ):
            clinical_score = None
        elif (
            raw_input.clinical_availability
            == EvidenceAvailability.NO_SUPPORT
        ):
            clinical_score = 0.0
        else:
            clinical_score = (
                score_clinical_category(
                    raw_input.clinical_category
                )
            )

        if (
            raw_input.cross_cancer_availability
            == EvidenceAvailability.UNAVAILABLE
        ):
            cross_score = None
        else:
            cross_score = (
                score_cross_cancer_recurrence(
                    raw_input.cross_cancer_cohort_count
                )
            )

        component_scores = ComponentScores(
            de_score=de_score,
            cancer_association_score=(
                cancer_score
            ),
            clinical_score=clinical_score,
            cross_cancer_score=cross_score,
        )

        final_score = calculate_final_score(
            component_scores
        )

        unranked.append(
            PrioritizedBiomarker(
                gene_id=raw_input.gene_id,
                cancer_cohort=(
                    configuration.cancer_cohort
                ),
                raw_input=raw_input,
                component_scores=(
                    component_scores
                ),
                final_score=final_score,
                rank=None,
            )
        )

    return rank_biomarkers(
        unranked
    )