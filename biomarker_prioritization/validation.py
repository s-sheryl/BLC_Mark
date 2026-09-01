"""Input validation for BLC Mark Phase 5 biomarker prioritization."""

import json
from pathlib import Path

import pandas as pd

from src.biomarker_prioritization.exceptions import (
    PrioritizationValidationError,
)
from src.biomarker_prioritization.models import (
    PrioritizationConfiguration,
)


VALIDATION_VERSION = "1.0"


PHASE3_REQUIRED_COLUMNS = {
    "gene_id",
    "tested",
    "effect_size",
    "effect_size_label",
    "adjusted_p_value",
    "significant",
}

PHASE4_REQUIRED_COLUMNS = {
    "gene_id",
    "cancer_cohort",
    "evidence_type",
    "source",
    "source_version",
    "evidence_id",
    "description",
    "retrieved_at",
    "source_url",
}


def _require_existing_file(
    path: Path,
    field_name: str,
) -> None:
    if not isinstance(path, Path):
        raise PrioritizationValidationError(
            f"{field_name} must be a pathlib.Path."
        )

    if not path.exists():
        raise PrioritizationValidationError(
            f"{field_name} does not exist: {path}"
        )

    if not path.is_file():
        raise PrioritizationValidationError(
            f"{field_name} must reference a file: {path}"
        )


def _load_phase3_results(
    path: Path,
) -> pd.DataFrame:
    try:
        dataframe = pd.read_csv(path)
    except Exception as error:
        raise PrioritizationValidationError(
            f"Could not read Phase 3 results: {error}"
        ) from error

    missing = PHASE3_REQUIRED_COLUMNS - set(dataframe.columns)

    if missing:
        raise PrioritizationValidationError(
            "Phase 3 results are missing required columns: "
            f"{sorted(missing)}."
        )

    return dataframe


def _load_phase4_evidence(
    path: Path,
) -> pd.DataFrame:
    try:
        dataframe = pd.read_csv(path)
    except Exception as error:
        raise PrioritizationValidationError(
            f"Could not read Phase 4 evidence: {error}"
        ) from error

    missing = PHASE4_REQUIRED_COLUMNS - set(dataframe.columns)

    if missing:
        raise PrioritizationValidationError(
            "Phase 4 evidence is missing required columns: "
            f"{sorted(missing)}."
        )

    return dataframe


def _load_phase4_metadata(
    path: Path,
) -> dict:
    try:
        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:
            metadata = json.load(handle)
    except Exception as error:
        raise PrioritizationValidationError(
            f"Could not read Phase 4 metadata: {error}"
        ) from error

    if not isinstance(metadata, dict):
        raise PrioritizationValidationError(
            "Phase 4 metadata must contain a JSON object."
        )

    return metadata


def validate_prioritization_inputs(
    configuration: PrioritizationConfiguration,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Validate Phase 3 and Phase 4 inputs for one Phase 5 cohort run.

    Returns validated, loaded inputs only after all structural and
    cross-phase consistency checks succeed.
    """

    if not isinstance(
        configuration,
        PrioritizationConfiguration,
    ):
        raise TypeError(
            "configuration must be a PrioritizationConfiguration."
        )

    _require_existing_file(
        configuration.phase3_results_path,
        "phase3_results_path",
    )
    _require_existing_file(
        configuration.phase4_evidence_path,
        "phase4_evidence_path",
    )
    _require_existing_file(
        configuration.phase4_metadata_path,
        "phase4_metadata_path",
    )

    phase3 = _load_phase3_results(
        configuration.phase3_results_path
    )

    phase4 = _load_phase4_evidence(
        configuration.phase4_evidence_path
    )

    metadata = _load_phase4_metadata(
        configuration.phase4_metadata_path
    )

    metadata_phase = metadata.get("phase")

    if metadata_phase != 4:
        raise PrioritizationValidationError(
            "Phase 4 metadata must report phase=4; "
            f"got {metadata_phase!r}."
        )

    metadata_cohort = metadata.get("cancer_cohort")

    if metadata_cohort != configuration.cancer_cohort:
        raise PrioritizationValidationError(
            "Configured cancer cohort does not match Phase 4 metadata: "
            f"{configuration.cancer_cohort!r} != "
            f"{metadata_cohort!r}."
        )

    if not phase4.empty:
        evidence_cohorts = set(
            phase4["cancer_cohort"]
            .dropna()
            .astype(str)
            .unique()
        )

        if evidence_cohorts != {
            configuration.cancer_cohort
        }:
            raise PrioritizationValidationError(
                "Phase 4 evidence contains unexpected cancer cohorts: "
                f"{sorted(evidence_cohorts)}."
            )

    significant_mask = (
        phase3["significant"]
        .astype(str)
        .str.lower()
        .eq("true")
    )

    candidates = phase3.loc[
        significant_mask
    ].copy()

    if candidates.empty:
        raise PrioritizationValidationError(
            "Phase 3 results contain no significant candidates."
        )

    if candidates["gene_id"].isna().any():
        raise PrioritizationValidationError(
            "Phase 3 significant candidates contain missing gene_id values."
        )

    duplicated_candidates = candidates[
        candidates["gene_id"].duplicated(
            keep=False
        )
    ]

    if not duplicated_candidates.empty:
        duplicate_ids = sorted(
            duplicated_candidates[
                "gene_id"
            ].astype(str).unique()
        )

        raise PrioritizationValidationError(
            "Phase 3 significant candidates contain duplicate gene IDs: "
            f"{duplicate_ids[:10]}."
        )

    candidate_count = metadata.get(
        "candidate_count"
    )

    if not isinstance(candidate_count, int):
        raise PrioritizationValidationError(
            "Phase 4 metadata candidate_count must be an integer."
        )

    if candidate_count != len(candidates):
        raise PrioritizationValidationError(
            "Phase 3 significant-candidate count does not match "
            "Phase 4 metadata candidate_count: "
            f"{len(candidates)} != {candidate_count}."
        )

    evidence_gene_ids = set(
        phase4["gene_id"]
        .dropna()
        .astype(str)
    )

    candidate_gene_ids = set(
        candidates["gene_id"]
        .astype(str)
    )

    unexpected_evidence_genes = (
        evidence_gene_ids
        - candidate_gene_ids
    )

    if unexpected_evidence_genes:
        raise PrioritizationValidationError(
            "Phase 4 evidence contains genes that are not significant "
            "Phase 3 candidates: "
            f"{sorted(unexpected_evidence_genes)[:10]}."
        )

    duplicate_evidence_mask = phase4.duplicated(
        subset=[
            "gene_id",
            "evidence_type",
            "source",
            "evidence_id",
        ],
        keep=False,
    )

    if duplicate_evidence_mask.any():
        raise PrioritizationValidationError(
            "Phase 4 evidence contains duplicate evidence records."
        )

    return (
        candidates.reset_index(drop=True),
        phase4.reset_index(drop=True),
        metadata,
    )