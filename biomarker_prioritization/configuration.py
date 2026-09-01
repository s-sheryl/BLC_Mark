"""Configuration construction for BLC Mark Phase 5 prioritization."""

from pathlib import Path
from typing import Any

from src.biomarker_prioritization.models import (
    PrioritizationConfiguration,
)


CONFIGURATION_VERSION = "1.0"


class InvalidPrioritizationConfigurationError(ValueError):
    """Raised when a Phase 5 prioritization configuration is invalid."""


def _resolve_path(value: Any, field_name: str) -> Path:
    """Resolve a raw path-like configuration value to pathlib.Path."""

    if isinstance(value, Path):
        return value

    if isinstance(value, str):
        if not value.strip():
            raise InvalidPrioritizationConfigurationError(
                f"'{field_name}' must not be an empty string."
            )
        return Path(value)

    raise InvalidPrioritizationConfigurationError(
        f"'{field_name}' must be a pathlib.Path or str, "
        f"got {type(value).__name__}."
    )


def build_configuration(
    *,
    analysis_id: str,
    cancer_cohort: str,
    phase3_results_path: Path | str,
    phase4_evidence_path: Path | str,
    phase4_metadata_path: Path | str,
    output_dir: Path | str,
) -> PrioritizationConfiguration:
    """
    Build and validate a Phase 5 PrioritizationConfiguration.

    Version 1 uses fixed equal component weights defined by the frozen
    biomarker prioritization specification. Callers cannot override them.
    """

    resolved_phase3_path = _resolve_path(
        phase3_results_path,
        "phase3_results_path",
    )

    resolved_phase4_evidence_path = _resolve_path(
        phase4_evidence_path,
        "phase4_evidence_path",
    )

    resolved_phase4_metadata_path = _resolve_path(
        phase4_metadata_path,
        "phase4_metadata_path",
    )

    resolved_output_dir = _resolve_path(
        output_dir,
        "output_dir",
    )

    try:
        return PrioritizationConfiguration(
            analysis_id=analysis_id,
            cancer_cohort=cancer_cohort,
            phase3_results_path=resolved_phase3_path,
            phase4_evidence_path=resolved_phase4_evidence_path,
            phase4_metadata_path=resolved_phase4_metadata_path,
            output_dir=resolved_output_dir,
        )
    except (TypeError, ValueError) as error:
        raise InvalidPrioritizationConfigurationError(
            f"Invalid biomarker prioritization configuration: {error}"
        ) from error