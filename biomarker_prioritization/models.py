"""Core data models for BLC Mark Phase 5 biomarker prioritization."""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path


PRIORITIZATION_MODELS_VERSION = "1.0"

DE_WEIGHT = 0.25
CANCER_ASSOCIATION_WEIGHT = 0.25
CLINICAL_WEIGHT = 0.25
CROSS_CANCER_WEIGHT = 0.25


class EvidenceAvailability(str, Enum):
    """Availability state of a Phase 5 evidence component."""

    AVAILABLE = "available"
    NO_SUPPORT = "no_support"
    UNAVAILABLE = "unavailable"


@dataclass(frozen=True)
class PrioritizationConfiguration:
    """Complete fixed Version 1 configuration for one cohort-level run."""

    analysis_id: str
    cancer_cohort: str

    phase3_results_path: Path
    phase4_evidence_path: Path
    phase4_metadata_path: Path
    output_dir: Path

    de_weight: float = DE_WEIGHT
    cancer_association_weight: float = CANCER_ASSOCIATION_WEIGHT
    clinical_weight: float = CLINICAL_WEIGHT
    cross_cancer_weight: float = CROSS_CANCER_WEIGHT

    def __post_init__(self) -> None:
        if (
            not isinstance(self.analysis_id, str)
            or not self.analysis_id.strip()
        ):
            raise ValueError(
                "analysis_id must be a non-empty string."
            )

        if (
            not isinstance(self.cancer_cohort, str)
            or not self.cancer_cohort.strip()
        ):
            raise ValueError(
                "cancer_cohort must be a non-empty string."
            )

        for field_name in (
            "phase3_results_path",
            "phase4_evidence_path",
            "phase4_metadata_path",
            "output_dir",
        ):
            if not isinstance(getattr(self, field_name), Path):
                raise TypeError(
                    f"{field_name} must be a pathlib.Path."
                )

        expected_weights = {
            "de_weight": DE_WEIGHT,
            "cancer_association_weight": CANCER_ASSOCIATION_WEIGHT,
            "clinical_weight": CLINICAL_WEIGHT,
            "cross_cancer_weight": CROSS_CANCER_WEIGHT,
        }

        for field_name, expected_value in expected_weights.items():
            value = getattr(self, field_name)

            if isinstance(value, bool) or not isinstance(
                value,
                (int, float),
            ):
                raise TypeError(
                    f"{field_name} must be numeric."
                )

            if value != expected_value:
                raise ValueError(
                    f"{field_name} is fixed at "
                    f"{expected_value} in Version 1."
                )

        if sum(expected_weights.values()) != 1.0:
            raise ValueError(
                "Version 1 prioritization weights must sum to 1.0."
            )


@dataclass(frozen=True)
class PrioritizationInput:
    """
    Raw Phase 3 and Phase 4 inputs for one candidate gene in one cohort.
    """

    gene_id: str
    cancer_cohort: str

    effect_size: float
    effect_size_label: str
    adjusted_p_value: float

    cancer_association_score: float | None
    clinical_category: str | None
    clinical_direction: str | None
    cross_cancer_cohort_count: int

    cancer_association_availability: EvidenceAvailability
    clinical_availability: EvidenceAvailability
    cross_cancer_availability: EvidenceAvailability

    functional_description: str | None = None
    pathway_count: int = 0

    def __post_init__(self) -> None:
        if not isinstance(self.gene_id, str) or not self.gene_id.strip():
            raise ValueError("gene_id must be a non-empty string.")

        if (
            not isinstance(self.cancer_cohort, str)
            or not self.cancer_cohort.strip()
        ):
            raise ValueError("cancer_cohort must be a non-empty string.")

        if isinstance(self.effect_size, bool) or not isinstance(
            self.effect_size, (int, float)
        ):
            raise TypeError("effect_size must be numeric.")

        if (
            not isinstance(self.effect_size_label, str)
            or not self.effect_size_label.strip()
        ):
            raise ValueError(
                "effect_size_label must be a non-empty string."
            )

        if isinstance(self.adjusted_p_value, bool) or not isinstance(
            self.adjusted_p_value, (int, float)
        ):
            raise TypeError("adjusted_p_value must be numeric.")

        if not (0.0 <= self.adjusted_p_value <= 1.0):
            raise ValueError(
                "adjusted_p_value must be between 0 and 1."
            )

        if self.cancer_association_score is not None:
            if (
                isinstance(self.cancer_association_score, bool)
                or not isinstance(
                    self.cancer_association_score,
                    (int, float),
                )
            ):
                raise TypeError(
                    "cancer_association_score must be numeric or None."
                )

            if not (
                0.0 <= self.cancer_association_score <= 1.0
            ):
                raise ValueError(
                    "cancer_association_score must be between 0 and 1."
                )

        if (
            isinstance(self.cross_cancer_cohort_count, bool)
            or not isinstance(
                self.cross_cancer_cohort_count,
                int,
            )
        ):
            raise TypeError(
                "cross_cancer_cohort_count must be an int."
            )

        if self.cross_cancer_cohort_count < 1:
            raise ValueError(
                "cross_cancer_cohort_count must be at least 1."
            )

        for field_name in (
            "cancer_association_availability",
            "clinical_availability",
            "cross_cancer_availability",
        ):
            if not isinstance(
                getattr(self, field_name),
                EvidenceAvailability,
            ):
                raise TypeError(
                    f"{field_name} must be an EvidenceAvailability."
                )

        if isinstance(self.pathway_count, bool) or not isinstance(
            self.pathway_count,
            int,
        ):
            raise TypeError("pathway_count must be an int.")

        if self.pathway_count < 0:
            raise ValueError("pathway_count cannot be negative.")


@dataclass(frozen=True)
class ComponentScores:
    """Normalized Phase 5 evidence scores for one candidate."""

    de_score: float
    cancer_association_score: float | None
    clinical_score: float | None
    cross_cancer_score: float | None

    def __post_init__(self) -> None:
        for field_name in (
            "de_score",
            "cancer_association_score",
            "clinical_score",
            "cross_cancer_score",
        ):
            value = getattr(self, field_name)

            if value is None:
                continue

            if isinstance(value, bool) or not isinstance(
                value,
                (int, float),
            ):
                raise TypeError(
                    f"{field_name} must be numeric or None."
                )

            if not (0.0 <= value <= 1.0):
                raise ValueError(
                    f"{field_name} must be between 0 and 1."
                )


@dataclass(frozen=True)
class PrioritizedBiomarker:
    """Final explainable Phase 5 result for one candidate gene."""

    gene_id: str
    cancer_cohort: str

    raw_input: PrioritizationInput
    component_scores: ComponentScores

    final_score: float | None
    rank: int | None

    scoring_version: str = PRIORITIZATION_MODELS_VERSION

    def __post_init__(self) -> None:
        if not isinstance(self.gene_id, str) or not self.gene_id.strip():
            raise ValueError("gene_id must be a non-empty string.")

        if (
            not isinstance(self.cancer_cohort, str)
            or not self.cancer_cohort.strip()
        ):
            raise ValueError("cancer_cohort must be a non-empty string.")

        if not isinstance(self.raw_input, PrioritizationInput):
            raise TypeError(
                "raw_input must be a PrioritizationInput."
            )

        if not isinstance(self.component_scores, ComponentScores):
            raise TypeError(
                "component_scores must be a ComponentScores."
            )

        if self.final_score is not None:
            if isinstance(self.final_score, bool) or not isinstance(
                self.final_score,
                (int, float),
            ):
                raise TypeError(
                    "final_score must be numeric or None."
                )

            if not (0.0 <= self.final_score <= 1.0):
                raise ValueError(
                    "final_score must be between 0 and 1."
                )

        if self.rank is not None:
            if isinstance(self.rank, bool) or not isinstance(
                self.rank,
                int,
            ):
                raise TypeError("rank must be an int or None.")

            if self.rank < 1:
                raise ValueError("rank must be at least 1.")