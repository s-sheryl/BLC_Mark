"""
Purpose:
    Construct and validate a DEAnalysisConfiguration from raw,
    caller-supplied values (as would arrive from a CLI, config file,
    or calling script), translating them into the typed models defined
    in models.py and raising a single, consistent DE exception type
    (InvalidConfigurationError) for every configuration problem.

Responsibilities:
    - Accept raw scalar/string values for enum-typed fields
      (cancer_cohort, expression_representation, statistical_method,
      multiple_testing_method) and resolve them to the corresponding
      enum members, or fail explicitly if the value is not a
      recognized member.
    - Delegate all scientific validation (identical groups, threshold
      ranges, minimum replication floor, etc.) to
      DEAnalysisConfiguration.__post_init__ in models.py, rather than
      re-implementing it here.
    - Re-raise any TypeError/ValueError raised during construction as
      InvalidConfigurationError, so every module downstream of
      configuration.py only needs to catch DE-specific exceptions.

Scope:
    This module performs no file I/O -- it does not check whether
    expression_matrix_path or metadata_path actually exist on disk
    (that is validation.py's responsibility, and validation.py's
    checks run against the paths only after configuration has
    resolved and validated them). This module also does not resolve
    the comparison itself (that is comparison.py's responsibility) --
    it only validates that reference_group and comparison_group are
    distinct, non-empty strings, which DEAnalysisConfiguration already
    enforces.
"""

from pathlib import Path
from typing import Any

from src.dataset_registry import CancerType
from src.differential_expression.exceptions import InvalidConfigurationError
from src.differential_expression.models import (
    DEAnalysisConfiguration,
    ExpressionRepresentation,
    GeneFilterConfiguration,
    MultipleTestingMethod,
    StatisticalMethod,
)

CONFIGURATION_VERSION = "1.0"

__all__ = [
    "CONFIGURATION_VERSION",
    "build_configuration",
]


def _resolve_enum(
    enum_cls: type,
    value: Any,
    field_name: str,
) -> Any:
    """Resolve a raw value (already an enum member, or its string
    value) to a member of `enum_cls`.

    Args:
        enum_cls: The enum class to resolve against.
        value: Either an existing member of `enum_cls`, or a string
            equal to one of its `.value`s.
        field_name: Name of the configuration field, used only in
            error messages.

    Returns:
        The resolved enum member.

    Raises:
        InvalidConfigurationError: If `value` is not a member of
            `enum_cls` and is not a string matching any member's
            value.
    """
    if isinstance(value, enum_cls):
        return value

    if isinstance(value, str):
        try:
            return enum_cls(value)
        except ValueError:
            valid_values = sorted(member.value for member in enum_cls)
            raise InvalidConfigurationError(
                f"'{field_name}' value {value!r} is not a recognized "
                f"{enum_cls.__name__}. Valid values: {valid_values}."
            ) from None

    raise InvalidConfigurationError(
        f"'{field_name}' must be a {enum_cls.__name__} or a matching "
        f"str, got {type(value).__name__}."
    )


def _resolve_path(value: Any, field_name: str) -> Path:
    """Resolve a raw value to a pathlib.Path.

    Args:
        value: Either a Path already, or a str/os.PathLike.
        field_name: Name of the configuration field, used only in
            error messages.

    Returns:
        The resolved Path.

    Raises:
        InvalidConfigurationError: If `value` cannot be interpreted
            as a path.
    """
    if isinstance(value, Path):
        return value

    if isinstance(value, str):
        if not value.strip():
            raise InvalidConfigurationError(
                f"'{field_name}' must not be an empty string."
            )
        return Path(value)

    raise InvalidConfigurationError(
        f"'{field_name}' must be a pathlib.Path or str, "
        f"got {type(value).__name__}."
    )


def build_configuration(
    *,
    analysis_id: str,
    cancer_cohort: CancerType | str,
    expression_matrix_path: Path | str,
    metadata_path: Path | str,
    gene_id_column: str,
    sample_id_column: str,
    group_column: str,
    reference_group: str,
    comparison_group: str,
    expression_representation: ExpressionRepresentation | str,
    statistical_method: StatisticalMethod | str,
    minimum_replicates_per_group: int,
    output_dir: Path | str,
    multiple_testing_method: MultipleTestingMethod | str = (
        MultipleTestingMethod.BENJAMINI_HOCHBERG
    ),
    significance_threshold: float = 0.05,
    effect_size_threshold: float | None = None,
    gene_filter: GeneFilterConfiguration | None = None,
    random_seed: int | None = None,
) -> DEAnalysisConfiguration:
    """Build and validate a DEAnalysisConfiguration from raw inputs.

    Every value is accepted either as its final typed form (enum
    member, Path) or as the raw form it would take arriving from a
    CLI/config file (matching str, str path). Every failure -- an
    unrecognized enum value, a bad path type, or any structural rule
    enforced by DEAnalysisConfiguration.__post_init__ (e.g. identical
    reference/comparison groups, an out-of-range threshold, a
    replication minimum below the mathematical floor) -- is reported
    as InvalidConfigurationError, so callers only need to catch one
    exception type for configuration problems.

    Args:
        See DEAnalysisConfiguration in models.py for the meaning of
        every field.

    Returns:
        A validated DEAnalysisConfiguration.

    Raises:
        InvalidConfigurationError: If any value cannot be resolved to
            its typed form, or if the resulting configuration is
            scientifically invalid.
    """
    resolved_cancer_cohort = _resolve_enum(CancerType, cancer_cohort, "cancer_cohort")
    resolved_representation = _resolve_enum(
        ExpressionRepresentation,
        expression_representation,
        "expression_representation",
    )
    resolved_method = _resolve_enum(
        StatisticalMethod, statistical_method, "statistical_method"
    )
    resolved_testing_method = _resolve_enum(
        MultipleTestingMethod,
        multiple_testing_method,
        "multiple_testing_method",
    )

    resolved_expression_path = _resolve_path(
        expression_matrix_path, "expression_matrix_path"
    )
    resolved_metadata_path = _resolve_path(metadata_path, "metadata_path")
    resolved_output_dir = _resolve_path(output_dir, "output_dir")

    resolved_gene_filter = (
        gene_filter if gene_filter is not None else GeneFilterConfiguration()
    )

    if not isinstance(resolved_gene_filter, GeneFilterConfiguration):
        raise InvalidConfigurationError(
            "'gene_filter' must be a GeneFilterConfiguration or None, "
            f"got {type(resolved_gene_filter).__name__}."
        )

    try:
        return DEAnalysisConfiguration(
            analysis_id=analysis_id,
            cancer_cohort=resolved_cancer_cohort,
            expression_matrix_path=resolved_expression_path,
            metadata_path=resolved_metadata_path,
            gene_id_column=gene_id_column,
            sample_id_column=sample_id_column,
            group_column=group_column,
            reference_group=reference_group,
            comparison_group=comparison_group,
            expression_representation=resolved_representation,
            statistical_method=resolved_method,
            minimum_replicates_per_group=minimum_replicates_per_group,
            output_dir=resolved_output_dir,
            multiple_testing_method=resolved_testing_method,
            significance_threshold=significance_threshold,
            effect_size_threshold=effect_size_threshold,
            gene_filter=resolved_gene_filter,
            random_seed=random_seed,
        )
    except (TypeError, ValueError) as error:
        raise InvalidConfigurationError(
            f"Invalid differential expression configuration: {error}"
        ) from error