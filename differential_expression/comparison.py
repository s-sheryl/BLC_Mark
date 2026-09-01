"""
Purpose:
    Resolve the explicit, caller-configured two-group comparison
    against validated, sample-matched metadata: determine which
    matched samples belong to the reference group, which belong to
    the comparison group, and which matched samples belong to neither
    (and must therefore be excluded, with the reason and stage
    recorded).

Responsibilities:
    - Consume the configured group column, reference group, and
      comparison group exactly as given -- never infer them from
      column order, sample names, or any other implicit signal.
    - Restrict group resolution to samples already confirmed to exist
      in both the expression matrix and the metadata (the matched
      samples produced by validation.match_samples).
    - Verify both configured groups actually have at least one
      matched sample.
    - Enforce the configured minimum replication requirement per
      group.
    - Preserve deterministic sample ordering (sorted) so the same
      inputs always produce the same GroupAssignment.
    - Record every excluded sample (whether excluded because it only
      appeared in one of the two input files, or because its group
      label matched neither configured group) with its reason and
      the stage at which exclusion occurred.

Scope:
    This module never drops, merges, or renames groups, and never
    swaps reference/comparison direction. It performs no statistical
    computation and no gene-level operations.
"""

from src.differential_expression.exceptions import (
    InsufficientReplicationError,
    InvalidConfigurationError,
    SampleMismatchError,
)
from src.differential_expression.models import (
    DEAnalysisConfiguration,
    ExcludedSample,
    GroupAssignment,
    SampleMatchResult,
    ValidatedMetadata,
)

COMPARISON_VERSION = "1.0"

__all__ = [
    "COMPARISON_VERSION",
    "resolve_comparison",
]


def resolve_comparison(
    configuration: DEAnalysisConfiguration,
    metadata: ValidatedMetadata,
    sample_match: SampleMatchResult,
) -> GroupAssignment:
    """Resolve the configured two-group comparison against matched
    samples.

    Args:
        configuration: The analysis configuration, specifying
            group_column, reference_group, comparison_group, and
            minimum_replicates_per_group.
        metadata: Validated sample metadata (must have been validated
            using the same group_column as `configuration`).
        sample_match: The result of matching expression-matrix samples
            against metadata samples.

    Returns:
        A GroupAssignment recording, with deterministic sample
        ordering, which matched samples belong to the reference group,
        which belong to the comparison group, and a complete record of
        every excluded sample (validation-stage: expression-only or
        metadata-only; comparison-stage: matched but belonging to
        neither configured group).

    Raises:
        InvalidConfigurationError: If `metadata.group_column` does not
            match `configuration.group_column` (an interface-usage
            error, not a scientific decision -- validation.py must be
            called with the same group_column that configuration
            specifies).
        SampleMismatchError: If either the configured reference group
            or comparison group has zero matched samples.
        InsufficientReplicationError: If either group's matched sample
            count is below configuration.minimum_replicates_per_group.
    """
    if metadata.group_column != configuration.group_column:
        raise InvalidConfigurationError(
            "Metadata was validated with group_column "
            f"'{metadata.group_column}' but the configuration specifies "
            f"group_column '{configuration.group_column}'; validation "
            "must be performed using the configured group_column before "
            "resolve_comparison() is called."
        )

    matched_samples = set(sample_match.matched_samples)

    group_labels_by_sample = (
        metadata.dataframe.set_index(metadata.sample_id_column)[
            metadata.group_column
        ]
        .astype(str)
        .to_dict()
    )

    reference_samples = sorted(
        sample_id
        for sample_id in matched_samples
        if group_labels_by_sample.get(sample_id) == configuration.reference_group
    )
    comparison_samples = sorted(
        sample_id
        for sample_id in matched_samples
        if group_labels_by_sample.get(sample_id) == configuration.comparison_group
    )

    if not reference_samples:
        raise SampleMismatchError(
            "The configured reference group "
            f"'{configuration.reference_group}' has zero matched samples "
            "(samples present in both the expression matrix and the "
            "metadata)."
        )

    if not comparison_samples:
        raise SampleMismatchError(
            "The configured comparison group "
            f"'{configuration.comparison_group}' has zero matched samples "
            "(samples present in both the expression matrix and the "
            "metadata)."
        )

    if len(reference_samples) < configuration.minimum_replicates_per_group:
        raise InsufficientReplicationError(
            f"Reference group '{configuration.reference_group}' has "
            f"{len(reference_samples)} matched sample(s), below the "
            "configured minimum of "
            f"{configuration.minimum_replicates_per_group}."
        )

    if len(comparison_samples) < configuration.minimum_replicates_per_group:
        raise InsufficientReplicationError(
            f"Comparison group '{configuration.comparison_group}' has "
            f"{len(comparison_samples)} matched sample(s), below the "
            "configured minimum of "
            f"{configuration.minimum_replicates_per_group}."
        )

    excluded_samples: list[ExcludedSample] = []

    for sample_id in sample_match.expression_only_samples:
        excluded_samples.append(
            ExcludedSample(
                sample_id=sample_id,
                reason=(
                    "Sample present in the expression matrix but not "
                    "found in the metadata."
                ),
                stage="validation",
            )
        )

    for sample_id in sample_match.metadata_only_samples:
        excluded_samples.append(
            ExcludedSample(
                sample_id=sample_id,
                reason=(
                    "Sample present in the metadata but not found in "
                    "the expression matrix."
                ),
                stage="validation",
            )
        )

    unassigned_matched_samples = (
        matched_samples - set(reference_samples) - set(comparison_samples)
    )

    for sample_id in sorted(unassigned_matched_samples):
        label = group_labels_by_sample.get(sample_id)
        excluded_samples.append(
            ExcludedSample(
                sample_id=sample_id,
                reason=(
                    f"Sample's group label ({label!r}) matches neither "
                    f"the configured reference group "
                    f"({configuration.reference_group!r}) nor the "
                    f"configured comparison group "
                    f"({configuration.comparison_group!r})."
                ),
                stage="comparison",
            )
        )

    excluded_samples.sort(key=lambda excluded: excluded.sample_id)

    return GroupAssignment(
        reference_group=configuration.reference_group,
        comparison_group=configuration.comparison_group,
        reference_sample_ids=tuple(reference_samples),
        comparison_sample_ids=tuple(comparison_samples),
        excluded_samples=tuple(excluded_samples),
    )