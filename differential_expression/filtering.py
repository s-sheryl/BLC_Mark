"""
Purpose:
    Apply the explicitly configured gene-level filter (if any) before
    statistical testing, and record exactly which genes were removed
    and why.

Responsibilities:
    - Apply GeneFilterConfiguration to the validated expression
      matrix, restricted to the samples included in the resolved
      comparison (reference + comparison group samples only -- a
      sample already excluded by comparison.py must not influence
      which genes are filtered).
    - When apply_filter is False (Version 1 default), filter zero
      genes and record that fact.
    - When apply_filter is True, remove genes whose mean expression
      across included samples is strictly below
      minimum_mean_expression, and record the criterion applied.
    - Never invent a filtering criterion beyond the single one defined
      in GeneFilterConfiguration (mean-expression threshold). If the
      configuration does not request filtering, none is applied.

Scope:
    This module performs no statistical testing. It only decides
    which genes proceed to methods.py and records that decision for
    QC and reproducibility purposes.
"""

from src.differential_expression.models import (
    GeneFilterConfiguration,
    GeneFilterResult,
    GroupAssignment,
    ValidatedExpressionMatrix,
)

FILTERING_VERSION = "1.0"

__all__ = [
    "FILTERING_VERSION",
    "apply_gene_filter",
]


def apply_gene_filter(
    gene_filter: GeneFilterConfiguration,
    expression: ValidatedExpressionMatrix,
    group_assignment: GroupAssignment,
) -> GeneFilterResult:
    """Apply the configured gene filter to the validated expression
    matrix, restricted to the samples included in the comparison.

    Args:
        gene_filter: The explicit filtering configuration.
        expression: The validated expression matrix.
        group_assignment: The resolved comparison, whose reference and
            comparison sample IDs define which columns are considered
            when computing mean expression for filtering.

    Returns:
        A GeneFilterResult recording which genes are tested and which
        were filtered, in original expression-matrix order, with the
        criterion applied.
    """
    ordered_gene_ids = list(expression.gene_ids)
    input_gene_count = len(ordered_gene_ids)

    if not gene_filter.apply_filter:
        return GeneFilterResult(
            tested_gene_ids=tuple(ordered_gene_ids),
            filtered_gene_ids=(),
            criterion_description=gene_filter.criterion_description,
            input_gene_count=input_gene_count,
        )

    included_samples = list(group_assignment.reference_sample_ids) + list(
        group_assignment.comparison_sample_ids
    )

    mean_expression = expression.dataframe[included_samples].mean(axis=1)

    threshold = gene_filter.minimum_mean_expression
    keep_mask = mean_expression >= threshold

    tested_gene_ids = tuple(
        gene_id
        for gene_id, keep in zip(ordered_gene_ids, keep_mask, strict=True)
        if keep
    )
    filtered_gene_ids = tuple(
        gene_id
        for gene_id, keep in zip(ordered_gene_ids, keep_mask, strict=True)
        if not keep
    )

    return GeneFilterResult(
        tested_gene_ids=tested_gene_ids,
        filtered_gene_ids=filtered_gene_ids,
        criterion_description=gene_filter.criterion_description,
        input_gene_count=input_gene_count,
    )