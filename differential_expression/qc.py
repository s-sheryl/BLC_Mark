"""
Purpose:
    Assemble the structured Quality Control report for one
    differential expression analysis, per specification Section 8.

Responsibilities:
    - Aggregate counts and decisions already produced by
      validation.py, comparison.py, filtering.py, and the
      configuration itself into a single QCReport.
    - Distinguish observed data facts (gene/sample counts, excluded
      samples) from analysis-configuration decisions (statistical
      method, multiple-testing method, significance threshold, gene
      filter criterion), per specification Section 8.5.
    - Perform the post-analysis structural checks required by
      Section 8.4: every tested gene has a corresponding result,
      gene identifiers are correctly associated with their results,
      and produced p-values/adjusted p-values are valid numeric
      values or explicitly missing.

Scope:
    This module does not repair or alter anything it inspects. If a
    post-analysis check fails, it raises QualityControlError rather
    than silently correcting the result set -- QC observes and
    reports, it never fixes.
"""

import math

from src.differential_expression.exceptions import QualityControlError
from src.differential_expression.models import (
    DEAnalysisConfiguration,
    GeneFilterResult,
    GeneResult,
    GroupAssignment,
    QCReport,
    SampleMatchResult,
)

QC_VERSION = "1.0"

__all__ = [
    "QC_VERSION",
    "build_qc_report",
    "verify_post_analysis_quality",
]


def build_qc_report(
    configuration: DEAnalysisConfiguration,
    sample_match: SampleMatchResult,
    group_assignment: GroupAssignment,
    gene_filter_result: GeneFilterResult,
) -> QCReport:
    """Assemble the QC report from data already produced upstream.

    Args:
        configuration: The analysis configuration.
        sample_match: The result of matching expression and metadata
            samples (validation.py).
        group_assignment: The resolved two-group comparison
            (comparison.py), including every excluded sample with its
            reason and stage.
        gene_filter_result: Which genes were tested vs. filtered
            (filtering.py).

    Returns:
        A QCReport with observed data facts and configuration
        decisions kept in clearly separate field groups.
    """
    initial_sample_count = len(sample_match.matched_samples) + len(
        sample_match.expression_only_samples
    ) + len(sample_match.metadata_only_samples)

    included_sample_count = len(group_assignment.reference_sample_ids) + len(
        group_assignment.comparison_sample_ids
    )

    return QCReport(
        input_gene_count=gene_filter_result.input_gene_count,
        tested_gene_count=len(gene_filter_result.tested_gene_ids),
        filtered_gene_count=len(gene_filter_result.filtered_gene_ids),
        initial_sample_count=initial_sample_count,
        included_sample_count=included_sample_count,
        excluded_sample_count=len(group_assignment.excluded_samples),
        excluded_samples=group_assignment.excluded_samples,
        reference_group_size=len(group_assignment.reference_sample_ids),
        comparison_group_size=len(group_assignment.comparison_sample_ids),
        gene_filter_criterion=gene_filter_result.criterion_description,
        statistical_method=configuration.statistical_method,
        multiple_testing_method=configuration.multiple_testing_method,
        significance_threshold=configuration.significance_threshold,
    )


def verify_post_analysis_quality(
    gene_results: tuple[GeneResult, ...],
    gene_filter_result: GeneFilterResult,
    configuration: DEAnalysisConfiguration,
) -> None:
    """Verify post-analysis QC conditions required by specification
    Section 8.4.

    Checks performed:
        - The result set contains exactly one GeneResult per gene in
          gene_filter_result.tested_gene_ids -- no gene silently
          dropped, none duplicated, none unexpectedly added.
        - Every gene identifier in the result set is one of the
          tested gene identifiers (correct gene/result association).
        - raw_p_value and adjusted_p_value, where present, are finite
          numeric values in [0, 1] (valid p-values).
        - adjusted_p_value is only present when raw_p_value is also
          present (a gene cannot be corrected without a raw p-value).
        - significant, where present, is consistent with
          adjusted_p_value and configuration.significance_threshold
          (effect-size-direction and configuration-vs-result
          consistency).

    Args:
        gene_results: The final, corrected gene-level results.
        gene_filter_result: Which genes were tested vs. filtered.
        configuration: The analysis configuration.

    Raises:
        QualityControlError: If any post-analysis condition is not
            satisfied.
    """
    expected_gene_ids = set(gene_filter_result.tested_gene_ids)
    actual_gene_ids = [result.gene_id for result in gene_results]

    if len(actual_gene_ids) != len(expected_gene_ids):
        raise QualityControlError(
            "Post-analysis QC failed: expected "
            f"{len(expected_gene_ids)} gene result(s) (one per tested "
            f"gene) but found {len(actual_gene_ids)}."
        )

    if set(actual_gene_ids) != expected_gene_ids:
        missing = expected_gene_ids - set(actual_gene_ids)
        unexpected = set(actual_gene_ids) - expected_gene_ids
        raise QualityControlError(
            "Post-analysis QC failed: gene identifiers in the result "
            f"set do not match the tested gene set. Missing: "
            f"{sorted(missing)}. Unexpected: {sorted(unexpected)}."
        )

    if len(actual_gene_ids) != len(set(actual_gene_ids)):
        raise QualityControlError(
            "Post-analysis QC failed: duplicate gene identifiers "
            "found in the result set."
        )

    for result in gene_results:
        if result.raw_p_value is not None:
            if not isinstance(result.raw_p_value, float) or not math.isfinite(
                result.raw_p_value
            ):
                raise QualityControlError(
                    "Post-analysis QC failed: gene "
                    f"'{result.gene_id}' has a non-finite raw_p_value "
                    f"({result.raw_p_value!r})."
                )
            if not (0.0 <= result.raw_p_value <= 1.0):
                raise QualityControlError(
                    "Post-analysis QC failed: gene "
                    f"'{result.gene_id}' has an out-of-range "
                    f"raw_p_value ({result.raw_p_value!r})."
                )

        if result.adjusted_p_value is not None:
            if result.raw_p_value is None:
                raise QualityControlError(
                    "Post-analysis QC failed: gene "
                    f"'{result.gene_id}' has an adjusted_p_value but "
                    "no raw_p_value; a gene cannot be corrected "
                    "without a raw p-value."
                )
            if not isinstance(
                result.adjusted_p_value, float
            ) or not math.isfinite(result.adjusted_p_value):
                raise QualityControlError(
                    "Post-analysis QC failed: gene "
                    f"'{result.gene_id}' has a non-finite "
                    f"adjusted_p_value ({result.adjusted_p_value!r})."
                )
            if not (0.0 <= result.adjusted_p_value <= 1.0):
                raise QualityControlError(
                    "Post-analysis QC failed: gene "
                    f"'{result.gene_id}' has an out-of-range "
                    f"adjusted_p_value ({result.adjusted_p_value!r})."
                )

            expected_significant = (
                result.adjusted_p_value < configuration.significance_threshold
            )
            if result.significant != expected_significant:
                raise QualityControlError(
                    "Post-analysis QC failed: gene "
                    f"'{result.gene_id}' has significant="
                    f"{result.significant!r}, inconsistent with "
                    f"adjusted_p_value={result.adjusted_p_value!r} and "
                    "the configured significance_threshold="
                    f"{configuration.significance_threshold!r}."
                )
        else:
            if result.significant is not None:
                raise QualityControlError(
                    "Post-analysis QC failed: gene "
                    f"'{result.gene_id}' has significant="
                    f"{result.significant!r} but no adjusted_p_value; "
                    "significance cannot be evaluated without an "
                    "adjusted p-value."
                )