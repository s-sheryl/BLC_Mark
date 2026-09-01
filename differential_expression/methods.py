"""
Purpose:
    Provide the statistical-method execution layer for Phase 3, with a
    clear per-method interface, explicit method/representation
    compatibility checks, and no silent method substitution.

Responsibilities:
    - Declare which expression representation(s) each supported
      statistical method is compatible with.
    - Verify method/representation compatibility before attempting to
      execute anything.
    - Execute the configured method and return one GeneResult per
      tested gene, with effect size, raw p-value, and a genuine
      per-method effect-size label -- never a "log2_fold_change" label
      unless the declared representation actually justifies it.
    - Fail explicitly (UnsupportedMethodError) whenever the configured
      method cannot genuinely execute in the current environment, or
      is incompatible with the declared expression representation.
    - Retain, rather than silently drop, genes for which the method
      cannot produce a valid statistical result (specification
      Section 6.4), recording a human-readable reason when available.

Scope and Version 1 method authorization:
    The frozen specification (Section 5.1) names a negative-binomial,
    count-based framework such as DESeq2 as the preferred method for
    raw RNA-seq count data. No genuine DESeq2-capable backend is
    installed in this environment (no R/Bioconductor bridge such as
    rpy2, and no Python-native equivalent such as pydeseq2 --
    verified absent; requirements.txt does not declare either).
    Configuring StatisticalMethod.DESEQ2 therefore always raises
    UnsupportedMethodError. DESeq2 execution is never faked and never
    silently replaced by another method.

    For already-normalized expression data, specification Section 5.2
    requires "a statistically appropriate method for that
    representation," explicitly justified, or explicit rejection of
    the dataset. This module implements Welch's two-sample t-test
    (two-sided, unequal-variance, via scipy.stats.ttest_ind) as that
    justified method for the NORMALIZED_LOG2 and NORMALIZED_LINEAR
    representations. Welch's t-test is the standard, well-established
    method for comparing two independent groups' means without
    assuming raw-count structure or equal variances, which is exactly
    the situation Section 5.2 describes. It is documented here, and
    recorded in analysis metadata, as the specific justification for
    this method choice (per Section 5.1's documentation requirement,
    applied to whichever method is actually used).

    Welch's t-test is NOT considered compatible with RAW_COUNTS: a
    t-test does not model count dispersion, so applying it to raw
    counts would not be "a statistically appropriate method" for that
    representation. RAW_COUNTS data can therefore currently only be
    analyzed with DESeq2, which is unavailable in this environment --
    RAW_COUNTS analyses fail explicitly rather than silently falling
    back to a t-test.
"""

import warnings

import numpy as np
import pandas as pd
from scipy import __version__ as SCIPY_VERSION
from scipy import stats as scipy_stats

from src.differential_expression.exceptions import (
    StatisticalMethodError,
    UnsupportedMethodError,
)
from src.differential_expression.models import (
    DEAnalysisConfiguration,
    ExpressionRepresentation,
    GeneFilterResult,
    GeneResult,
    GroupAssignment,
    StatisticalMethod,
    ValidatedExpressionMatrix,
)

METHODS_VERSION = "1.0"

# Expression representations each statistical method is compatible
# with, per specification Section 5.1/5.2 and the module docstring
# above. DESeq2 is a count-based, negative-binomial framework and
# must not be applied to already normalized or log-transformed
# values. Welch's t-test is justified only for already-normalized
# representations, not for raw counts.
_COMPATIBLE_REPRESENTATIONS: dict[
    StatisticalMethod, frozenset[ExpressionRepresentation]
] = {
    StatisticalMethod.DESEQ2: frozenset({ExpressionRepresentation.RAW_COUNTS}),
    StatisticalMethod.WELCH_T_TEST: frozenset(
        {
            ExpressionRepresentation.NORMALIZED_LOG2,
            ExpressionRepresentation.NORMALIZED_LINEAR,
        }
    ),
}

# Effect-size label is only "log2_fold_change" when the declared
# representation is already on the log2 scale (specification Section
# 6.3: the result must include the log2 fold change "for methods that
# provide" one). For any other representation, a mean difference is
# reported under a non-fold-change label so it is never mislabeled.
_LOG2_FOLD_CHANGE_LABEL = "log2_fold_change"
_MEAN_DIFFERENCE_LABEL = "mean_difference"

# Mathematical floor for computing a two-sample variance-based
# statistic for one gene: at least two non-missing values are
# required in each group. This mirrors
# models.ABSOLUTE_MINIMUM_REPLICATES_PER_GROUP but is applied per
# gene (a gene may have missing values in some samples even when the
# comparison's overall replication requirement is satisfied).
_MINIMUM_NON_MISSING_PER_GROUP = 2

__all__ = [
    "METHODS_VERSION",
    "check_method_compatibility",
    "execute_statistical_method",
    "get_method_version",
]


def check_method_compatibility(
    statistical_method: StatisticalMethod,
    expression_representation: ExpressionRepresentation,
) -> None:
    """Verify that a configured method is compatible with the
    declared expression representation.

    Args:
        statistical_method: The configured method.
        expression_representation: The declared expression
            representation.

    Raises:
        UnsupportedMethodError: If the method is incompatible with
            the representation.
    """
    compatible_representations = _COMPATIBLE_REPRESENTATIONS.get(
        statistical_method, frozenset()
    )

    if expression_representation not in compatible_representations:
        raise UnsupportedMethodError(
            f"Statistical method '{statistical_method.value}' is not "
            "compatible with expression representation "
            f"'{expression_representation.value}'. Compatible "
            f"representation(s) for this method: "
            f"{sorted(r.value for r in compatible_representations) or 'none'}."
        )


def get_method_version(statistical_method: StatisticalMethod) -> str | None:
    """Return the software version relevant to a statistical method,
    for reproducibility.py to record in analysis metadata.

    Args:
        statistical_method: The configured method.

    Returns:
        The scipy version string for WELCH_T_TEST. None for DESEQ2,
        since no backend is available in this environment to report
        a version for (fabricating one is not permitted).
    """
    if statistical_method == StatisticalMethod.WELCH_T_TEST:
        return f"scipy {SCIPY_VERSION}"
    return None


def _effect_size_label(
    expression_representation: ExpressionRepresentation,
) -> str:
    """Determine the correct effect-size label for a declared
    expression representation, per specification Section 6.3.

    Args:
        expression_representation: The declared expression
            representation.

    Returns:
        "log2_fold_change" only if the representation is already on
        the log2 scale; otherwise "mean_difference".
    """
    if expression_representation == ExpressionRepresentation.NORMALIZED_LOG2:
        return _LOG2_FOLD_CHANGE_LABEL
    return _MEAN_DIFFERENCE_LABEL


def _run_welch_t_test(
    expression: ValidatedExpressionMatrix,
    group_assignment: GroupAssignment,
    gene_filter_result: GeneFilterResult,
    expression_representation: ExpressionRepresentation,
) -> tuple[GeneResult, ...]:
    """Execute Welch's two-sample t-test (two-sided, equal_var=False)
    for every tested gene.

    Args:
        expression: The validated expression matrix.
        group_assignment: The resolved two-group comparison.
        gene_filter_result: Which genes were tested vs. filtered.
        expression_representation: The declared expression
            representation, used only to select the effect-size
            label.

    Returns:
        One GeneResult per gene in gene_filter_result.tested_gene_ids.
        A gene for which fewer than
        _MINIMUM_NON_MISSING_PER_GROUP non-missing values are
        available in either group, or for which scipy reports a
        non-finite result (e.g. near-identical values in both
        groups), is retained with effect_size/raw_p_value set to
        None and missing_reason explaining why -- never silently
        dropped and never a fabricated value (specification Section
        6.4).

        adjusted_p_value and significant are left as None here; they
        are populated by multiple_testing.py, which is the single
        place multiple-testing correction happens.
    """
    dataframe = expression.dataframe
    gene_id_column = expression.gene_id_column
    label = _effect_size_label(expression_representation)

    reference_columns = list(group_assignment.reference_sample_ids)
    comparison_columns = list(group_assignment.comparison_sample_ids)

    indexed = dataframe.set_index(gene_id_column)

    results: list[GeneResult] = []

    for gene_id in gene_filter_result.tested_gene_ids:
        row = indexed.loc[gene_id]

        reference_values = pd.to_numeric(
            row[reference_columns], errors="raise"
        ).dropna().to_numpy(dtype=float)
        comparison_values = pd.to_numeric(
            row[comparison_columns], errors="raise"
        ).dropna().to_numpy(dtype=float)

        if (
            len(reference_values) < _MINIMUM_NON_MISSING_PER_GROUP
            or len(comparison_values) < _MINIMUM_NON_MISSING_PER_GROUP
        ):
            results.append(
                GeneResult(
                    gene_id=gene_id,
                    tested=True,
                    effect_size=None,
                    effect_size_label=None,
                    raw_p_value=None,
                    adjusted_p_value=None,
                    significant=None,
                    missing_reason=(
                        "Fewer than "
                        f"{_MINIMUM_NON_MISSING_PER_GROUP} non-missing "
                        "expression value(s) available for this gene in "
                        "the reference group "
                        f"({len(reference_values)} available) or the "
                        "comparison group "
                        f"({len(comparison_values)} available); Welch's "
                        "t-test cannot be computed."
                    ),
                )
            )
            continue

        effect_size = float(
            np.mean(comparison_values) - np.mean(reference_values)
        )

        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always")
            test_result = scipy_stats.ttest_ind(
                comparison_values,
                reference_values,
                equal_var=False,
                nan_policy="raise",
            )

        raw_p_value = float(test_result.pvalue)

        if not np.isfinite(raw_p_value):
            warning_messages = "; ".join(
                str(w.message) for w in caught_warnings
            ) or "scipy.stats.ttest_ind returned a non-finite p-value."
            results.append(
                GeneResult(
                    gene_id=gene_id,
                    tested=True,
                    effect_size=effect_size,
                    effect_size_label=label,
                    raw_p_value=None,
                    adjusted_p_value=None,
                    significant=None,
                    missing_reason=(
                        "Welch's t-test could not produce a valid "
                        f"p-value for this gene: {warning_messages}"
                    ),
                )
            )
            continue

        results.append(
            GeneResult(
                gene_id=gene_id,
                tested=True,
                effect_size=effect_size,
                effect_size_label=label,
                raw_p_value=raw_p_value,
                adjusted_p_value=None,
                significant=None,
                missing_reason=None,
            )
        )

    return tuple(results)


def execute_statistical_method(
    configuration: DEAnalysisConfiguration,
    expression: ValidatedExpressionMatrix,
    group_assignment: GroupAssignment,
    gene_filter_result: GeneFilterResult,
) -> tuple[GeneResult, ...]:
    """Execute the configured statistical method against every tested
    gene.

    Args:
        configuration: The analysis configuration (statistical_method,
            expression_representation).
        expression: The validated expression matrix.
        group_assignment: The resolved two-group comparison.
        gene_filter_result: Which genes were tested vs. filtered.

    Returns:
        One GeneResult per gene in gene_filter_result.tested_gene_ids
        when the configured method is WELCH_T_TEST and is compatible
        with the declared representation.

    Raises:
        UnsupportedMethodError: If the configured method is
            incompatible with the declared expression representation,
            or if the configured method is DESEQ2 (no genuine backend
            available in this environment -- see module docstring).
        StatisticalMethodError: If the statistical framework reports
            a dataset-level failure that prevents the method from
            running at all (as opposed to an individual gene
            producing a missing result, which is not an error).
    """
    check_method_compatibility(
        configuration.statistical_method,
        configuration.expression_representation,
    )

    if configuration.statistical_method == StatisticalMethod.DESEQ2:
        raise UnsupportedMethodError(
            "Statistical method 'deseq2' is configured and compatible "
            "with the declared 'raw_counts' expression representation, "
            "but no genuine DESeq2-capable backend is available in "
            "this environment (no R/Bioconductor bridge such as rpy2, "
            "and no Python-native equivalent such as pydeseq2 is "
            "installed or declared in requirements.txt). Per "
            "specification Section 5.4, the analysis fails explicitly "
            "rather than substituting a different method or faking "
            "DESeq2 execution."
        )

    if configuration.statistical_method == StatisticalMethod.WELCH_T_TEST:
        try:
            return _run_welch_t_test(
                expression=expression,
                group_assignment=group_assignment,
                gene_filter_result=gene_filter_result,
                expression_representation=configuration.expression_representation,
            )
        except KeyError as error:
            raise StatisticalMethodError(
                "Welch's t-test failed at the dataset level: a tested "
                f"gene identifier could not be located in the "
                f"expression matrix index ({error})."
            ) from error

    # Unreachable while StatisticalMethod has exactly the two members
    # handled above. This branch exists so that adding a future
    # authorized method to the enum without also adding its execution
    # branch here fails loudly rather than silently falling through.
    raise UnsupportedMethodError(
        "No execution path is implemented for statistical method "
        f"'{configuration.statistical_method.value}'."
    )