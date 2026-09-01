"""Tests for src.differential_expression.methods.

Version 1 authorizes two statistical methods: DESEQ2 (specification
Section 5.1), which always fails explicitly because no genuine
DESeq2-capable backend is installed in this environment, and
WELCH_T_TEST, the justified method (Section 5.2) for already-normalized
expression data, executed via scipy.stats.ttest_ind. These tests verify
both the explicit-failure behavior for DESeq2/incompatible
configurations and genuine, correct execution for Welch's t-test,
including missing-value and degenerate-gene handling.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.differential_expression.configuration import build_configuration
from src.differential_expression.exceptions import UnsupportedMethodError
from src.differential_expression.methods import (
    check_method_compatibility,
    execute_statistical_method,
    get_method_version,
)
from src.differential_expression.models import (
    ExpressionRepresentation,
    GeneFilterResult,
    GroupAssignment,
    StatisticalMethod,
    ValidatedExpressionMatrix,
)


def _expression_matrix(extra_rows=None):
    data = {
        "gene_id": ["A", "B"],
        "S1": [10.0, 12.0],
        "S2": [11.0, 13.0],
        "S3": [50.0, 55.0],
        "S4": [52.0, 58.0],
    }
    df = pd.DataFrame(data)
    if extra_rows is not None:
        df = pd.concat([df, pd.DataFrame(extra_rows)], ignore_index=True)
    return ValidatedExpressionMatrix(
        file_path=Path("/tmp/expression.csv"),
        gene_id_column="gene_id",
        sample_columns=("S1", "S2", "S3", "S4"),
        gene_ids=tuple(df["gene_id"]),
        dataframe=df,
    )


def _group_assignment():
    return GroupAssignment(
        reference_group="normal",
        comparison_group="tumor",
        reference_sample_ids=("S1", "S2"),
        comparison_sample_ids=("S3", "S4"),
    )


def _gene_filter_result(tested_gene_ids=("A", "B")):
    return GeneFilterResult(
        tested_gene_ids=tuple(tested_gene_ids),
        filtered_gene_ids=(),
        criterion_description="No gene filtering configured.",
        input_gene_count=len(tested_gene_ids),
    )


def _config(**overrides):
    kwargs = dict(
        analysis_id="a1",
        cancer_cohort="TCGA-BRCA",
        expression_matrix_path="/tmp/expression.csv",
        metadata_path="/tmp/metadata.csv",
        gene_id_column="gene_id",
        sample_id_column="sample_id",
        group_column="group",
        reference_group="normal",
        comparison_group="tumor",
        expression_representation="normalized_log2",
        statistical_method="welch_t_test",
        minimum_replicates_per_group=2,
        output_dir="/tmp/output",
    )
    kwargs.update(overrides)
    return build_configuration(**kwargs)


# --- Compatibility checks ---


def test_deseq2_compatible_with_raw_counts_representation():
    check_method_compatibility(
        StatisticalMethod.DESEQ2, ExpressionRepresentation.RAW_COUNTS
    )


@pytest.mark.parametrize(
    "representation",
    [ExpressionRepresentation.NORMALIZED_LOG2, ExpressionRepresentation.NORMALIZED_LINEAR],
)
def test_deseq2_incompatible_with_normalized_representations(representation):
    with pytest.raises(UnsupportedMethodError):
        check_method_compatibility(StatisticalMethod.DESEQ2, representation)


@pytest.mark.parametrize(
    "representation",
    [ExpressionRepresentation.NORMALIZED_LOG2, ExpressionRepresentation.NORMALIZED_LINEAR],
)
def test_welch_t_test_compatible_with_normalized_representations(representation):
    check_method_compatibility(StatisticalMethod.WELCH_T_TEST, representation)


def test_welch_t_test_incompatible_with_raw_counts():
    with pytest.raises(UnsupportedMethodError):
        check_method_compatibility(
            StatisticalMethod.WELCH_T_TEST, ExpressionRepresentation.RAW_COUNTS
        )


# --- DESeq2 explicit unavailability (never faked, never substituted) ---


def test_execute_raises_for_raw_counts_due_to_missing_backend():
    config = _config(expression_representation="raw_counts", statistical_method="deseq2")
    with pytest.raises(UnsupportedMethodError, match="backend"):
        execute_statistical_method(
            config, _expression_matrix(), _group_assignment(), _gene_filter_result()
        )


def test_execute_raises_for_normalized_deseq2_due_to_incompatibility():
    config = _config(expression_representation="normalized_log2", statistical_method="deseq2")
    with pytest.raises(UnsupportedMethodError, match="not compatible"):
        execute_statistical_method(
            config, _expression_matrix(), _group_assignment(), _gene_filter_result()
        )


def test_no_silent_fallback_deseq2_always_raises_same_exception_type():
    for representation in ExpressionRepresentation:
        config = _config(
            expression_representation=representation.value,
            statistical_method="deseq2",
        )
        with pytest.raises(UnsupportedMethodError):
            execute_statistical_method(
                config, _expression_matrix(), _group_assignment(), _gene_filter_result()
            )


def test_get_method_version_deseq2_is_none_not_fabricated():
    assert get_method_version(StatisticalMethod.DESEQ2) is None


# --- Welch's t-test genuine execution ---


def test_welch_t_test_executes_and_returns_one_result_per_tested_gene():
    config = _config()
    results = execute_statistical_method(
        config, _expression_matrix(), _group_assignment(), _gene_filter_result()
    )
    assert len(results) == 2
    assert {result.gene_id for result in results} == {"A", "B"}


def test_welch_t_test_effect_size_direction_is_comparison_minus_reference():
    # comparison group (S3, S4) has higher values than reference (S1, S2)
    # for both genes, so effect size must be positive.
    config = _config()
    results = execute_statistical_method(
        config, _expression_matrix(), _group_assignment(), _gene_filter_result()
    )
    for result in results:
        assert result.effect_size is not None
        assert result.effect_size > 0


def test_welch_t_test_effect_size_label_log2_fold_change_when_representation_is_log2():
    config = _config(expression_representation="normalized_log2")
    results = execute_statistical_method(
        config, _expression_matrix(), _group_assignment(), _gene_filter_result()
    )
    for result in results:
        assert result.effect_size_label == "log2_fold_change"


def test_welch_t_test_effect_size_label_mean_difference_when_representation_is_linear():
    config = _config(expression_representation="normalized_linear")
    results = execute_statistical_method(
        config, _expression_matrix(), _group_assignment(), _gene_filter_result()
    )
    for result in results:
        assert result.effect_size_label == "mean_difference"


def test_welch_t_test_raw_p_value_is_populated_and_finite():
    config = _config()
    results = execute_statistical_method(
        config, _expression_matrix(), _group_assignment(), _gene_filter_result()
    )
    for result in results:
        assert result.raw_p_value is not None
        assert np.isfinite(result.raw_p_value)
        assert 0.0 <= result.raw_p_value <= 1.0


def test_welch_t_test_does_not_populate_adjusted_p_value_or_significance():
    # adjusted_p_value/significant belong to multiple_testing.py, not
    # methods.py -- methods.py must leave them None.
    config = _config()
    results = execute_statistical_method(
        config, _expression_matrix(), _group_assignment(), _gene_filter_result()
    )
    for result in results:
        assert result.adjusted_p_value is None
        assert result.significant is None


def test_welch_t_test_matches_scipy_reference_value():
    from scipy import stats as scipy_stats

    config = _config()
    results = execute_statistical_method(
        config, _expression_matrix(), _group_assignment(), _gene_filter_result()
    )
    result_a = next(r for r in results if r.gene_id == "A")

    expected = scipy_stats.ttest_ind(
        [50.0, 52.0], [10.0, 11.0], equal_var=False
    )
    assert result_a.raw_p_value == pytest.approx(float(expected.pvalue))
    assert result_a.effect_size == pytest.approx(51.0 - 10.5)


def test_get_method_version_welch_t_test_reports_scipy():
    version = get_method_version(StatisticalMethod.WELCH_T_TEST)
    assert version is not None
    assert "scipy" in version


# --- Missing values / degenerate genes retained, never dropped ---


def test_gene_with_insufficient_non_missing_values_is_retained_with_missing_reason():
    extra = {"gene_id": ["C"], "S1": [5.0], "S2": [np.nan], "S3": [9.0], "S4": [9.5]}
    matrix = _expression_matrix(extra_rows=extra)
    config = _config()
    results = execute_statistical_method(
        config, matrix, _group_assignment(), _gene_filter_result(("A", "B", "C"))
    )
    gene_c = next(r for r in results if r.gene_id == "C")
    assert gene_c.tested is True
    assert gene_c.raw_p_value is None
    assert gene_c.effect_size is None
    assert gene_c.missing_reason is not None
    assert len(results) == 3  # never silently dropped


def test_gene_with_identical_values_in_both_groups_is_retained_with_missing_reason():
    # Near-identical values across both groups make scipy's t-test
    # numerically degenerate (nan p-value); this must be recorded as
    # missing, not silently coerced into a fake result.
    extra = {"gene_id": ["D"], "S1": [7.0], "S2": [7.0], "S3": [7.0], "S4": [7.0]}
    matrix = _expression_matrix(extra_rows=extra)
    config = _config()
    results = execute_statistical_method(
        config, matrix, _group_assignment(), _gene_filter_result(("A", "B", "D"))
    )
    gene_d = next(r for r in results if r.gene_id == "D")
    assert gene_d.tested is True
    assert gene_d.raw_p_value is None
    assert gene_d.missing_reason is not None
    assert len(results) == 3