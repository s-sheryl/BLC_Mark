"""Tests for src.differential_expression.qc."""

from pathlib import Path

import pytest

from src.differential_expression.configuration import build_configuration
from src.differential_expression.exceptions import QualityControlError
from src.differential_expression.models import (
    ExcludedSample,
    GeneFilterResult,
    GeneResult,
    GroupAssignment,
    SampleMatchResult,
)
from src.differential_expression.qc import build_qc_report, verify_post_analysis_quality


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
        significance_threshold=0.05,
    )
    kwargs.update(overrides)
    return build_configuration(**kwargs)


def _sample_match():
    return SampleMatchResult(
        matched_samples=("S1", "S2", "S3", "S4"),
        expression_only_samples=("S5",),
        metadata_only_samples=("S6",),
    )


def _group_assignment():
    return GroupAssignment(
        reference_group="normal",
        comparison_group="tumor",
        reference_sample_ids=("S1", "S2"),
        comparison_sample_ids=("S3", "S4"),
        excluded_samples=(
            ExcludedSample(sample_id="S5", reason="expression-only", stage="validation"),
            ExcludedSample(sample_id="S6", reason="metadata-only", stage="validation"),
        ),
    )


def _gene_filter_result():
    return GeneFilterResult(
        tested_gene_ids=("A", "B"),
        filtered_gene_ids=("C",),
        criterion_description="No filtering.",
        input_gene_count=3,
    )


def test_build_qc_report_counts():
    config = _config()
    report = build_qc_report(config, _sample_match(), _group_assignment(), _gene_filter_result())

    assert report.input_gene_count == 3
    assert report.tested_gene_count == 2
    assert report.filtered_gene_count == 1
    assert report.initial_sample_count == 6  # 4 matched + 1 expr-only + 1 meta-only
    assert report.included_sample_count == 4
    assert report.excluded_sample_count == 2
    assert report.reference_group_size == 2
    assert report.comparison_group_size == 2


def test_build_qc_report_distinguishes_observations_from_configuration():
    config = _config(significance_threshold=0.1)
    report = build_qc_report(config, _sample_match(), _group_assignment(), _gene_filter_result())

    # Configuration decisions must reflect what was configured...
    assert report.significance_threshold == 0.1
    assert report.statistical_method == config.statistical_method
    # ...while observed facts must reflect what was actually found.
    assert report.tested_gene_count == 2


def _gene_result(gene_id, raw_p=0.01, adj_p=0.02, significant=True, threshold=0.05):
    return GeneResult(
        gene_id=gene_id,
        tested=True,
        effect_size=1.0,
        effect_size_label="log2_fold_change",
        raw_p_value=raw_p,
        adjusted_p_value=adj_p,
        significant=significant,
        missing_reason=None,
    )


def test_post_analysis_qc_passes_for_valid_results():
    config = _config()
    results = (
        _gene_result("A", 0.01, 0.02, True),
        _gene_result("B", 0.5, 0.5, False),
    )
    verify_post_analysis_quality(results, _gene_filter_result(), config)


def test_post_analysis_qc_fails_on_missing_gene():
    config = _config()
    results = (_gene_result("A"),)  # missing "B"
    with pytest.raises(QualityControlError):
        verify_post_analysis_quality(results, _gene_filter_result(), config)


def test_post_analysis_qc_fails_on_unexpected_gene():
    config = _config()
    results = (_gene_result("A"), _gene_result("B"), _gene_result("Z"))
    with pytest.raises(QualityControlError):
        verify_post_analysis_quality(results, _gene_filter_result(), config)


def test_post_analysis_qc_fails_on_duplicate_gene():
    config = _config()
    results = (_gene_result("A"), _gene_result("A"))
    gfr = GeneFilterResult(
        tested_gene_ids=("A",),
        filtered_gene_ids=(),
        criterion_description="No filtering.",
        input_gene_count=1,
    )
    with pytest.raises(QualityControlError):
        verify_post_analysis_quality(results, gfr, config)


def test_post_analysis_qc_fails_on_out_of_range_p_value():
    config = _config()
    results = (
        GeneResult(
            gene_id="A", tested=True, effect_size=1.0, effect_size_label="log2_fold_change",
            raw_p_value=1.5, adjusted_p_value=None, significant=None, missing_reason=None,
        ),
        _gene_result("B"),
    )
    with pytest.raises(QualityControlError):
        verify_post_analysis_quality(results, _gene_filter_result(), config)


def test_post_analysis_qc_fails_on_adjusted_without_raw():
    config = _config()
    results = (
        GeneResult(
            gene_id="A", tested=True, effect_size=1.0, effect_size_label="log2_fold_change",
            raw_p_value=None, adjusted_p_value=0.02, significant=True, missing_reason=None,
        ),
        _gene_result("B"),
    )
    with pytest.raises(QualityControlError):
        verify_post_analysis_quality(results, _gene_filter_result(), config)


def test_post_analysis_qc_fails_on_inconsistent_significance():
    config = _config(significance_threshold=0.05)
    results = (
        GeneResult(
            gene_id="A", tested=True, effect_size=1.0, effect_size_label="log2_fold_change",
            raw_p_value=0.01, adjusted_p_value=0.02, significant=False,  # wrong: 0.02 < 0.05
            missing_reason=None,
        ),
        _gene_result("B"),
    )
    with pytest.raises(QualityControlError):
        verify_post_analysis_quality(results, _gene_filter_result(), config)


def test_post_analysis_qc_passes_with_missing_gene_results():
    config = _config()
    results = (
        GeneResult(
            gene_id="A", tested=True, effect_size=None, effect_size_label=None,
            raw_p_value=None, adjusted_p_value=None, significant=None,
            missing_reason="Insufficient values.",
        ),
        _gene_result("B"),
    )
    verify_post_analysis_quality(results, _gene_filter_result(), config)