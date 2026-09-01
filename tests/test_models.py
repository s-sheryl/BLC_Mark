"""Tests for src.differential_expression.models."""

from datetime import datetime, timezone
from pathlib import Path

import pytest

from src.dataset_registry import CancerType
from src.differential_expression.models import (
    ABSOLUTE_MINIMUM_REPLICATES_PER_GROUP,
    AnalysisResult,
    AnalysisStatus,
    ComparisonDefinition,
    DEAnalysisConfiguration,
    ExcludedSample,
    ExpressionRepresentation,
    FailureInfo,
    GeneFilterConfiguration,
    GroupAssignment,
    MultipleTestingMethod,
    QCReport,
    StatisticalMethod,
)


def _valid_config_kwargs(**overrides):
    kwargs = dict(
        analysis_id="analysis-001",
        cancer_cohort=CancerType.BRCA,
        expression_matrix_path=Path("/tmp/expression.csv"),
        metadata_path=Path("/tmp/metadata.csv"),
        gene_id_column="gene_id",
        sample_id_column="sample_id",
        group_column="group",
        reference_group="normal",
        comparison_group="tumor",
        expression_representation=ExpressionRepresentation.NORMALIZED_LOG2,
        statistical_method=StatisticalMethod.DESEQ2,
        minimum_replicates_per_group=3,
        output_dir=Path("/tmp/output"),
    )
    kwargs.update(overrides)
    return kwargs


def test_valid_configuration_constructs():
    config = DEAnalysisConfiguration(**_valid_config_kwargs())
    assert config.reference_group == "normal"
    assert config.multiple_testing_method == MultipleTestingMethod.BENJAMINI_HOCHBERG
    assert config.significance_threshold == 0.05
    assert config.effect_size_threshold is None


def test_identical_reference_and_comparison_group_rejected():
    with pytest.raises(ValueError):
        DEAnalysisConfiguration(**_valid_config_kwargs(comparison_group="normal"))


def test_invalid_significance_threshold_rejected():
    with pytest.raises(ValueError):
        DEAnalysisConfiguration(**_valid_config_kwargs(significance_threshold=1.5))

    with pytest.raises(ValueError):
        DEAnalysisConfiguration(**_valid_config_kwargs(significance_threshold=0.0))


def test_minimum_replicates_below_floor_rejected():
    with pytest.raises(ValueError):
        DEAnalysisConfiguration(
            **_valid_config_kwargs(
                minimum_replicates_per_group=ABSOLUTE_MINIMUM_REPLICATES_PER_GROUP - 1
            )
        )


def test_minimum_replicates_at_floor_accepted():
    config = DEAnalysisConfiguration(
        **_valid_config_kwargs(
            minimum_replicates_per_group=ABSOLUTE_MINIMUM_REPLICATES_PER_GROUP
        )
    )
    assert config.minimum_replicates_per_group == ABSOLUTE_MINIMUM_REPLICATES_PER_GROUP


def test_negative_effect_size_threshold_rejected():
    with pytest.raises(ValueError):
        DEAnalysisConfiguration(**_valid_config_kwargs(effect_size_threshold=-1.0))


def test_wrong_type_cancer_cohort_rejected():
    with pytest.raises(TypeError):
        DEAnalysisConfiguration(**_valid_config_kwargs(cancer_cohort="TCGA-BRCA"))


def test_gene_filter_configuration_defaults_to_no_filtering():
    filter_config = GeneFilterConfiguration()
    assert filter_config.apply_filter is False
    assert filter_config.minimum_mean_expression is None


def test_gene_filter_configuration_requires_threshold_when_enabled():
    with pytest.raises(ValueError):
        GeneFilterConfiguration(apply_filter=True)


def test_gene_filter_configuration_rejects_threshold_when_disabled():
    with pytest.raises(ValueError):
        GeneFilterConfiguration(apply_filter=False, minimum_mean_expression=1.0)


def test_gene_filter_configuration_valid_enabled():
    filter_config = GeneFilterConfiguration(
        apply_filter=True, minimum_mean_expression=1.0
    )
    assert filter_config.apply_filter is True
    assert filter_config.minimum_mean_expression == 1.0


def test_comparison_definition_rejects_identical_groups():
    with pytest.raises(ValueError):
        ComparisonDefinition(
            group_column="group", reference_group="a", comparison_group="a"
        )


def test_group_assignment_rejects_overlapping_samples():
    with pytest.raises(ValueError):
        GroupAssignment(
            reference_group="normal",
            comparison_group="tumor",
            reference_sample_ids=("s1", "s2"),
            comparison_sample_ids=("s2", "s3"),
        )


def test_group_assignment_valid():
    assignment = GroupAssignment(
        reference_group="normal",
        comparison_group="tumor",
        reference_sample_ids=("s1", "s2"),
        comparison_sample_ids=("s3", "s4"),
        excluded_samples=(ExcludedSample("s5", "missing group label", "validation"),),
    )
    assert assignment.reference_sample_ids == ("s1", "s2")
    assert len(assignment.excluded_samples) == 1


def test_analysis_result_succeeded_requires_full_payload():
    with pytest.raises(ValueError):
        AnalysisResult(analysis_id="a1", status=AnalysisStatus.SUCCEEDED)


def test_analysis_result_failed_requires_failure_info():
    with pytest.raises(ValueError):
        AnalysisResult(analysis_id="a1", status=AnalysisStatus.FAILED)


def test_analysis_result_failed_with_failure_info_is_valid():
    result = AnalysisResult(
        analysis_id="a1",
        status=AnalysisStatus.FAILED,
        failure=FailureInfo(
            stage="validation", category="InvalidExpressionMatrixError", message="boom"
        ),
    )
    assert result.failure.stage == "validation"


def test_analysis_result_succeeded_with_full_payload_is_valid():
    qc = QCReport(
        input_gene_count=10,
        tested_gene_count=10,
        filtered_gene_count=0,
        initial_sample_count=6,
        included_sample_count=6,
        excluded_sample_count=0,
        excluded_samples=(),
        reference_group_size=3,
        comparison_group_size=3,
        gene_filter_criterion="No gene filtering configured.",
        statistical_method=StatisticalMethod.DESEQ2,
        multiple_testing_method=MultipleTestingMethod.BENJAMINI_HOCHBERG,
        significance_threshold=0.05,
    )
    result = AnalysisResult(
        analysis_id="a1",
        status=AnalysisStatus.SUCCEEDED,
        gene_results=(),
        qc_report=qc,
        metadata=None if False else _minimal_metadata(),
    )
    assert result.status == AnalysisStatus.SUCCEEDED


def _minimal_metadata():
    from src.differential_expression.models import AnalysisMetadata

    return AnalysisMetadata(
        analysis_id="a1",
        timestamp_utc=datetime.now(timezone.utc),
        cancer_cohort=CancerType.BRCA,
        expression_matrix_path=Path("/tmp/expression.csv"),
        metadata_path=Path("/tmp/metadata.csv"),
        expression_matrix_sha256="0" * 64,
        metadata_sha256="0" * 64,
        expression_representation=ExpressionRepresentation.NORMALIZED_LOG2,
        statistical_method=StatisticalMethod.DESEQ2,
        statistical_method_version=None,
        design="two-group comparison",
        reference_group="normal",
        comparison_group="tumor",
        included_sample_ids=("s1", "s2", "s3", "s4", "s5", "s6"),
        reference_group_size=3,
        comparison_group_size=3,
        gene_filter_criterion="No gene filtering configured.",
        multiple_testing_method=MultipleTestingMethod.BENJAMINI_HOCHBERG,
        significance_threshold=0.05,
        effect_size_threshold=None,
        python_version="3.12.3",
        package_versions={"pandas": "2.0.0"},
        de_package_version="1.0",
        blc_mark_version=None,
        analysis_status=AnalysisStatus.SUCCEEDED,
    )