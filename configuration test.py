"""Tests for src.differential_expression.configuration."""

from pathlib import Path

import pytest

from src.dataset_registry import CancerType
from src.differential_expression.configuration import build_configuration
from src.differential_expression.exceptions import InvalidConfigurationError
from src.differential_expression.models import (
    ExpressionRepresentation,
    StatisticalMethod,
)


def _valid_kwargs(**overrides):
    kwargs = dict(
        analysis_id="analysis-001",
        cancer_cohort="TCGA-BRCA",
        expression_matrix_path="/tmp/expression.csv",
        metadata_path="/tmp/metadata.csv",
        gene_id_column="gene_id",
        sample_id_column="sample_id",
        group_column="group",
        reference_group="normal",
        comparison_group="tumor",
        expression_representation="normalized_log2",
        statistical_method="deseq2",
        minimum_replicates_per_group=3,
        output_dir="/tmp/output",
    )
    kwargs.update(overrides)
    return kwargs


def test_build_configuration_from_strings():
    config = build_configuration(**_valid_kwargs())
    assert config.cancer_cohort == CancerType.BRCA
    assert config.expression_representation == ExpressionRepresentation.NORMALIZED_LOG2
    assert config.statistical_method == StatisticalMethod.DESEQ2
    assert config.expression_matrix_path == Path("/tmp/expression.csv")


def test_build_configuration_from_typed_values():
    config = build_configuration(
        **_valid_kwargs(
            cancer_cohort=CancerType.LUAD,
            expression_matrix_path=Path("/tmp/e.csv"),
            metadata_path=Path("/tmp/m.csv"),
            expression_representation=ExpressionRepresentation.RAW_COUNTS,
            statistical_method=StatisticalMethod.DESEQ2,
            output_dir=Path("/tmp/out"),
        )
    )
    assert config.cancer_cohort == CancerType.LUAD
    assert config.expression_representation == ExpressionRepresentation.RAW_COUNTS


def test_invalid_cancer_cohort_string_rejected():
    with pytest.raises(InvalidConfigurationError):
        build_configuration(**_valid_kwargs(cancer_cohort="TCGA-PANCREAS"))


def test_invalid_statistical_method_string_rejected():
    with pytest.raises(InvalidConfigurationError):
        build_configuration(**_valid_kwargs(statistical_method="ancova"))


def test_welch_t_test_is_a_valid_statistical_method_string():
    # Version 1 authorizes both DESEQ2 (Section 5.1, unavailable in
    # this environment -- see methods.py) and WELCH_T_TEST (the
    # justified method for normalized data under Section 5.2). Both
    # must be accepted as valid configuration values.
    config = build_configuration(**_valid_kwargs(statistical_method="welch_t_test"))
    assert config.statistical_method.value == "welch_t_test"


def test_invalid_expression_representation_string_rejected():
    with pytest.raises(InvalidConfigurationError):
        build_configuration(**_valid_kwargs(expression_representation="raw"))


def test_empty_path_string_rejected():
    with pytest.raises(InvalidConfigurationError):
        build_configuration(**_valid_kwargs(expression_matrix_path="   "))


def test_wrong_path_type_rejected():
    with pytest.raises(InvalidConfigurationError):
        build_configuration(**_valid_kwargs(expression_matrix_path=12345))


def test_identical_reference_and_comparison_group_wrapped_as_configuration_error():
    with pytest.raises(InvalidConfigurationError):
        build_configuration(**_valid_kwargs(comparison_group="normal"))


def test_replication_below_floor_wrapped_as_configuration_error():
    with pytest.raises(InvalidConfigurationError):
        build_configuration(**_valid_kwargs(minimum_replicates_per_group=1))


def test_invalid_significance_threshold_wrapped_as_configuration_error():
    with pytest.raises(InvalidConfigurationError):
        build_configuration(**_valid_kwargs(significance_threshold=2.0))