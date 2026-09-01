"""Tests for src.differential_expression.comparison."""

from pathlib import Path

import pandas as pd
import pytest

from src.differential_expression.comparison import resolve_comparison
from src.differential_expression.configuration import build_configuration
from src.differential_expression.exceptions import (
    InsufficientReplicationError,
    InvalidConfigurationError,
    SampleMismatchError,
)
from src.differential_expression.models import SampleMatchResult, ValidatedMetadata


def _config(**overrides):
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
        minimum_replicates_per_group=2,
        output_dir="/tmp/output",
    )
    kwargs.update(overrides)
    return build_configuration(**kwargs)


def _metadata(rows, group_column="group"):
    df = pd.DataFrame(rows)
    return ValidatedMetadata(
        file_path=Path("/tmp/metadata.csv"),
        sample_id_column="sample_id",
        group_column=group_column,
        dataframe=df,
    )


def test_valid_two_group_comparison():
    metadata = _metadata(
        {
            "sample_id": ["S1", "S2", "S3", "S4"],
            "group": ["normal", "normal", "tumor", "tumor"],
        }
    )
    match = SampleMatchResult(
        matched_samples=("S1", "S2", "S3", "S4"),
        expression_only_samples=(),
        metadata_only_samples=(),
    )
    assignment = resolve_comparison(_config(), metadata, match)
    assert assignment.reference_sample_ids == ("S1", "S2")
    assert assignment.comparison_sample_ids == ("S3", "S4")
    assert assignment.excluded_samples == ()


def test_missing_reference_group_raises():
    metadata = _metadata(
        {"sample_id": ["S1", "S2"], "group": ["tumor", "tumor"]}
    )
    match = SampleMatchResult(
        matched_samples=("S1", "S2"),
        expression_only_samples=(),
        metadata_only_samples=(),
    )
    with pytest.raises(SampleMismatchError):
        resolve_comparison(_config(), metadata, match)


def test_missing_comparison_group_raises():
    metadata = _metadata(
        {"sample_id": ["S1", "S2"], "group": ["normal", "normal"]}
    )
    match = SampleMatchResult(
        matched_samples=("S1", "S2"),
        expression_only_samples=(),
        metadata_only_samples=(),
    )
    with pytest.raises(SampleMismatchError):
        resolve_comparison(_config(), metadata, match)


def test_insufficient_replication_raises():
    metadata = _metadata(
        {
            "sample_id": ["S1", "S2", "S3"],
            "group": ["normal", "tumor", "tumor"],
        }
    )
    match = SampleMatchResult(
        matched_samples=("S1", "S2", "S3"),
        expression_only_samples=(),
        metadata_only_samples=(),
    )
    with pytest.raises(InsufficientReplicationError):
        resolve_comparison(_config(minimum_replicates_per_group=2), metadata, match)


def test_group_column_mismatch_raises_configuration_error():
    metadata = _metadata(
        {"sample_id": ["S1", "S2"], "cohort": ["normal", "tumor"]},
        group_column="cohort",
    )
    match = SampleMatchResult(
        matched_samples=("S1", "S2"),
        expression_only_samples=(),
        metadata_only_samples=(),
    )
    with pytest.raises(InvalidConfigurationError):
        resolve_comparison(_config(), metadata, match)


def test_deterministic_ordering():
    metadata = _metadata(
        {
            "sample_id": ["S4", "S1", "S3", "S2"],
            "group": ["tumor", "normal", "tumor", "normal"],
        }
    )
    match = SampleMatchResult(
        matched_samples=("S4", "S1", "S3", "S2"),
        expression_only_samples=(),
        metadata_only_samples=(),
    )
    assignment = resolve_comparison(_config(), metadata, match)
    assert assignment.reference_sample_ids == ("S1", "S2")
    assert assignment.comparison_sample_ids == ("S3", "S4")


def test_third_group_and_mismatched_samples_are_recorded_as_excluded():
    metadata = _metadata(
        {
            "sample_id": ["S1", "S2", "S3", "S4", "S5"],
            "group": ["normal", "normal", "tumor", "tumor", "metastatic"],
        }
    )
    match = SampleMatchResult(
        matched_samples=("S1", "S2", "S3", "S4", "S5"),
        expression_only_samples=("S6",),
        metadata_only_samples=("S7",),
    )
    assignment = resolve_comparison(_config(), metadata, match)
    excluded_ids = {excluded.sample_id for excluded in assignment.excluded_samples}
    assert excluded_ids == {"S5", "S6", "S7"}

    stages = {excluded.sample_id: excluded.stage for excluded in assignment.excluded_samples}
    assert stages["S5"] == "comparison"
    assert stages["S6"] == "validation"
    assert stages["S7"] == "validation"