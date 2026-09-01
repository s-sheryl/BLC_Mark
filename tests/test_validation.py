"""Tests for src.differential_expression.validation."""

from pathlib import Path

import pandas as pd
import pytest

from src.differential_expression.exceptions import (
    InvalidExpressionMatrixError,
    InvalidMetadataError,
    SampleMismatchError,
)
from src.differential_expression.validation import (
    match_samples,
    validate_expression_matrix,
    validate_metadata,
)


@pytest.fixture
def valid_expression_csv(tmp_path: Path) -> Path:
    df = pd.DataFrame(
        {
            "gene_id": ["GENE_A", "GENE_B", "GENE_C"],
            "S1": [1.0, 2.0, 3.0],
            "S2": [1.5, 2.5, 3.5],
            "S3": [10.0, 20.0, 30.0],
            "S4": [11.0, 21.0, 31.0],
        }
    )
    path = tmp_path / "expression.csv"
    df.to_csv(path, index=False)
    return path


@pytest.fixture
def valid_metadata_csv(tmp_path: Path) -> Path:
    df = pd.DataFrame(
        {
            "sample_id": ["S1", "S2", "S3", "S4"],
            "group": ["normal", "normal", "tumor", "tumor"],
        }
    )
    path = tmp_path / "metadata.csv"
    df.to_csv(path, index=False)
    return path


# --- Expression matrix validation ---


def test_valid_expression_matrix(valid_expression_csv):
    result = validate_expression_matrix(valid_expression_csv, gene_id_column="gene_id")
    assert result.gene_id_column == "gene_id"
    assert result.sample_columns == ("S1", "S2", "S3", "S4")
    assert result.gene_ids == ("GENE_A", "GENE_B", "GENE_C")


def test_missing_expression_file(tmp_path):
    with pytest.raises(InvalidExpressionMatrixError):
        validate_expression_matrix(tmp_path / "does_not_exist.csv", gene_id_column="gene_id")


def test_missing_gene_id_column(tmp_path):
    df = pd.DataFrame({"not_gene_id": ["A"], "S1": [1.0]})
    path = tmp_path / "expression.csv"
    df.to_csv(path, index=False)
    with pytest.raises(InvalidExpressionMatrixError):
        validate_expression_matrix(path, gene_id_column="gene_id")


def test_no_sample_columns(tmp_path):
    df = pd.DataFrame({"gene_id": ["A", "B"]})
    path = tmp_path / "expression.csv"
    df.to_csv(path, index=False)
    with pytest.raises(InvalidExpressionMatrixError):
        validate_expression_matrix(path, gene_id_column="gene_id")


def test_duplicate_gene_ids(tmp_path):
    df = pd.DataFrame({"gene_id": ["A", "A"], "S1": [1.0, 2.0]})
    path = tmp_path / "expression.csv"
    df.to_csv(path, index=False)
    with pytest.raises(InvalidExpressionMatrixError):
        validate_expression_matrix(path, gene_id_column="gene_id")


def test_non_numeric_expression_values(tmp_path):
    df = pd.DataFrame({"gene_id": ["A", "B"], "S1": ["not_a_number", "2.0"]})
    path = tmp_path / "expression.csv"
    df.to_csv(path, index=False)
    with pytest.raises(InvalidExpressionMatrixError):
        validate_expression_matrix(path, gene_id_column="gene_id")


def test_missing_gene_identifier_value(tmp_path):
    df = pd.DataFrame({"gene_id": ["A", None], "S1": [1.0, 2.0]})
    path = tmp_path / "expression.csv"
    df.to_csv(path, index=False)
    with pytest.raises(InvalidExpressionMatrixError):
        validate_expression_matrix(path, gene_id_column="gene_id")


# --- Metadata validation ---


def test_valid_metadata(valid_metadata_csv):
    result = validate_metadata(
        valid_metadata_csv, sample_id_column="sample_id", group_column="group"
    )
    assert result.sample_id_column == "sample_id"
    assert list(result.dataframe["sample_id"]) == ["S1", "S2", "S3", "S4"]


def test_missing_metadata_file(tmp_path):
    with pytest.raises(InvalidMetadataError):
        validate_metadata(
            tmp_path / "missing.csv", sample_id_column="sample_id", group_column="group"
        )


def test_missing_sample_id_column(tmp_path):
    df = pd.DataFrame({"group": ["normal"]})
    path = tmp_path / "metadata.csv"
    df.to_csv(path, index=False)
    with pytest.raises(InvalidMetadataError):
        validate_metadata(path, sample_id_column="sample_id", group_column="group")


def test_missing_group_column(tmp_path):
    df = pd.DataFrame({"sample_id": ["S1"]})
    path = tmp_path / "metadata.csv"
    df.to_csv(path, index=False)
    with pytest.raises(InvalidMetadataError):
        validate_metadata(path, sample_id_column="sample_id", group_column="group")


def test_duplicate_sample_ids_in_metadata(tmp_path):
    df = pd.DataFrame({"sample_id": ["S1", "S1"], "group": ["normal", "tumor"]})
    path = tmp_path / "metadata.csv"
    df.to_csv(path, index=False)
    with pytest.raises(InvalidMetadataError):
        validate_metadata(path, sample_id_column="sample_id", group_column="group")


def test_missing_group_label(tmp_path):
    df = pd.DataFrame({"sample_id": ["S1", "S2"], "group": ["normal", None]})
    path = tmp_path / "metadata.csv"
    df.to_csv(path, index=False)
    with pytest.raises(InvalidMetadataError):
        validate_metadata(path, sample_id_column="sample_id", group_column="group")


def test_empty_string_group_label(tmp_path):
    df = pd.DataFrame({"sample_id": ["S1", "S2"], "group": ["normal", "   "]})
    path = tmp_path / "metadata.csv"
    df.to_csv(path, index=False)
    with pytest.raises(InvalidMetadataError):
        validate_metadata(path, sample_id_column="sample_id", group_column="group")


# --- Sample matching ---


def test_complete_sample_overlap(valid_expression_csv, valid_metadata_csv):
    expression = validate_expression_matrix(valid_expression_csv, gene_id_column="gene_id")
    metadata = validate_metadata(
        valid_metadata_csv, sample_id_column="sample_id", group_column="group"
    )
    result = match_samples(expression, metadata)
    assert result.matched_samples == ("S1", "S2", "S3", "S4")
    assert result.expression_only_samples == ()
    assert result.metadata_only_samples == ()


def test_partial_sample_overlap(tmp_path):
    expr_df = pd.DataFrame({"gene_id": ["A"], "S1": [1.0], "S2": [2.0], "S5": [3.0]})
    expr_path = tmp_path / "expression.csv"
    expr_df.to_csv(expr_path, index=False)

    meta_df = pd.DataFrame(
        {"sample_id": ["S1", "S2", "S6"], "group": ["normal", "normal", "tumor"]}
    )
    meta_path = tmp_path / "metadata.csv"
    meta_df.to_csv(meta_path, index=False)

    expression = validate_expression_matrix(expr_path, gene_id_column="gene_id")
    metadata = validate_metadata(meta_path, sample_id_column="sample_id", group_column="group")

    result = match_samples(expression, metadata)
    assert result.matched_samples == ("S1", "S2")
    assert result.expression_only_samples == ("S5",)
    assert result.metadata_only_samples == ("S6",)


def test_zero_sample_overlap_raises(tmp_path):
    expr_df = pd.DataFrame({"gene_id": ["A"], "S1": [1.0], "S2": [2.0]})
    expr_path = tmp_path / "expression.csv"
    expr_df.to_csv(expr_path, index=False)

    meta_df = pd.DataFrame({"sample_id": ["S9", "S10"], "group": ["normal", "tumor"]})
    meta_path = tmp_path / "metadata.csv"
    meta_df.to_csv(meta_path, index=False)

    expression = validate_expression_matrix(expr_path, gene_id_column="gene_id")
    metadata = validate_metadata(meta_path, sample_id_column="sample_id", group_column="group")

    with pytest.raises(SampleMismatchError):
        match_samples(expression, metadata)