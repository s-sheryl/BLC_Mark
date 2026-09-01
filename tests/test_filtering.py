"""Tests for src.differential_expression.filtering."""

from pathlib import Path

import pandas as pd
import pytest

from src.differential_expression.filtering import apply_gene_filter
from src.differential_expression.models import (
    GeneFilterConfiguration,
    GroupAssignment,
    ValidatedExpressionMatrix,
)


def _expression_matrix():
    df = pd.DataFrame(
        {
            "gene_id": ["LOW", "MID", "HIGH"],
            "S1": [0.1, 5.0, 100.0],
            "S2": [0.2, 5.5, 110.0],
            "S3": [0.1, 5.2, 105.0],
            "S4": [0.3, 4.8, 95.0],
        }
    )
    return ValidatedExpressionMatrix(
        file_path=Path("/tmp/expression.csv"),
        gene_id_column="gene_id",
        sample_columns=("S1", "S2", "S3", "S4"),
        gene_ids=("LOW", "MID", "HIGH"),
        dataframe=df,
    )


def _group_assignment():
    return GroupAssignment(
        reference_group="normal",
        comparison_group="tumor",
        reference_sample_ids=("S1", "S2"),
        comparison_sample_ids=("S3", "S4"),
    )


def test_no_filter_keeps_all_genes():
    result = apply_gene_filter(
        GeneFilterConfiguration(), _expression_matrix(), _group_assignment()
    )
    assert result.tested_gene_ids == ("LOW", "MID", "HIGH")
    assert result.filtered_gene_ids == ()
    assert result.input_gene_count == 3


def test_filter_removes_low_expression_genes():
    gene_filter = GeneFilterConfiguration(apply_filter=True, minimum_mean_expression=1.0)
    result = apply_gene_filter(gene_filter, _expression_matrix(), _group_assignment())
    assert result.tested_gene_ids == ("MID", "HIGH")
    assert result.filtered_gene_ids == ("LOW",)
    assert result.input_gene_count == 3


def test_filter_is_deterministic():
    gene_filter = GeneFilterConfiguration(apply_filter=True, minimum_mean_expression=1.0)
    first = apply_gene_filter(gene_filter, _expression_matrix(), _group_assignment())
    second = apply_gene_filter(gene_filter, _expression_matrix(), _group_assignment())
    assert first.tested_gene_ids == second.tested_gene_ids
    assert first.filtered_gene_ids == second.filtered_gene_ids


def test_filter_only_uses_included_samples():
    # A sample not included in the comparison (e.g. excluded during
    # comparison resolution) must not influence the mean used for
    # filtering. Here S4 is deliberately excluded from the group
    # assignment and has a very different value from the others.
    df = pd.DataFrame(
        {
            "gene_id": ["GENE_A"],
            "S1": [2.0],
            "S2": [2.0],
            "S3": [2.0],
            "S4": [1000.0],
        }
    )
    expression = ValidatedExpressionMatrix(
        file_path=Path("/tmp/expression.csv"),
        gene_id_column="gene_id",
        sample_columns=("S1", "S2", "S3", "S4"),
        gene_ids=("GENE_A",),
        dataframe=df,
    )
    group_assignment = GroupAssignment(
        reference_group="normal",
        comparison_group="tumor",
        reference_sample_ids=("S1",),
        comparison_sample_ids=("S2", "S3"),
    )
    gene_filter = GeneFilterConfiguration(apply_filter=True, minimum_mean_expression=1.5)
    result = apply_gene_filter(gene_filter, expression, group_assignment)
    assert result.tested_gene_ids == ("GENE_A",)


def test_filtered_and_tested_are_disjoint_and_complete():
    gene_filter = GeneFilterConfiguration(apply_filter=True, minimum_mean_expression=50.0)
    result = apply_gene_filter(gene_filter, _expression_matrix(), _group_assignment())
    assert set(result.tested_gene_ids) & set(result.filtered_gene_ids) == set()
    assert len(result.tested_gene_ids) + len(result.filtered_gene_ids) == result.input_gene_count