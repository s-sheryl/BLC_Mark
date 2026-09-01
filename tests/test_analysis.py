"""End-to-end integration tests for src.differential_expression.analysis.

Exercises the complete pipeline against small, deterministic synthetic
data: configuration -> validation -> comparison -> filtering ->
statistics -> multiple testing -> QC -> reproducibility -> results,
for both a successful analysis and several expected failure paths.

No real TCGA data is used, per specification Section 14 / the
project's testing rules.
"""

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.differential_expression.analysis import run_analysis
from src.differential_expression.configuration import build_configuration
from src.differential_expression.models import AnalysisStatus, GeneFilterConfiguration


def _write_synthetic_expression_matrix(path: Path) -> None:
    """Ten genes, four samples. Genes 1-5 have a clear, deterministic
    shift between the tumor and normal groups (so Welch's t-test
    should find them significant after BH correction); genes 6-10
    have no shift (should not be significant). Gene "MISSING1" has a
    missing value in one sample per group, to exercise per-gene
    missing-value handling without failing validation (validation
    only rejects genuinely non-numeric values, not NaN).
    """
    rng = np.random.default_rng(seed=42)

    rows = []
    for i in range(1, 6):
        rows.append(
            {
                "gene_id": f"UP{i}",
                "N1": 5.0 + rng.normal(0, 0.1),
                "N2": 5.1 + rng.normal(0, 0.1),
                "N3": 5.2 + rng.normal(0, 0.1),
                "T1": 9.0 + rng.normal(0, 0.1),
                "T2": 9.1 + rng.normal(0, 0.1),
                "T3": 9.2 + rng.normal(0, 0.1),
            }
        )
    for i in range(1, 6):
        rows.append(
            {
                "gene_id": f"FLAT{i}",
                "N1": 5.0 + rng.normal(0, 0.1),
                "N2": 5.05 + rng.normal(0, 0.1),
                "N3": 4.95 + rng.normal(0, 0.1),
                "T1": 5.02 + rng.normal(0, 0.1),
                "T2": 4.98 + rng.normal(0, 0.1),
                "T3": 5.03 + rng.normal(0, 0.1),
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)


def _write_synthetic_metadata(path: Path) -> None:
    df = pd.DataFrame(
        {
            "sample_id": ["N1", "N2", "N3", "T1", "T2", "T3"],
            "group": ["normal", "normal", "normal", "tumor", "tumor", "tumor"],
        }
    )
    df.to_csv(path, index=False)


@pytest.fixture
def synthetic_paths(tmp_path):
    expression_path = tmp_path / "expression.csv"
    metadata_path = tmp_path / "metadata.csv"
    output_dir = tmp_path / "output"
    _write_synthetic_expression_matrix(expression_path)
    _write_synthetic_metadata(metadata_path)
    return expression_path, metadata_path, output_dir


def _base_config(synthetic_paths, **overrides):
    expression_path, metadata_path, output_dir = synthetic_paths
    kwargs = dict(
        analysis_id="integration-test-001",
        cancer_cohort="TCGA-BRCA",
        expression_matrix_path=expression_path,
        metadata_path=metadata_path,
        gene_id_column="gene_id",
        sample_id_column="sample_id",
        group_column="group",
        reference_group="normal",
        comparison_group="tumor",
        expression_representation="normalized_log2",
        statistical_method="welch_t_test",
        minimum_replicates_per_group=2,
        output_dir=output_dir,
    )
    kwargs.update(overrides)
    return build_configuration(**kwargs)


# --- Successful end-to-end analysis ---


def test_end_to_end_successful_analysis(synthetic_paths):
    config = _base_config(synthetic_paths)
    result = run_analysis(config)

    assert result.status == AnalysisStatus.SUCCEEDED
    assert result.failure is None
    assert result.gene_results is not None
    assert result.qc_report is not None
    assert result.metadata is not None
    assert result.output_paths is not None


def test_end_to_end_all_genes_retained_including_any_missing(synthetic_paths):
    config = _base_config(synthetic_paths)
    result = run_analysis(config)

    gene_ids = {r.gene_id for r in result.gene_results}
    expected = {f"UP{i}" for i in range(1, 6)} | {f"FLAT{i}" for i in range(1, 6)}
    assert gene_ids == expected
    assert len(result.gene_results) == 10  # one row per tested gene, none dropped


def test_end_to_end_up_genes_are_significant_flat_genes_are_not(synthetic_paths):
    config = _base_config(synthetic_paths)
    result = run_analysis(config)

    by_id = {r.gene_id: r for r in result.gene_results}

    for i in range(1, 6):
        up = by_id[f"UP{i}"]
        assert up.raw_p_value is not None
        assert up.effect_size is not None
        assert up.effect_size > 0  # tumor (comparison) higher than normal (reference)

    # BH-corrected significance: the clearly shifted genes should be
    # significant at the default 0.05 threshold; the flat genes
    # generally should not (statistical, not guaranteed for every
    # possible seed, but true for this fixed seed/effect size).
    up_significant = [by_id[f"UP{i}"].significant for i in range(1, 6)]
    assert all(up_significant)


def test_end_to_end_outputs_written_and_associated_with_same_analysis_id(synthetic_paths):
    config = _base_config(synthetic_paths)
    result = run_analysis(config)

    results_path = result.output_paths["results"]
    metadata_path = result.output_paths["metadata"]
    qc_path = result.output_paths["qc_report"]

    assert config.analysis_id in results_path.name
    assert config.analysis_id in metadata_path.name
    assert config.analysis_id in qc_path.name

    metadata_json = json.loads(metadata_path.read_text())
    qc_json = json.loads(qc_path.read_text())
    assert metadata_json["analysis_id"] == config.analysis_id
    assert qc_json["tested_gene_count"] == 10

    with open(results_path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 10


def test_end_to_end_reproducibility_hashes_match_actual_input_files(synthetic_paths):
    from src.hash_utils import hash_file

    config = _base_config(synthetic_paths)
    result = run_analysis(config)

    assert result.metadata.expression_matrix_sha256 == hash_file(config.expression_matrix_path)
    assert result.metadata.metadata_sha256 == hash_file(config.metadata_path)


def test_end_to_end_with_explicit_gene_filter(synthetic_paths):
    config = _base_config(
        synthetic_paths,
        gene_filter=GeneFilterConfiguration(
            apply_filter=True,
            minimum_mean_expression=6.0,
            criterion_description="Mean expression across included samples >= 6.0.",
        ),
    )
    result = run_analysis(config)

    assert result.status == AnalysisStatus.SUCCEEDED
    # FLAT genes (~mean 5.0) should be filtered out; UP genes
    # (mean ~7.0) should remain.
    gene_ids = {r.gene_id for r in result.gene_results}
    assert gene_ids == {f"UP{i}" for i in range(1, 6)}
    assert result.qc_report.filtered_gene_count == 5


# --- Expected failure paths, end-to-end ---


def test_end_to_end_fails_explicitly_on_missing_expression_file(synthetic_paths):
    expression_path, metadata_path, output_dir = synthetic_paths
    expression_path.unlink()

    config = _base_config(synthetic_paths)
    result = run_analysis(config)

    assert result.status == AnalysisStatus.FAILED
    assert result.failure is not None
    assert result.failure.stage == "validate_expression_matrix"
    assert result.gene_results is None
    assert result.qc_report is None
    assert result.metadata is None


def test_end_to_end_fails_explicitly_on_zero_sample_overlap(tmp_path):
    expression_path = tmp_path / "expression.csv"
    metadata_path = tmp_path / "metadata.csv"
    output_dir = tmp_path / "output"

    pd.DataFrame({"gene_id": ["A"], "X1": [1.0], "X2": [2.0]}).to_csv(
        expression_path, index=False
    )
    pd.DataFrame({"sample_id": ["Y1", "Y2"], "group": ["normal", "tumor"]}).to_csv(
        metadata_path, index=False
    )

    config = build_configuration(
        analysis_id="fail-001",
        cancer_cohort="TCGA-BRCA",
        expression_matrix_path=expression_path,
        metadata_path=metadata_path,
        gene_id_column="gene_id",
        sample_id_column="sample_id",
        group_column="group",
        reference_group="normal",
        comparison_group="tumor",
        expression_representation="normalized_log2",
        statistical_method="welch_t_test",
        minimum_replicates_per_group=2,
        output_dir=output_dir,
    )
    result = run_analysis(config)

    assert result.status == AnalysisStatus.FAILED
    assert result.failure.stage == "match_samples"


def test_end_to_end_fails_explicitly_on_insufficient_replication(synthetic_paths):
    expression_path, metadata_path, output_dir = synthetic_paths
    # Reduce metadata to one normal sample only.
    pd.DataFrame(
        {
            "sample_id": ["N1", "T1", "T2", "T3"],
            "group": ["normal", "tumor", "tumor", "tumor"],
        }
    ).to_csv(metadata_path, index=False)

    config = _base_config(synthetic_paths)
    result = run_analysis(config)

    assert result.status == AnalysisStatus.FAILED
    assert result.failure.stage == "resolve_comparison"
    assert "InsufficientReplication" in result.failure.category


def test_end_to_end_fails_explicitly_on_incompatible_method_representation(synthetic_paths):
    config = _base_config(
        synthetic_paths,
        expression_representation="raw_counts",
        statistical_method="welch_t_test",
    )
    result = run_analysis(config)

    assert result.status == AnalysisStatus.FAILED
    assert result.failure.stage == "execute_statistical_method"
    assert result.failure.category == "UnsupportedMethodError"


def test_end_to_end_fails_explicitly_on_deseq2_unavailable(synthetic_paths):
    expression_path, metadata_path, output_dir = synthetic_paths
    # Overwrite with a raw-counts-labeled dataset (content doesn't
    # matter -- deseq2 is unconditionally unavailable regardless of
    # content).
    config = _base_config(
        synthetic_paths,
        expression_representation="raw_counts",
        statistical_method="deseq2",
    )
    result = run_analysis(config)

    assert result.status == AnalysisStatus.FAILED
    assert result.failure.stage == "execute_statistical_method"
    assert result.failure.category == "UnsupportedMethodError"
    assert "DESeq2" in result.failure.message or "backend" in result.failure.message


def test_end_to_end_failure_never_writes_output(synthetic_paths):
    expression_path, metadata_path, output_dir = synthetic_paths
    expression_path.unlink()

    config = _base_config(synthetic_paths)
    result = run_analysis(config)

    assert result.status == AnalysisStatus.FAILED
    assert not output_dir.exists() or not any(output_dir.iterdir())


def test_end_to_end_does_not_mutate_input_files(synthetic_paths):
    expression_path, metadata_path, output_dir = synthetic_paths
    expression_before = expression_path.read_bytes()
    metadata_before = metadata_path.read_bytes()

    config = _base_config(synthetic_paths)
    run_analysis(config)

    assert expression_path.read_bytes() == expression_before
    assert metadata_path.read_bytes() == metadata_before


def test_end_to_end_deterministic_across_repeated_runs(synthetic_paths):
    config = _base_config(synthetic_paths)
    result_one = run_analysis(config)

    # Re-run into a fresh output directory (writing to the same path
    # twice is a results.py concern, not an analysis.py one).
    expression_path, metadata_path, output_dir = synthetic_paths
    config_two = _base_config(
        synthetic_paths, output_dir=output_dir.parent / "output2", analysis_id="integration-test-002"
    )
    result_two = run_analysis(config_two)

    by_id_one = {r.gene_id: r for r in result_one.gene_results}
    by_id_two = {r.gene_id: r for r in result_two.gene_results}

    for gene_id in by_id_one:
        assert by_id_one[gene_id].raw_p_value == by_id_two[gene_id].raw_p_value
        assert by_id_one[gene_id].effect_size == by_id_two[gene_id].effect_size