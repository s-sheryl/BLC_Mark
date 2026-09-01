"""Tests for src.differential_expression.results."""

import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import pytest

from src.dataset_registry import CancerType
from src.differential_expression.exceptions import ResultWritingError
from src.differential_expression.models import (
    AnalysisMetadata,
    AnalysisStatus,
    ExpressionRepresentation,
    GeneResult,
    MultipleTestingMethod,
    QCReport,
    StatisticalMethod,
)
from src.differential_expression.results import write_analysis_outputs


def _gene_results():
    return (
        GeneResult(
            gene_id="A", tested=True, effect_size=1.5, effect_size_label="log2_fold_change",
            raw_p_value=0.01, adjusted_p_value=0.02, significant=True, missing_reason=None,
        ),
        GeneResult(
            gene_id="B", tested=True, effect_size=None, effect_size_label=None,
            raw_p_value=None, adjusted_p_value=None, significant=None,
            missing_reason="Insufficient non-missing values.",
        ),
    )


def _metadata():
    return AnalysisMetadata(
        analysis_id="a1",
        timestamp_utc=datetime.now(timezone.utc),
        cancer_cohort=CancerType.BRCA,
        expression_matrix_path=Path("/tmp/expression.csv"),
        metadata_path=Path("/tmp/metadata.csv"),
        expression_matrix_sha256="a" * 64,
        metadata_sha256="b" * 64,
        expression_representation=ExpressionRepresentation.NORMALIZED_LOG2,
        statistical_method=StatisticalMethod.WELCH_T_TEST,
        statistical_method_version="scipy 1.17.1",
        design="Two-group comparison.",
        reference_group="normal",
        comparison_group="tumor",
        included_sample_ids=("S1", "S2", "S3", "S4"),
        reference_group_size=2,
        comparison_group_size=2,
        gene_filter_criterion="No filtering.",
        multiple_testing_method=MultipleTestingMethod.BENJAMINI_HOCHBERG,
        significance_threshold=0.05,
        effect_size_threshold=None,
        python_version="3.12.3",
        package_versions={"scipy": "1.17.1"},
        de_package_version="1.0",
        blc_mark_version=None,
        analysis_status=AnalysisStatus.SUCCEEDED,
    )


def _qc_report():
    return QCReport(
        input_gene_count=2,
        tested_gene_count=2,
        filtered_gene_count=0,
        initial_sample_count=4,
        included_sample_count=4,
        excluded_sample_count=0,
        excluded_samples=(),
        reference_group_size=2,
        comparison_group_size=2,
        gene_filter_criterion="No filtering.",
        statistical_method=StatisticalMethod.WELCH_T_TEST,
        multiple_testing_method=MultipleTestingMethod.BENJAMINI_HOCHBERG,
        significance_threshold=0.05,
    )


def test_write_analysis_outputs_creates_all_three_files(tmp_path):
    output_dir = tmp_path / "output"
    paths = write_analysis_outputs(
        "a1", output_dir, _gene_results(), _metadata(), _qc_report()
    )
    assert paths["results"].exists()
    assert paths["metadata"].exists()
    assert paths["qc_report"].exists()


def test_results_csv_preserves_missing_statistical_values(tmp_path):
    output_dir = tmp_path / "output"
    paths = write_analysis_outputs(
        "a1", output_dir, _gene_results(), _metadata(), _qc_report()
    )
    with open(paths["results"], newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 2  # gene B retained despite missing values
    gene_a = next(r for r in rows if r["gene_id"] == "A")
    gene_b = next(r for r in rows if r["gene_id"] == "B")
    assert gene_a["raw_p_value"] == "0.01"
    assert gene_b["raw_p_value"] == ""
    assert gene_b["missing_reason"] == "Insufficient non-missing values."


def test_metadata_json_is_valid_and_complete(tmp_path):
    output_dir = tmp_path / "output"
    paths = write_analysis_outputs(
        "a1", output_dir, _gene_results(), _metadata(), _qc_report()
    )
    data = json.loads(paths["metadata"].read_text())
    assert data["analysis_id"] == "a1"
    assert data["statistical_method"] == "welch_t_test"
    assert data["cancer_cohort"] == "TCGA-BRCA"
    assert set(data["included_sample_ids"]) == {"S1", "S2", "S3", "S4"}


def test_qc_json_is_valid_and_complete(tmp_path):
    output_dir = tmp_path / "output"
    paths = write_analysis_outputs(
        "a1", output_dir, _gene_results(), _metadata(), _qc_report()
    )
    data = json.loads(paths["qc_report"].read_text())
    assert data["tested_gene_count"] == 2
    assert data["statistical_method"] == "welch_t_test"


def test_outputs_named_with_analysis_id(tmp_path):
    output_dir = tmp_path / "output"
    paths = write_analysis_outputs(
        "my-analysis-42", output_dir, _gene_results(), _metadata(), _qc_report()
    )
    assert "my-analysis-42" in paths["results"].name
    assert "my-analysis-42" in paths["metadata"].name
    assert "my-analysis-42" in paths["qc_report"].name


def test_staging_directory_removed_after_success(tmp_path):
    output_dir = tmp_path / "output"
    write_analysis_outputs("a1", output_dir, _gene_results(), _metadata(), _qc_report())
    assert not (output_dir / ".a1.staging").exists()


def test_unwritable_output_location_raises_result_writing_error(tmp_path):
    # Point output_dir at a path whose parent does not exist and
    # cannot be created (a file exists where a directory is needed).
    blocking_file = tmp_path / "blocking"
    blocking_file.write_text("not a directory")
    output_dir = blocking_file / "output"

    with pytest.raises(ResultWritingError):
        write_analysis_outputs("a1", output_dir, _gene_results(), _metadata(), _qc_report())


def test_failure_does_not_leave_partial_output_in_final_directory(tmp_path, monkeypatch):
    output_dir = tmp_path / "output"

    # Force a failure partway through writing (after the CSV write
    # succeeds, JSON metadata write fails) by making metadata
    # non-serializable via monkeypatching write_json_file.
    import src.differential_expression.results as results_module

    original_write_json = results_module.write_json_file

    call_count = {"n": 0}

    def _flaky_write_json(path, data, **kwargs):
        call_count["n"] += 1
        if call_count["n"] == 1:
            raise OSError("Simulated disk failure.")
        return original_write_json(path, data, **kwargs)

    monkeypatch.setattr(results_module, "write_json_file", _flaky_write_json)

    with pytest.raises(ResultWritingError):
        write_analysis_outputs("a1", output_dir, _gene_results(), _metadata(), _qc_report())

    # No file for this analysis should exist in the final output
    # directory -- a partial write must never look like success.
    if output_dir.exists():
        remaining = list(output_dir.iterdir())
        assert not any("a1" in path.name for path in remaining)