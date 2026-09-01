"""Tests for src.differential_expression.reproducibility."""

from pathlib import Path

import pytest

from src.differential_expression.configuration import build_configuration
from src.differential_expression.exceptions import ReproducibilityError
from src.differential_expression.models import AnalysisStatus, GroupAssignment
from src.differential_expression.reproducibility import build_analysis_metadata
from src.hash_utils import hash_file


def _write_temp_file(tmp_path: Path, name: str, content: str) -> Path:
    path = tmp_path / name
    path.write_text(content)
    return path


def _config(tmp_path, **overrides):
    expr_path = _write_temp_file(tmp_path, "expression.csv", "gene_id,S1,S2,S3,S4\nA,1,2,3,4\n")
    meta_path = _write_temp_file(
        tmp_path, "metadata.csv", "sample_id,group\nS1,normal\nS2,normal\nS3,tumor\nS4,tumor\n"
    )
    kwargs = dict(
        analysis_id="a1",
        cancer_cohort="TCGA-BRCA",
        expression_matrix_path=expr_path,
        metadata_path=meta_path,
        gene_id_column="gene_id",
        sample_id_column="sample_id",
        group_column="group",
        reference_group="normal",
        comparison_group="tumor",
        expression_representation="normalized_log2",
        statistical_method="welch_t_test",
        minimum_replicates_per_group=2,
        output_dir=tmp_path / "output",
    )
    kwargs.update(overrides)
    return build_configuration(**kwargs)


def _group_assignment():
    return GroupAssignment(
        reference_group="normal",
        comparison_group="tumor",
        reference_sample_ids=("S1", "S2"),
        comparison_sample_ids=("S3", "S4"),
    )


def test_input_hashes_match_hash_utils(tmp_path):
    config = _config(tmp_path)
    metadata_record = build_analysis_metadata(
        config, _group_assignment(), "No filtering.", AnalysisStatus.SUCCEEDED
    )
    assert metadata_record.expression_matrix_sha256 == hash_file(config.expression_matrix_path)
    assert metadata_record.metadata_sha256 == hash_file(config.metadata_path)


def test_configuration_captured(tmp_path):
    config = _config(tmp_path)
    metadata_record = build_analysis_metadata(
        config, _group_assignment(), "No filtering.", AnalysisStatus.SUCCEEDED
    )
    assert metadata_record.analysis_id == "a1"
    assert metadata_record.reference_group == "normal"
    assert metadata_record.comparison_group == "tumor"
    assert metadata_record.statistical_method == config.statistical_method
    assert metadata_record.multiple_testing_method == config.multiple_testing_method
    assert metadata_record.significance_threshold == config.significance_threshold
    assert metadata_record.gene_filter_criterion == "No filtering."


def test_included_sample_ids_captured(tmp_path):
    config = _config(tmp_path)
    metadata_record = build_analysis_metadata(
        config, _group_assignment(), "No filtering.", AnalysisStatus.SUCCEEDED
    )
    assert set(metadata_record.included_sample_ids) == {"S1", "S2", "S3", "S4"}
    assert metadata_record.reference_group_size == 2
    assert metadata_record.comparison_group_size == 2


def test_software_version_recorded_not_fabricated(tmp_path):
    config = _config(tmp_path)
    metadata_record = build_analysis_metadata(
        config, _group_assignment(), "No filtering.", AnalysisStatus.SUCCEEDED
    )
    assert metadata_record.statistical_method_version is not None
    assert "scipy" in metadata_record.statistical_method_version
    assert "scipy" in metadata_record.package_versions
    assert metadata_record.python_version  # non-empty


def test_deseq2_method_version_is_none_not_fabricated(tmp_path):
    config = _config(tmp_path, statistical_method="deseq2", expression_representation="raw_counts")
    metadata_record = build_analysis_metadata(
        config, _group_assignment(), "No filtering.", AnalysisStatus.SUCCEEDED
    )
    assert metadata_record.statistical_method_version is None
    assert metadata_record.package_versions == {}


def test_missing_expression_file_raises_reproducibility_error(tmp_path):
    config = _config(tmp_path)
    config.expression_matrix_path.unlink()
    with pytest.raises(ReproducibilityError):
        build_analysis_metadata(
            config, _group_assignment(), "No filtering.", AnalysisStatus.SUCCEEDED
        )


def test_deterministic_metadata_construction_hashes(tmp_path):
    config = _config(tmp_path)
    first = build_analysis_metadata(
        config, _group_assignment(), "No filtering.", AnalysisStatus.SUCCEEDED
    )
    second = build_analysis_metadata(
        config, _group_assignment(), "No filtering.", AnalysisStatus.SUCCEEDED
    )
    # Hashes and configuration-derived fields must be identical across
    # repeated construction; only the timestamp legitimately varies.
    assert first.expression_matrix_sha256 == second.expression_matrix_sha256
    assert first.metadata_sha256 == second.metadata_sha256
    assert first.included_sample_ids == second.included_sample_ids
    assert first.statistical_method_version == second.statistical_method_version