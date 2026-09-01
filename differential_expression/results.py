"""
Purpose:
    Assemble the final Phase 3 outputs (differential expression result
    table, analysis metadata, QC report) and write them safely to
    disk, per specification Section 12.

Responsibilities:
    - Serialize one row per tested gene into a CSV result table
      (specification Section 12.1/12.5), preserving genes with
      missing statistical values (Section 6.4) rather than dropping
      them.
    - Serialize AnalysisMetadata and QCReport to JSON (Section 12.5).
    - Write every output under a single per-analysis staging
      directory first, and only move the completed files into their
      final location once every write has succeeded -- so a failure
      partway through writing never leaves a partial, misleading
      output visible at the final path (specification Section 10.4 /
      12.4: partially written output must never be presented as a
      completed analysis, and every output must be associated with
      the same analysis identifier).

Scope:
    This module performs no statistical computation and makes no
    scientific decisions. It uses the existing project I/O
    infrastructure (src.io_utils, src.file_utils) rather than
    re-implementing file writing.
"""

import csv
import io
import shutil
from dataclasses import asdict
from pathlib import Path

from src.differential_expression.exceptions import ResultWritingError
from src.differential_expression.models import (
    AnalysisMetadata,
    GeneResult,
    QCReport,
)
from src.file_utils import ensure_directory_exists, is_writable
from src.io_utils import move_file, write_json_file, write_text_file

RESULTS_VERSION = "1.0"

RESULTS_FILENAME_TEMPLATE = "{analysis_id}_differential_expression_results.csv"
METADATA_FILENAME_TEMPLATE = "{analysis_id}_analysis_metadata.json"
QC_FILENAME_TEMPLATE = "{analysis_id}_qc_report.json"

_GENE_RESULT_COLUMNS = (
    "gene_id",
    "tested",
    "effect_size",
    "effect_size_label",
    "raw_p_value",
    "adjusted_p_value",
    "significant",
    "missing_reason",
)

__all__ = [
    "RESULTS_VERSION",
    "RESULTS_FILENAME_TEMPLATE",
    "METADATA_FILENAME_TEMPLATE",
    "QC_FILENAME_TEMPLATE",
    "write_analysis_outputs",
]


def _gene_results_to_csv(gene_results: tuple[GeneResult, ...]) -> str:
    """Serialize gene-level results to CSV text, one row per tested
    gene, preserving genes with missing statistical values.

    Args:
        gene_results: The complete, corrected gene-level results.

    Returns:
        CSV text with a header row and one row per gene in
        `gene_results`, in the same order.
    """
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(_GENE_RESULT_COLUMNS)

    for result in gene_results:
        writer.writerow(
            [
                result.gene_id,
                result.tested,
                result.effect_size if result.effect_size is not None else "",
                result.effect_size_label or "",
                result.raw_p_value if result.raw_p_value is not None else "",
                (
                    result.adjusted_p_value
                    if result.adjusted_p_value is not None
                    else ""
                ),
                (
                    ""
                    if result.significant is None
                    else str(result.significant)
                ),
                result.missing_reason or "",
            ]
        )

    return buffer.getvalue()


def _metadata_to_json_dict(metadata: AnalysisMetadata) -> dict:
    """Convert AnalysisMetadata to a JSON-serializable dict.

    Path, enum, and datetime fields are not natively JSON-serializable
    and are converted explicitly rather than relying on a generic
    default= callback that could silently stringify an unexpected
    type.

    Args:
        metadata: The analysis metadata to serialize.

    Returns:
        A dict containing only JSON-serializable values.
    """
    raw = asdict(metadata)

    raw["expression_matrix_path"] = str(metadata.expression_matrix_path)
    raw["metadata_path"] = str(metadata.metadata_path)
    raw["timestamp_utc"] = metadata.timestamp_utc.isoformat()
    raw["cancer_cohort"] = metadata.cancer_cohort.value
    raw["expression_representation"] = metadata.expression_representation.value
    raw["statistical_method"] = metadata.statistical_method.value
    raw["multiple_testing_method"] = metadata.multiple_testing_method.value
    raw["analysis_status"] = metadata.analysis_status.value
    raw["included_sample_ids"] = list(metadata.included_sample_ids)

    return raw


def _qc_report_to_json_dict(qc_report: QCReport) -> dict:
    """Convert QCReport to a JSON-serializable dict.

    Args:
        qc_report: The QC report to serialize.

    Returns:
        A dict containing only JSON-serializable values.
    """
    raw = asdict(qc_report)

    raw["excluded_samples"] = [
        {
            "sample_id": excluded.sample_id,
            "reason": excluded.reason,
            "stage": excluded.stage,
        }
        for excluded in qc_report.excluded_samples
    ]
    raw["statistical_method"] = qc_report.statistical_method.value
    raw["multiple_testing_method"] = qc_report.multiple_testing_method.value
    raw["notes"] = list(qc_report.notes)

    return raw


def write_analysis_outputs(
    analysis_id: str,
    output_dir: Path,
    gene_results: tuple[GeneResult, ...],
    metadata: AnalysisMetadata,
    qc_report: QCReport,
) -> dict[str, Path]:
    """Write the complete set of Phase 3 outputs for one analysis.

    All three outputs are written to a per-analysis staging directory
    first; only if every write succeeds are the files moved into
    `output_dir` under their final names. If any step fails, the
    staging directory (and therefore any partial output) is removed,
    and ResultWritingError is raised -- `output_dir` never contains a
    partial result set for this analysis_id.

    Args:
        analysis_id: The analysis identifier, used to build output
            filenames and the staging directory name.
        output_dir: The directory the final output files belong in.
        gene_results: The complete, corrected gene-level results.
        metadata: The analysis's reproducibility metadata.
        qc_report: The analysis's QC report.

    Returns:
        A dict mapping "results", "metadata", and "qc_report" to the
        final Path of each written file.

    Raises:
        ResultWritingError: If `output_dir` is not writable, or if
            any output cannot be safely written.
    """
    if not isinstance(output_dir, Path):
        raise ResultWritingError(
            f"'output_dir' must be a pathlib.Path, got {type(output_dir).__name__}."
        )

    try:
        ensure_directory_exists(output_dir)
    except (NotADirectoryError, PermissionError, OSError) as error:
        raise ResultWritingError(
            f"Could not create or access output directory {output_dir}: {error}"
        ) from error

    if not is_writable(output_dir):
        raise ResultWritingError(
            f"Output directory is not writable: {output_dir}"
        )

    staging_dir = output_dir / f".{analysis_id}.staging"

    try:
        ensure_directory_exists(staging_dir)
    except (NotADirectoryError, PermissionError, OSError) as error:
        raise ResultWritingError(
            f"Could not create staging directory {staging_dir}: {error}"
        ) from error

    results_filename = RESULTS_FILENAME_TEMPLATE.format(analysis_id=analysis_id)
    metadata_filename = METADATA_FILENAME_TEMPLATE.format(analysis_id=analysis_id)
    qc_filename = QC_FILENAME_TEMPLATE.format(analysis_id=analysis_id)

    staged_results_path = staging_dir / results_filename
    staged_metadata_path = staging_dir / metadata_filename
    staged_qc_path = staging_dir / qc_filename

    try:
        write_text_file(staged_results_path, _gene_results_to_csv(gene_results))
        write_json_file(staged_metadata_path, _metadata_to_json_dict(metadata))
        write_json_file(staged_qc_path, _qc_report_to_json_dict(qc_report))
    except (OSError, TypeError) as error:
        shutil.rmtree(staging_dir, ignore_errors=True)
        raise ResultWritingError(
            f"Failed to write Phase 3 outputs for analysis "
            f"'{analysis_id}': {error}. No partial output has been "
            "left in the final output directory."
        ) from error

    final_results_path = output_dir / results_filename
    final_metadata_path = output_dir / metadata_filename
    final_qc_path = output_dir / qc_filename

    try:
        move_file(staged_results_path, final_results_path)
        move_file(staged_metadata_path, final_metadata_path)
        move_file(staged_qc_path, final_qc_path)
    except (OSError, FileNotFoundError, FileExistsError, IsADirectoryError) as error:
        shutil.rmtree(staging_dir, ignore_errors=True)
        raise ResultWritingError(
            f"Failed to finalize Phase 3 outputs for analysis "
            f"'{analysis_id}': {error}."
        ) from error

    shutil.rmtree(staging_dir, ignore_errors=True)

    return {
        "results": final_results_path,
        "metadata": final_metadata_path,
        "qc_report": final_qc_path,
    }