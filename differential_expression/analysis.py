"""
Purpose:
    Orchestrate the complete Phase 3 differential expression pipeline:
    validation -> comparison -> filtering -> statistical method ->
    multiple-testing correction -> QC -> reproducibility metadata ->
    result writing, per specification Sections 3-13.

Responsibilities:
    - Execute every stage in the deterministic order defined by the
      specification, delegating all actual logic to the module that
      owns it (this module contains no validation, comparison,
      filtering, statistical, correction, QC, reproducibility, or
      result-writing logic of its own).
    - Log every major lifecycle event via the existing project logging
      infrastructure (src.logging_manager).
    - Catch every DE-specific exception at its stage boundary, convert
      it into a structured FailureInfo, and return a FAILED
      AnalysisResult -- never letting a failure look like a completed
      analysis, never swallowing the underlying exception silently
      (the original exception is always logged and chained via
      FailureInfo.cause).
    - Never retry with a different method, comparison, threshold, or
      sample/gene set after a failure. A failure ends the analysis.

Scope:
    This module never downloads data, never mutates the input
    expression matrix or metadata files, and never performs any
    analysis beyond what is defined by the specification (no
    biomarker ranking, no pathway analysis, no visualization).
"""

from src.differential_expression.comparison import resolve_comparison
from src.differential_expression.exceptions import DEError
from src.differential_expression.filtering import apply_gene_filter
from src.differential_expression.methods import execute_statistical_method
from src.differential_expression.models import (
    AnalysisResult,
    AnalysisStatus,
    DEAnalysisConfiguration,
    FailureInfo,
)
from src.differential_expression.multiple_testing import (
    apply_multiple_testing_correction,
)
from src.differential_expression.qc import build_qc_report, verify_post_analysis_quality
from src.differential_expression.reproducibility import build_analysis_metadata
from src.differential_expression.results import write_analysis_outputs
from src.differential_expression.validation import (
    match_samples,
    validate_expression_matrix,
    validate_metadata,
)
from src.logging_manager import get_logger

ANALYSIS_VERSION = "1.0"

__all__ = [
    "ANALYSIS_VERSION",
    "run_analysis",
]

_logger = get_logger(__name__)

# Stage names used consistently in log messages and FailureInfo, in
# the deterministic order the pipeline actually executes them.
_STAGE_VALIDATE_EXPRESSION = "validate_expression_matrix"
_STAGE_VALIDATE_METADATA = "validate_metadata"
_STAGE_MATCH_SAMPLES = "match_samples"
_STAGE_RESOLVE_COMPARISON = "resolve_comparison"
_STAGE_APPLY_FILTER = "apply_gene_filter"
_STAGE_EXECUTE_METHOD = "execute_statistical_method"
_STAGE_MULTIPLE_TESTING = "apply_multiple_testing_correction"
_STAGE_POST_ANALYSIS_QC = "verify_post_analysis_quality"
_STAGE_BUILD_QC = "build_qc_report"
_STAGE_REPRODUCIBILITY = "build_analysis_metadata"
_STAGE_WRITE_OUTPUTS = "write_analysis_outputs"


def _failure_result(
    analysis_id: str, stage: str, error: Exception
) -> AnalysisResult:
    """Build a FAILED AnalysisResult from a caught exception, logging
    the failure before returning.

    Args:
        analysis_id: The analysis identifier.
        stage: The pipeline stage at which the failure occurred.
        error: The exception that was raised.

    Returns:
        A FAILED AnalysisResult carrying structured FailureInfo.
    """
    category = type(error).__name__
    message = str(error)

    _logger.error(
        "Analysis '%s' failed at stage '%s' (%s): %s",
        analysis_id,
        stage,
        category,
        message,
    )

    return AnalysisResult(
        analysis_id=analysis_id,
        status=AnalysisStatus.FAILED,
        failure=FailureInfo(
            stage=stage,
            category=category,
            message=message,
            cause=repr(error.__cause__) if error.__cause__ is not None else None,
        ),
    )


def run_analysis(configuration: DEAnalysisConfiguration) -> AnalysisResult:
    """Run the complete Phase 3 differential expression pipeline for
    one analysis.

    Args:
        configuration: A validated DEAnalysisConfiguration (typically
            produced by configuration.build_configuration()).

    Returns:
        An AnalysisResult. status == AnalysisStatus.SUCCEEDED only
        when gene_results, qc_report, metadata, and output_paths are
        all populated; status == AnalysisStatus.FAILED with
        structured FailureInfo for any stage failure. This function
        never raises -- every DE-specific exception is caught at its
        stage boundary and converted into a FAILED AnalysisResult, so
        callers always receive a structured result rather than having
        to catch exceptions themselves.
    """
    analysis_id = configuration.analysis_id

    _logger.info(
        "Starting differential expression analysis '%s' (cohort=%s, "
        "method=%s, representation=%s).",
        analysis_id,
        configuration.cancer_cohort.value,
        configuration.statistical_method.value,
        configuration.expression_representation.value,
    )

    # --- Stage 1: input validation ---
    try:
        _logger.info(
            "Analysis '%s': validating expression matrix at %s.",
            analysis_id,
            configuration.expression_matrix_path,
        )
        expression = validate_expression_matrix(
            configuration.expression_matrix_path,
            configuration.gene_id_column,
        )
    except DEError as error:
        return _failure_result(analysis_id, _STAGE_VALIDATE_EXPRESSION, error)

    try:
        _logger.info(
            "Analysis '%s': validating metadata at %s.",
            analysis_id,
            configuration.metadata_path,
        )
        metadata = validate_metadata(
            configuration.metadata_path,
            configuration.sample_id_column,
            configuration.group_column,
        )
    except DEError as error:
        return _failure_result(analysis_id, _STAGE_VALIDATE_METADATA, error)

    # --- Stage 2: cross-check samples ---
    try:
        _logger.info("Analysis '%s': matching expression and metadata samples.", analysis_id)
        sample_match = match_samples(expression, metadata)
    except DEError as error:
        return _failure_result(analysis_id, _STAGE_MATCH_SAMPLES, error)

    # --- Stage 3: resolve the explicit two-group comparison ---
    try:
        _logger.info(
            "Analysis '%s': resolving comparison '%s' vs '%s' on column '%s'.",
            analysis_id,
            configuration.comparison_group,
            configuration.reference_group,
            configuration.group_column,
        )
        group_assignment = resolve_comparison(configuration, metadata, sample_match)
    except DEError as error:
        return _failure_result(analysis_id, _STAGE_RESOLVE_COMPARISON, error)

    _logger.info(
        "Analysis '%s': comparison resolved (reference=%d sample(s), "
        "comparison=%d sample(s), excluded=%d sample(s)).",
        analysis_id,
        len(group_assignment.reference_sample_ids),
        len(group_assignment.comparison_sample_ids),
        len(group_assignment.excluded_samples),
    )

    # --- Stage 4: apply explicit gene filtering ---
    try:
        _logger.info("Analysis '%s': applying gene filter.", analysis_id)
        gene_filter_result = apply_gene_filter(
            configuration.gene_filter, expression, group_assignment
        )
    except DEError as error:
        return _failure_result(analysis_id, _STAGE_APPLY_FILTER, error)

    _logger.info(
        "Analysis '%s': gene filtering complete (%d tested, %d filtered "
        "of %d input gene(s)).",
        analysis_id,
        len(gene_filter_result.tested_gene_ids),
        len(gene_filter_result.filtered_gene_ids),
        gene_filter_result.input_gene_count,
    )

    # --- Stage 5: execute the configured statistical method ---
    try:
        _logger.info(
            "Analysis '%s': executing statistical method '%s'.",
            analysis_id,
            configuration.statistical_method.value,
        )
        gene_results = execute_statistical_method(
            configuration, expression, group_assignment, gene_filter_result
        )
    except DEError as error:
        return _failure_result(analysis_id, _STAGE_EXECUTE_METHOD, error)

    _logger.info(
        "Analysis '%s': statistical method execution complete (%d "
        "gene(s) processed).",
        analysis_id,
        len(gene_results),
    )

    # --- Stage 6: multiple-testing correction ---
    try:
        _logger.info(
            "Analysis '%s': applying multiple-testing correction (%s).",
            analysis_id,
            configuration.multiple_testing_method.value,
        )
        gene_results = apply_multiple_testing_correction(
            gene_results,
            configuration.multiple_testing_method,
            configuration.significance_threshold,
        )
    except DEError as error:
        return _failure_result(analysis_id, _STAGE_MULTIPLE_TESTING, error)

    # --- Stage 7: post-analysis QC ---
    try:
        _logger.info("Analysis '%s': running post-analysis QC checks.", analysis_id)
        verify_post_analysis_quality(gene_results, gene_filter_result, configuration)
    except DEError as error:
        return _failure_result(analysis_id, _STAGE_POST_ANALYSIS_QC, error)

    try:
        qc_report = build_qc_report(
            configuration, sample_match, group_assignment, gene_filter_result
        )
    except DEError as error:
        return _failure_result(analysis_id, _STAGE_BUILD_QC, error)

    # --- Stage 8: reproducibility metadata ---
    try:
        _logger.info("Analysis '%s': constructing reproducibility metadata.", analysis_id)
        metadata_record = build_analysis_metadata(
            configuration,
            group_assignment,
            gene_filter_result.criterion_description,
            AnalysisStatus.SUCCEEDED,
        )
    except DEError as error:
        return _failure_result(analysis_id, _STAGE_REPRODUCIBILITY, error)

    # --- Stage 9: write outputs ---
    try:
        _logger.info(
            "Analysis '%s': writing outputs to %s.",
            analysis_id,
            configuration.output_dir,
        )
        output_paths = write_analysis_outputs(
            analysis_id,
            configuration.output_dir,
            gene_results,
            metadata_record,
            qc_report,
        )
    except DEError as error:
        return _failure_result(analysis_id, _STAGE_WRITE_OUTPUTS, error)

    _logger.info(
        "Analysis '%s' completed successfully. Outputs: %s.",
        analysis_id,
        output_paths,
    )

    return AnalysisResult(
        analysis_id=analysis_id,
        status=AnalysisStatus.SUCCEEDED,
        gene_results=gene_results,
        qc_report=qc_report,
        metadata=metadata_record,
        output_paths=output_paths,
    )