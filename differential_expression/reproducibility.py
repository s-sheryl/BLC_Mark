"""
Purpose:
    Construct the complete reproducibility/analysis metadata record
    for one differential expression analysis, per specification
    Section 9 and Section 12.2.

Responsibilities:
    - Hash the expression matrix and metadata input files using the
      existing project hashing infrastructure (src.hash_utils), never
      re-implementing SHA-256.
    - Record the complete analysis configuration: comparison, method,
      multiple-testing method, thresholds, filtering criteria,
      included sample identifiers.
    - Record the software environment: Python version, and the
      version of any package actually used by the statistical method
      (via methods.get_method_version(), which itself never
      fabricates a version for an unavailable backend).
    - Construct a single AnalysisMetadata record associated with the
      analysis identifier.

Scope:
    This module performs no statistical computation and makes no
    QC/statistical decisions. It only records what already happened
    elsewhere in the pipeline. If an input file cannot be hashed
    (e.g. it no longer exists at the recorded path), this module
    raises ReproducibilityError rather than fabricating a hash or
    silently omitting it.
"""
from src.version import BLC_MARK_VERSION
import sys
from datetime import datetime, timezone

from src.differential_expression.exceptions import ReproducibilityError
from src.differential_expression.methods import get_method_version
from src.differential_expression.models import (
    AnalysisMetadata,
    AnalysisStatus,
    DEAnalysisConfiguration,
    GroupAssignment,
)
from src.hash_utils import hash_file

REPRODUCIBILITY_VERSION = "1.0"

DE_PACKAGE_VERSION = "1.0"


__all__ = [
    "REPRODUCIBILITY_VERSION",
    "DE_PACKAGE_VERSION",
    "BLC_MARK_VERSION",
    "build_analysis_metadata",
]


def _hash_input_file(path, label: str) -> str:
    """Hash an input file via src.hash_utils, translating any failure
    into ReproducibilityError.

    Args:
        path: Path to the file to hash.
        label: Human-readable label used in error messages.

    Returns:
        The SHA-256 hex digest of the file.

    Raises:
        ReproducibilityError: If the file cannot be hashed (e.g. it
            does not exist at the recorded path -- this should not
            happen if validation.py already ran successfully against
            the same path, but reproducibility.py does not assume
            that and fails explicitly rather than silently omitting
            the hash).
    """
    try:
        return hash_file(path)
    except (FileNotFoundError, IsADirectoryError, OSError) as error:
        raise ReproducibilityError(
            f"Could not hash {label} at {path} while constructing "
            f"reproducibility metadata: {error}"
        ) from error


def build_analysis_metadata(
    configuration: DEAnalysisConfiguration,
    group_assignment: GroupAssignment,
    gene_filter_criterion: str,
    analysis_status: AnalysisStatus,
) -> AnalysisMetadata:
    """Construct the complete reproducibility metadata for an
    analysis.

    Args:
        configuration: The analysis configuration.
        group_assignment: The resolved two-group comparison.
        gene_filter_criterion: The exact filtering criterion applied
            (from GeneFilterResult.criterion_description).
        analysis_status: The status to record -- typically SUCCEEDED
            when this is called after a completed analysis.

    Returns:
        A complete AnalysisMetadata record.

    Raises:
        ReproducibilityError: If either input file cannot be hashed.
    """
    expression_sha256 = _hash_input_file(
        configuration.expression_matrix_path, "the expression matrix"
    )
    metadata_sha256 = _hash_input_file(
        configuration.metadata_path, "the sample metadata"
    )

    included_sample_ids = tuple(
        sorted(
            set(group_assignment.reference_sample_ids)
            | set(group_assignment.comparison_sample_ids)
        )
    )

    method_version = get_method_version(configuration.statistical_method)

    package_versions: dict[str, str] = {}
    if method_version is not None:
        # get_method_version() returns "<package> <version>"; split
        # once so package_versions is keyed by package name.
        package_name, _, version_value = method_version.partition(" ")
        if package_name and version_value:
            package_versions[package_name] = version_value

    design = (
        f"Two-group comparison: '{configuration.comparison_group}' "
        f"(comparison) vs '{configuration.reference_group}' "
        "(reference), on column "
        f"'{configuration.group_column}'."
    )

    return AnalysisMetadata(
        analysis_id=configuration.analysis_id,
        timestamp_utc=datetime.now(timezone.utc),
        cancer_cohort=configuration.cancer_cohort,
        expression_matrix_path=configuration.expression_matrix_path,
        metadata_path=configuration.metadata_path,
        expression_matrix_sha256=expression_sha256,
        metadata_sha256=metadata_sha256,
        expression_representation=configuration.expression_representation,
        statistical_method=configuration.statistical_method,
        statistical_method_version=method_version,
        design=design,
        reference_group=configuration.reference_group,
        comparison_group=configuration.comparison_group,
        included_sample_ids=included_sample_ids,
        reference_group_size=len(group_assignment.reference_sample_ids),
        comparison_group_size=len(group_assignment.comparison_sample_ids),
        gene_filter_criterion=gene_filter_criterion,
        multiple_testing_method=configuration.multiple_testing_method,
        significance_threshold=configuration.significance_threshold,
        effect_size_threshold=configuration.effect_size_threshold,
        python_version=sys.version,
        package_versions=package_versions,
        de_package_version=DE_PACKAGE_VERSION,
        blc_mark_version=BLC_MARK_VERSION,
        analysis_status=analysis_status,
    )
