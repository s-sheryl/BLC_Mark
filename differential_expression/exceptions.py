"""
Purpose:
    Define the exception hierarchy specific to Phase 3 (Differential
    Expression). Every exception here ultimately derives from the
    existing project-wide root, src.exceptions.BLCMarkError, so
    callers that already catch BLCMarkError continue to catch
    every DE failure without modification.

Responsibilities:
    - Provide one DE root exception, DEError, that every other
      DE-specific exception subclasses.
    - Provide named exceptions for each distinct DE failure category
      required by the Version 1 differential expression specification
      (input validation, configuration, comparison, statistical
      method, multiple testing, quality control, reproducibility,
      and result writing).

Scope:
    This module only defines exception classes. It performs no I/O,
    no validation logic, no statistical computation, and does not
    modify src/exceptions.py. It is purely additive to the existing
    root hierarchy.

Version:
    Phase 3 Version 1.0 DE exception hierarchy.
"""

from src.exceptions import BLCMarkError

DE_EXCEPTIONS_VERSION = "1.0"

__all__ = [
    "DE_EXCEPTIONS_VERSION",
    "DEError",
    "DEValidationError",
    "InvalidExpressionMatrixError",
    "InvalidMetadataError",
    "SampleMismatchError",
    "InsufficientReplicationError",
    "InvalidConfigurationError",
    "UnsupportedMethodError",
    "StatisticalMethodError",
    "MultipleTestingError",
    "QualityControlError",
    "ReproducibilityError",
    "ResultWritingError",
]


class DEError(BLCMarkError):
    """Root exception for every Phase 3 (Differential Expression) failure.

    Subclasses BLCMarkError so existing code that catches the
    project-wide root exception continues to catch DE failures
    without modification.
    """


class DEValidationError(DEError):
    """Raised for structural validation failures of DE inputs.

    This is the general validation-failure category used when a more
    specific subclass (InvalidExpressionMatrixError,
    InvalidMetadataError, SampleMismatchError) does not apply.
    """


class InvalidExpressionMatrixError(DEValidationError):
    """Raised when the expression matrix fails structural validation.

    Examples: missing file, unreadable file, missing gene-identifier
    column, no sample columns, duplicate gene identifiers, duplicate
    sample identifiers, non-numeric expression values.
    """


class InvalidMetadataError(DEValidationError):
    """Raised when the sample metadata fails structural validation.

    Examples: missing file, unreadable file, missing sample-identifier
    column, missing group column, duplicate sample identifiers,
    empty/missing group labels.
    """


class SampleMismatchError(DEValidationError):
    """Raised when expression and metadata samples cannot be
    reconciled sufficiently to proceed.

    This is raised specifically for the zero-overlap case described
    by the specification (Section 3.4): if expression samples and
    metadata samples cannot be matched at all, differential
    expression analysis must fail explicitly rather than silently
    proceeding with an implicit subset.
    """


class InsufficientReplicationError(DEValidationError):
    """Raised when a comparison group does not contain enough
    biological replicates for the configured statistical method,
    per specification Section 4.3.
    """


class InvalidConfigurationError(DEError):
    """Raised when the DE analysis configuration itself is
    scientifically invalid or incomplete.

    Examples: identical reference/comparison groups, an unsupported
    statistical method or multiple-testing method, an invalid
    significance threshold, or missing required configuration values.
    Per specification Section 10.2, configuration errors are
    distinguished from failures caused by the input dataset (which
    raise DEValidationError subclasses instead).
    """


class UnsupportedMethodError(DEError):
    """Raised when a configured statistical method is not available
    for execution in the current environment or is incompatible with
    the declared expression representation.

    This exception is raised instead of silently substituting a
    different method or silently proceeding, per specification
    Sections 5.4 and 10.3.
    """


class StatisticalMethodError(DEError):
    """Raised when the statistical framework itself reports a
    dataset-level failure while executing the configured method.

    Per specification Section 10.3, an individual gene that cannot
    produce a valid statistical result is retained in the output with
    missing values (this is not an error); this exception is reserved
    for failures that prevent the statistical method from running at
    all for the dataset.
    """


class MultipleTestingError(DEError):
    """Raised when multiple-testing correction cannot be performed
    or completed as configured.
    """


class QualityControlError(DEError):
    """Raised when a required pre-analysis quality-control condition
    fails, per specification Section 8.1.

    Note this is intentionally distinct from
    src.exceptions.QualityControlError, which belongs to the
    project-wide, dataset-download-oriented QC concept. This class
    lives in the DE package and is never confused with that class by
    virtue of its distinct module path.
    """


class ReproducibilityError(DEError):
    """Raised when reproducibility metadata cannot be constructed,
    e.g. an input file cannot be hashed or a required path is
    missing.
    """


class ResultWritingError(DEError):
    """Raised when DE outputs (results table, metadata, or QC report)
    cannot be safely written.

    Per specification Section 10.4, partially written output must
    never be presented as a completed analysis; raising this
    exception is how the analysis is marked unsuccessful.
    """


if __name__ == "__main__":
    assert issubclass(DEError, BLCMarkError)
    assert issubclass(DEValidationError, DEError)
    assert issubclass(InvalidExpressionMatrixError, DEValidationError)
    assert issubclass(InvalidMetadataError, DEValidationError)
    assert issubclass(SampleMismatchError, DEValidationError)
    assert issubclass(InsufficientReplicationError, DEValidationError)
    assert issubclass(InvalidConfigurationError, DEError)
    assert issubclass(UnsupportedMethodError, DEError)
    assert issubclass(StatisticalMethodError, DEError)
    assert issubclass(MultipleTestingError, DEError)
    assert issubclass(QualityControlError, DEError)
    assert issubclass(ReproducibilityError, DEError)
    assert issubclass(ResultWritingError, DEError)

    print(
        f"BLC Mark differential_expression.exceptions "
        f"(version {DE_EXCEPTIONS_VERSION}) hierarchy verified."
    )