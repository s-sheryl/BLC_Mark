"""
Purpose:
    Define the centralized exception hierarchy shared across every
    BLC Mark module.

Responsibilities:
    - Provide a single root exception, BLCMarkError, that every
      project-specific exception inherits from.
    - Group related failures under stable, meaningful parent exceptions
      so calling code can catch failures at the appropriate level.

Scope:
    This module only defines exception classes. It performs no I/O,
    creates no files, reads no configuration, uses no logging, makes
    no network requests, validates nothing, downloads nothing, and
    executes no project logic on import.

Version:
    This is the frozen Version 1.0 exception hierarchy. Its class names
    and inheritance structure are intended to remain stable so
    dependent modules can rely on them without future breaking changes.
"""

EXCEPTIONS_VERSION = "1.0"

__all__ = [
    "EXCEPTIONS_VERSION",
    "BLCMarkError",
    "ConfigurationError",
    "DownloadError",
    "DownloadVerificationError",
    "ChecksumMismatchError",
    "DatasetError",
    "MetadataError",
    "PreprocessingError",
    "QualityControlError",
    "XenaError",
]


class BLCMarkError(Exception):
    """Root exception for every BLC Mark-specific failure."""


class ConfigurationError(BLCMarkError):
    """Raised when project configuration is missing or invalid."""


class DownloadError(BLCMarkError):
    """Raised for failures related to downloading files."""


class DownloadVerificationError(DownloadError):
    """Raised when a downloaded file fails post-download verification."""


class ChecksumMismatchError(DownloadVerificationError):
    """Raised when a downloaded file's checksum does not match expectations."""


class DatasetError(BLCMarkError):
    """Raised for failures related to dataset management."""


class MetadataError(BLCMarkError):
    """Raised for failures related to metadata management."""


class PreprocessingError(BLCMarkError):
    """Raised for failures during dataset preprocessing."""


class QualityControlError(BLCMarkError):
    """Raised for failures during dataset quality control evaluation."""


class XenaError(BLCMarkError):
    """Raised for failures related to communicating with a Xena data hub."""


if __name__ == "__main__":
    assert issubclass(ConfigurationError, BLCMarkError)

    assert issubclass(DownloadError, BLCMarkError)
    assert issubclass(DownloadVerificationError, DownloadError)
    assert issubclass(ChecksumMismatchError, DownloadVerificationError)

    assert issubclass(DatasetError, BLCMarkError)
    assert issubclass(MetadataError, BLCMarkError)
    assert issubclass(PreprocessingError, BLCMarkError)
    assert issubclass(QualityControlError, BLCMarkError)
    assert issubclass(XenaError, BLCMarkError)
    assert issubclass(BLCMarkError, Exception)

    print(
        f"BLC Mark exceptions.py "
        f"(version {EXCEPTIONS_VERSION}) hierarchy verified."
    )
