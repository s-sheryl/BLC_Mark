"""
Purpose:
    Verify that files already downloaded into BLC Mark's
    data/raw/ directory are structurally fit to enter the
    preprocessing pipeline -- present, non-empty, an expected file
    type, and accompanied by complete metadata.

Responsibilities:
    - Confirm a downloaded file actually exists on disk.
    - Confirm it isn't empty.
    - Confirm its extension matches what the pipeline expects to
      read.
    - Confirm its checksum is well-formed (a plausible SHA-256, not
      that it matches any particular expected value -- that
      comparison belongs to whoever recorded the expected checksum).
    - Confirm its DatasetMetadata record is complete enough to be
      trusted downstream.
    - Aggregate all of the above into a single ValidationResult per
      dataset, and summarize many results into a report.

Scope:
    This module validates structure, not biology. It never opens a
    file to inspect gene expression values, sample identifiers, or
    any other biological content -- that belongs to preprocessing,
    which runs only after a dataset has passed here. It does not
    download files (download_manager.py, xena_client.py), and it does
    not compute checksums itself (metadata_manager.py does that) --
    it only checks that a checksum it's given looks well-formed.
"""

import hashlib
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path

from metadata_manager import DatasetMetadata

VALIDATION_MANAGER_VERSION = "1.0"

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_ALLOWED_EXTENSIONS = (
    ".txt",
    ".tsv",
    ".csv",
    ".gz",
)

# Extensions the preprocessing pipeline is expected to be able to
# read. A file arriving with any other suffix is treated as a
# structural failure, not something preprocessing should have to
# discover on its own.

MINIMUM_FILE_SIZE_BYTES = 1

# A SHA-256 hex digest is always exactly 64 hexadecimal characters.
# This is a structural fact about the algorithm, not something that
# depends on any particular file's contents.
SHA256_HEX_LENGTH = 64

__all__ = [
    "VALIDATION_MANAGER_VERSION",
    "DEFAULT_ALLOWED_EXTENSIONS",
    "MINIMUM_FILE_SIZE_BYTES",
    "SHA256_HEX_LENGTH",
    "ValidationStatus",
    "ValidationResult",
    "ValidationManager",
]


class ValidationStatus(str, Enum):
    """
    Overall outcome of validating one dataset.

    WARNING is distinct from FAILED: a dataset in WARNING state has
    no structural errors (it's safe to hand to preprocessing) but has
    something worth a human's attention -- e.g. an unusually small
    file size. FAILED means at least one hard requirement was not
    met and the dataset should not proceed.
    """

    PENDING = "pending"
    PASSED = "passed"
    FAILED = "failed"
    WARNING = "warning"


@dataclass(frozen=True)
class ValidationResult:
    """
    Outcome of running every validation check against one dataset.

    Frozen for the same reason DatasetMetadata is frozen in
    metadata_manager.py: a validation result is a record of what was
    true at the moment validation ran, and should not be quietly
    edited afterward. errors and warnings are plain lists rather than
    a single message, because validate_dataset() is required to keep
    checking after the first failure -- a dataset can fail more than
    one way at once, and a researcher fixing it deserves the full
    picture in one pass.
    """

    dataset_id: str
    passed: bool
    status: ValidationStatus
    validation_time: datetime
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    notes: str | None = None


class ValidationManager:
    """
    Runs structural validation checks against downloaded datasets.

    Every check here answers "is this file safe to hand to
    preprocessing", never "is this file biologically correct".
    """

    VERSION = VALIDATION_MANAGER_VERSION

    def __init__(
        self,
        allowed_extensions: tuple[str, ...] = DEFAULT_ALLOWED_EXTENSIONS,
    ) -> None:

        """
        Store validation configuration. No file is touched here --
        validation only happens when one of the validate_* methods is
        called explicitly.

        Args:
            allowed_extensions: File suffixes considered acceptable
                for a downloaded dataset. Defaults to
                DEFAULT_ALLOWED_EXTENSIONS.
        """
        self.allowed_extensions = allowed_extensions

    def validate_file_exists(self, file_path: Path) -> bool:
        """
        Confirm a file exists at the given path.

        Args:
            file_path: Path to the downloaded file.

        Returns:
            True if the file exists.

        Raises:
            FileNotFoundError: If no file exists at `file_path`.
        """
        if not file_path.exists():
            raise FileNotFoundError(f"File does not exist: {file_path}")
        return True

    def validate_file_size(self, file_path: Path) -> bool:
        """
        Confirm a file meets the minimum size requirement.

        A zero-byte file most often indicates a failed or interrupted
        download that nonetheless left a file behind -- this catches
        that case before it reaches preprocessing.

        Args:
            file_path: Path to the downloaded file.

        Returns:
            True if the file's size is >= MINIMUM_FILE_SIZE_BYTES.

        Raises:
            ValueError: If the file is empty (or below the minimum).
        """
        file_size = file_path.stat().st_size
        if file_size < MINIMUM_FILE_SIZE_BYTES:
            raise ValueError(
                f"File is empty or below the minimum size "
                f"({file_size} bytes < {MINIMUM_FILE_SIZE_BYTES}): {file_path}"
            )
        return True

    def validate_extension(self, file_path: Path) -> bool:
        """
        Confirm a file's suffix is one preprocessing knows how to
        read.

        Args:
            file_path: Path to the downloaded file.

        Returns:
            True if the suffix is in self.allowed_extensions.

        Raises:
            ValueError: If the suffix is not recognized. `.gz` files
                are accepted at their outer suffix only -- this
                method does not look inside a compressed file to
                validate what it contains.
        """
        suffix = file_path.suffix.lower()
        if suffix not in self.allowed_extensions:
            raise ValueError(
                f"Unsupported file extension '{suffix}' for {file_path}. "
                f"Allowed extensions: {self.allowed_extensions}"
            )
        return True

    def validate_checksum(self, checksum: str) -> bool:
        """
        Confirm a checksum string is a well-formed SHA-256 hex digest.

        This checks *shape*, not correctness -- it does not
        recompute the checksum or compare it against an expected
        value. Confirming a checksum matches the file it claims to
        describe is metadata_manager.py's responsibility at the point
        the checksum is generated.

        Args:
            checksum: The checksum string to check.

        Returns:
            True if `checksum` is present, exactly 64 characters, and
            valid hexadecimal.

        Raises:
            ValueError: If the checksum is missing, the wrong length,
                or contains non-hexadecimal characters.
        """
        if not checksum:
            raise ValueError("Checksum is missing or empty.")

        if len(checksum) != SHA256_HEX_LENGTH:
            raise ValueError(
                f"Checksum has length {len(checksum)}, expected "
                f"{SHA256_HEX_LENGTH} for a SHA-256 hex digest: {checksum!r}"
            )

        try:
            int(checksum, 16)
        except ValueError as error:
            raise ValueError(
                f"Checksum is not valid hexadecimal: {checksum!r}"
            ) from error

        return True

    def validate_metadata(self, record: DatasetMetadata) -> bool:
        """
        Confirm a DatasetMetadata record carries the minimum fields
        required for a dataset to be trustworthy downstream.

        This checks that required fields are populated -- it does not
        re-derive or re-check their values (e.g. it does not
        recompute the checksum or re-stat the file; validate_checksum()
        and validate_file_size() handle those independently).

        Args:
            record: The metadata record to check.

        Returns:
            True if dataset_id, sha256_checksum, local_path, and
            file_size_bytes are all present and non-empty.

        Raises:
            ValueError: If any required field is missing, empty, or
                otherwise clearly invalid (e.g. a negative file size).
        """
        missing_fields: list[str] = []

        if not record.dataset_id:
            missing_fields.append("dataset_id")
        if not record.sha256_checksum:
            missing_fields.append("sha256_checksum")
        if record.local_path is None:
            missing_fields.append("local_path")
        if record.file_size_bytes is None or record.file_size_bytes < 0:
            missing_fields.append("file_size_bytes")

        if missing_fields:
            raise ValueError(
                f"DatasetMetadata for '{record.dataset_id or '<unknown>'}' "
                f"is missing required field(s): {', '.join(missing_fields)}."
            )

        return True

    def validate_dataset(self, record: DatasetMetadata) -> ValidationResult:
        """
        Run every validation check against one dataset and collect
        the full outcome.

        Every check below runs regardless of whether an earlier one
        failed -- a dataset with three problems should be reported
        with all three at once, not fixed and re-run three separate
        times to discover them one at a time.

        Args:
            record: The metadata record describing the dataset to
                validate.

        Returns:
            A ValidationResult summarizing every check performed.
        """
        errors: list[str] = []
        warnings: list[str] = []

        try:
            self.validate_metadata(record)
        except ValueError as error:
            errors.append(str(error))

        # The remaining checks depend on local_path/checksum existing,
        # but we still attempt each independently so a missing field
        # doesn't hide problems with the fields that ARE present.
        if record.local_path is not None:
            try:
                self.validate_file_exists(record.local_path)
            except FileNotFoundError as error:
                errors.append(str(error))
            else:
                try:
                    self.validate_file_size(record.local_path)
                except ValueError as error:
                    errors.append(str(error))

                try:
                    self.validate_extension(record.local_path)
                except ValueError as error:
                    errors.append(str(error))
        else:
            errors.append("Cannot check file existence, size, or extension: local_path is missing.")

        if record.sha256_checksum:
            try:
                self.validate_checksum(record.sha256_checksum)
            except ValueError as error:
                errors.append(str(error))
        else:
            errors.append("Cannot validate checksum: sha256_checksum is missing.")

        # Soft, non-fatal observation: a very small file is often
        # legitimate (e.g. a small clinical annotation file) but is
        # also a common symptom of a truncated download, so it's
        # surfaced as a warning rather than silently ignored.
        if record.file_size_bytes is not None and 0 < record.file_size_bytes < 1024:
            warnings.append(
                f"File size is unusually small ({record.file_size_bytes} bytes). "
                "Confirm this is expected rather than a truncated download."
            )

        if errors:
            status = ValidationStatus.FAILED
        elif warnings:
            status = ValidationStatus.WARNING
        else:
            status = ValidationStatus.PASSED

        return ValidationResult(
            dataset_id=record.dataset_id,
            passed=not errors,
            status=status,
            validation_time=datetime.now(timezone.utc),
            errors=errors,
            warnings=warnings,
        )

    def build_validation_report(self, results: list[ValidationResult]) -> dict:
        """
        Summarize a batch of ValidationResult objects.

        Pure aggregation over results already computed -- no file
        access happens here, so this can be called freely (e.g. from
        a dashboard) without re-running any checks.

        Args:
            results: Validation results to summarize, typically the
                output of several validate_dataset() calls.

        Returns:
            A dictionary with:
                - "total": count of results
                - "passed": count with status == PASSED
                - "failed": count with status == FAILED
                - "warnings": count with status == WARNING
                - "datasets": list of {"dataset_id", "status"} entries
        """
        return {
            "total": len(results),
            "passed": sum(1 for r in results if r.status == ValidationStatus.PASSED),
            "failed": sum(1 for r in results if r.status == ValidationStatus.FAILED),
            "pending": sum(1 for r in results if r.status == ValidationStatus.PENDING),
            "warnings": sum(1 for r in results if r.status == ValidationStatus.WARNING),
            "datasets": [
                {"dataset_id": r.dataset_id, "status": r.status.value} for r in results
            ],
        }

    def render_validation_summary(self, report: dict) -> str:
        """
        Render a validation report as a human-readable, multi-line
        string suitable for console output.

        Kept separate from build_validation_report() so the same
        report data could later feed a dashboard widget instead of
        text, without recomputing anything.

        Args:
            report: The dictionary produced by build_validation_report().

        Returns:
            A formatted multi-line summary string.
        """
        lines = [
            "Validation Summary",
            "-------------------",
            f"Total datasets:  {report['total']}",
            f"Passed:          {report['passed']}",
            f"Failed:          {report['failed']}",
            f"Warnings:        {report['warnings']}",
            "",
            "Per-dataset status:",
        ]

        for dataset in report["datasets"]:
            lines.append(f"  [{dataset['status'].upper():>7}] {dataset['dataset_id']}")

        return "\n".join(lines)

    # -- Export: architecture only, not yet implemented --

    def export_report_json(self, report: dict, destination: Path) -> None:
        """
        Export a validation report to a JSON file.

        TODO(implementation): serialize `report` (already plain
        dicts/lists/strings, so this should be a direct json.dump)
        and write it to `destination`.

        Raises:
            NotImplementedError: Always, until implemented.
        """
        raise NotImplementedError(
            "ValidationManager.export_report_json() is not implemented "
            f"yet. It will write the validation report to {destination}."
        )

    def export_report_csv(self, results: list[ValidationResult], destination: Path) -> None:
        """
        Export per-dataset validation results to a CSV file.

        TODO(implementation): flatten each ValidationResult to a row
        (joining errors/warnings lists into delimited strings) and
        write to `destination` using the standard library csv module.

        Raises:
            NotImplementedError: Always, until implemented.
        """
        raise NotImplementedError(
            "ValidationManager.export_report_csv() is not implemented "
            f"yet. It will write {len(results)} result(s) to {destination}."
        )

# ==========================================================
# Future Validation Pipeline
#
# ValidationManager
#        │
#        ▼
# PreprocessingManager
#        │
#        ▼
# Differential Expression
#        │
#        ▼
# Biomarker Discovery
#
# Validation guarantees that only structurally valid files
# enter the computational pipeline.
# ==========================================================

if __name__ == "__main__":
    import tempfile

    from dataset_registry import CancerType, DataType

    print("BLC Mark Validation Manager")
    print(f"Version: {VALIDATION_MANAGER_VERSION}\n")

    manager = ValidationManager()

    # Sanity check only -- this creates one small synthetic text file
    # (not biological data of any kind) purely to exercise every
    # validation method end to end. No network access, no TCGA
    # download, nothing biological is inspected.
    with tempfile.TemporaryDirectory() as scratch_dir:
        scratch_file = Path(scratch_dir) / "synthetic_test_file.txt"
        # Padded past the small-file warning threshold (1024 bytes) so
        # this sanity check demonstrates a clean PASS. The content
        # itself remains meaningless placeholder text, not biology.
        scratch_file.write_text("BLC Mark validation manager sanity check.\n" * 30)

        synthetic_checksum = hashlib.sha256(scratch_file.read_bytes()).hexdigest()

        synthetic_record = DatasetMetadata(
            dataset_id="SYNTHETIC.TEST.FILE",
            cancer_type=CancerType.BRCA,
            data_type=DataType.GENE_EXPRESSION_RNASEQ,
            local_path=scratch_file,
            sha256_checksum=synthetic_checksum,
            file_size_bytes=scratch_file.stat().st_size,
            notes="Synthetic file used only to verify validation logic.",
        )

        result = manager.validate_dataset(synthetic_record)

        print(f"Synthetic dataset: {result.dataset_id}")
        print(f"Status: {result.status.value.upper()}")
        print(f"Passed: {result.passed}")
        if result.errors:
            print(f"Errors: {result.errors}")
        if result.warnings:
            print(f"Warnings: {result.warnings}")

        report = manager.build_validation_report([result])
        print()
        print(manager.render_validation_summary(report))

       