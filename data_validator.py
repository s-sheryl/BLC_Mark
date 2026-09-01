"""
Purpose:
    Manage the catalog of transcriptomic datasets available to
    BLC Mark before they enter preprocessing -- registering a
    dataset's identity and file locations, verifying its files are
    present and readable, and persisting that record as JSON.

Responsibilities:
    - Register a new dataset's metadata under a unique dataset_id.
    - Read, list, and remove registered dataset records.
    - Verify a registered dataset's referenced files actually exist
      and are readable, and that its metadata is complete.
    - Persist and load dataset records as JSON under a configurable
      metadata directory.

Scope:
    This module never downloads data -- that is download_manager.py
    and xena_client.py's responsibility. It never inspects the
    contents of an expression matrix or interprets biology -- that is
    preprocessing_manager.py and quality_control.py's responsibility.
    This module only manages the catalog entry describing a dataset:
    where its files are, what it is, and whether that record is
    complete and consistent with what's on disk.
"""

import csv
import hashlib
import json
import logging
import re
from dataclasses import asdict, dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

DATASET_MANAGER_VERSION = "1.1"

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIRECTORY = PROJECT_ROOT / "data" / "processed" / "datasets"

# Streaming chunk size for SHA-256 hashing -- large expression matrices
# should never need to be read fully into memory just to fingerprint
# them. Matches the chunk size used elsewhere in the project (see
# metadata_manager.py and download_manager.py) for consistency.
CHECKSUM_CHUNK_BYTES = 1024 * 1024  # 1 MB

REQUIRED_METADATA_FIELDS: tuple[str, ...] = (
    "dataset_id",
    "dataset_name",
    "source",
    "cancer_type",
    "organism",
    "platform",
    "expression_file",
)

__all__ = [
    "DATASET_MANAGER_VERSION",
    "DEFAULT_METADATA_DIRECTORY",
    "DatasetMetadata",
    "DatasetManager",
]


@dataclass(frozen=True)
class DatasetMetadata:
    """
    Catalog record describing one transcriptomic dataset.

    This record describes a dataset's identity and file locations,
    not its downloaded bytes or biological content. Frozen so a
    record can never be edited in place -- updates go through
    DatasetManager, which produces a new record via
    dataclasses.replace() and persists it, preserving the invariant
    that whatever is on disk always matches an object that was
    deliberately saved rather than mutated ad hoc.

    Attributes:
        dataset_id: Unique identifier for this dataset within the
            catalog. Enforced unique by DatasetManager.register_dataset().
        dataset_name: Human-readable name for the dataset.
        source: Originating repository, e.g. "TCGA", "GEO". Not
            restricted to a fixed set of values -- any string
            describing the source is accepted, keeping this module
            usable for datasets outside the sources known today.
        cancer_type: Cancer type the dataset represents. Stored as a
            plain string rather than a fixed enum so this module
            never has to be modified to support a new cancer type.
        organism: Source organism, e.g. "Homo sapiens".
        platform: Sequencing or array platform used to generate the
            data.
        expression_file: Path to the expression matrix file.
        metadata_file: Path to an accompanying sample metadata file,
            if one exists.
        number_of_samples: Sample count, if known at registration
            time. None if not yet determined.
        number_of_genes: Gene count, if known at registration time.
            None if not yet determined.
        date_added: Timestamp when this dataset was registered.
        version: Version label for this dataset record, allowing the
            same dataset_id to be re-registered under a new version
            without losing the concept of "which version is this".
        notes: Free-text notes about the dataset.
        expression_file_sha256: SHA-256 hex digest of expression_file
            as it existed at registration time. Computed automatically
            by DatasetManager.register_dataset(); None for records
            created before this field existed. Used to detect
            accidental modification and to verify reproducibility of
            downloaded datasets.
        metadata_file_sha256: SHA-256 hex digest of metadata_file as
            it existed at registration time, if metadata_file was
            supplied. None if no metadata_file was given, or for
            records created before this field existed.
        expression_file_size: Size of expression_file in bytes at
            registration time. Used alongside expression_file_sha256
            to detect corrupted or truncated downloads.
        metadata_file_size: Size of metadata_file in bytes at
            registration time, if metadata_file was supplied.
    """

    dataset_id: str
    dataset_name: str
    source: str
    cancer_type: str
    organism: str
    platform: str
    expression_file: Path

    metadata_file: Path | None = None
    number_of_samples: int | None = None
    number_of_genes: int | None = None
    date_added: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    version: str = "1.0"
    notes: str | None = None

    expression_file_sha256: str | None = None
    metadata_file_sha256: str | None = None
    expression_file_size: int | None = None
    metadata_file_size: int | None = None


class DatasetManager:
    """
    Registers, persists, and validates dataset catalog records.

    This class owns the JSON files under its configured metadata
    directory -- one file per dataset_id -- and is the single place
    in the project that reads or writes them. It never downloads a
    file, and it never opens an expression matrix to read its
    contents; validate_registered_dataset() confirms files exist and
    are readable, nothing more.
    """

    def __init__(self, metadata_directory: Path = DEFAULT_METADATA_DIRECTORY) -> None:
        """
        Store where dataset records are persisted. No file or
        directory is created here -- that happens lazily the first
        time save_metadata() is called.

        Args:
            metadata_directory: Directory that will hold one JSON
                file per registered dataset.
        """
        self.metadata_directory = metadata_directory

    def _metadata_path_for(self, dataset_id: str) -> Path:
        """
        Compute the JSON file path for a given dataset_id.

        Args:
            dataset_id: The dataset identifier.

        Returns:
            Path to the dataset's metadata JSON file.
        """
        return self.metadata_directory / f"{dataset_id}.json"

    def dataset_exists(self, dataset_id: str) -> bool:
        """
        Check whether a dataset is already registered.

        Args:
            dataset_id: The dataset identifier to check.

        Returns:
            True if a metadata record exists for this dataset_id.
        """
        return self._metadata_path_for(dataset_id).exists()

    def register_dataset(
        self,
        dataset_id: str,
        dataset_name: str,
        source: str,
        cancer_type: str,
        organism: str,
        platform: str,
        expression_file: Path,
        *,
        metadata_file: Path | None = None,
        number_of_samples: int | None = None,
        number_of_genes: int | None = None,
        version: str = "1.0",
        notes: str | None = None,
    ) -> DatasetMetadata:
        """
        Register a new dataset in the catalog.

        Args:
            dataset_id: Unique identifier for the dataset. Must not
                already be registered.
            dataset_name: Human-readable name for the dataset.
            source: Originating repository, e.g. "TCGA", "GEO".
            cancer_type: Cancer type the dataset represents.
            organism: Source organism.
            platform: Sequencing or array platform.
            expression_file: Path to the expression matrix file.
            metadata_file: Path to an accompanying sample metadata
                file, if one exists.
            number_of_samples: Sample count, if known.
            number_of_genes: Gene count, if known.
            version: Version label for this dataset record.
            notes: Free-text notes about the dataset.

        Returns:
            The newly created DatasetMetadata record, including
            computed SHA-256 fingerprints and file sizes for
            expression_file (and metadata_file, if supplied).

        Raises:
            ValueError: If dataset_id is already registered.
            FileNotFoundError: If expression_file does not exist, or
                metadata_file is supplied but does not exist.
            OSError: If expression_file or metadata_file exists but
                cannot be opened for reading.
        """
        if self.dataset_exists(dataset_id):
            raise ValueError(
                f"Dataset '{dataset_id}' is already registered. "
                "Use a different dataset_id or remove the existing "
                "registration first."
            )

        self._validate_files_before_registration(
            dataset_id=dataset_id,
            expression_file=expression_file,
            metadata_file=metadata_file,
        )

        expression_file_sha256 = _calculate_sha256(expression_file)
        expression_file_size = _get_file_size(expression_file)

        metadata_file_sha256 = (
            _calculate_sha256(metadata_file) if metadata_file is not None else None
        )
        metadata_file_size = (
            _get_file_size(metadata_file) if metadata_file is not None else None
        )

        record = DatasetMetadata(
            dataset_id=dataset_id,
            dataset_name=dataset_name,
            source=source,
            cancer_type=cancer_type,
            organism=organism,
            platform=platform,
            expression_file=expression_file,
            metadata_file=metadata_file,
            number_of_samples=number_of_samples,
            number_of_genes=number_of_genes,
            version=version,
            notes=notes,
            expression_file_sha256=expression_file_sha256,
            metadata_file_sha256=metadata_file_sha256,
            expression_file_size=expression_file_size,
            metadata_file_size=metadata_file_size,
        )

        self.save_metadata(record)
        logger.info(
            "Registered dataset '%s' (%s). Expression file SHA-256: %s",
            dataset_id,
            dataset_name,
            expression_file_sha256,
        )

        return record

    def _validate_files_before_registration(
        self,
        dataset_id: str,
        expression_file: Path,
        metadata_file: Path | None,
    ) -> None:
        """
        Verify referenced files exist and are readable before a
        dataset is registered.

        Registration is refused rather than allowed to proceed with a
        dangling reference -- a catalog entry pointing at a file that
        was never there is worse than no entry at all, since it looks
        registered while being unusable.

        Args:
            dataset_id: The dataset identifier being registered, used
                only for error messages.
            expression_file: Path to the expression matrix file.
            metadata_file: Path to the sample metadata file, if
                supplied.

        Raises:
            FileNotFoundError: If expression_file does not exist, or
                metadata_file is supplied but does not exist.
            OSError: If expression_file or metadata_file exists but
                cannot be opened for reading.
        """
        if not expression_file.exists():
            raise FileNotFoundError(
                f"Cannot register dataset '{dataset_id}': expression file "
                f"does not exist: {expression_file}"
            )
        if not _is_readable(expression_file):
            raise OSError(
                f"Cannot register dataset '{dataset_id}': expression file "
                f"exists but is not readable: {expression_file}"
            )

        if metadata_file is not None:
            if not metadata_file.exists():
                raise FileNotFoundError(
                    f"Cannot register dataset '{dataset_id}': metadata file "
                    f"does not exist: {metadata_file}"
                )
            if not _is_readable(metadata_file):
                raise OSError(
                    f"Cannot register dataset '{dataset_id}': metadata file "
                    f"exists but is not readable: {metadata_file}"
                )

    def remove_dataset(self, dataset_id: str) -> None:
        """
        Remove a dataset's metadata record from the catalog.

        This deletes only the catalog record, not the underlying
        expression or metadata files it references -- removing a
        dataset from the catalog does not delete data from disk.

        Args:
            dataset_id: The dataset identifier to remove.

        Raises:
            FileNotFoundError: If no dataset is registered under this
                dataset_id.
        """
        metadata_path = self._metadata_path_for(dataset_id)

        if not metadata_path.exists():
            raise FileNotFoundError(f"No registered dataset found for '{dataset_id}'.")

        metadata_path.unlink()
        logger.info("Removed dataset '%s' from the catalog.", dataset_id)

    def get_dataset(self, dataset_id: str) -> DatasetMetadata:
        """
        Retrieve a registered dataset's metadata record.

        Args:
            dataset_id: The dataset identifier to retrieve.

        Returns:
            The dataset's DatasetMetadata record.

        Raises:
            FileNotFoundError: If no dataset is registered under this
                dataset_id.
        """
        return self.load_metadata(dataset_id)

    def list_datasets(self) -> list[DatasetMetadata]:
        """
        List every dataset currently registered in the catalog.

        Returns:
            A list of DatasetMetadata records, one per registered
            dataset. Empty if the metadata directory does not exist
            or contains no records.
        """
        if not self.metadata_directory.exists():
            return []

        records: list[DatasetMetadata] = []
        for metadata_path in sorted(self.metadata_directory.glob("*.json")):
            dataset_id = metadata_path.stem
            try:
                records.append(self.load_metadata(dataset_id))
            except (json.JSONDecodeError, KeyError, ValueError) as error:
                logger.warning(
                    "Skipping unreadable metadata file '%s': %s", metadata_path, error
                )

        return records

    def save_metadata(self, record: DatasetMetadata) -> None:
        """
        Persist a DatasetMetadata record as JSON.

        Args:
            record: The dataset record to save.
        """
        self.metadata_directory.mkdir(parents=True, exist_ok=True)

        serialized = asdict(record)
        serialized["expression_file"] = str(record.expression_file)
        serialized["metadata_file"] = (
            str(record.metadata_file) if record.metadata_file is not None else None
        )
        serialized["date_added"] = record.date_added.isoformat()

        metadata_path = self._metadata_path_for(record.dataset_id)
        with open(metadata_path, "w") as metadata_file_handle:
            json.dump(serialized, metadata_file_handle, indent=2)

        logger.debug("Saved metadata for dataset '%s' to %s", record.dataset_id, metadata_path)

    def load_metadata(self, dataset_id: str) -> DatasetMetadata:
        """
        Load a dataset's metadata record from disk.

        Args:
            dataset_id: The dataset identifier to load.

        Returns:
            The dataset's DatasetMetadata record.

        Raises:
            FileNotFoundError: If no metadata file exists for this
                dataset_id.
            KeyError: If the metadata file is missing a required
                field.
            ValueError: If the metadata file's JSON cannot be parsed
                or a date field cannot be interpreted.
        """
        metadata_path = self._metadata_path_for(dataset_id)

        if not metadata_path.exists():
            raise FileNotFoundError(f"No registered dataset found for '{dataset_id}'.")

        with open(metadata_path, "r") as metadata_file_handle:
            try:
                raw_data = json.load(metadata_file_handle)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"Metadata file for '{dataset_id}' is not valid JSON: {error}"
                ) from error

        try:
            return DatasetMetadata(
                dataset_id=raw_data["dataset_id"],
                dataset_name=raw_data["dataset_name"],
                source=raw_data["source"],
                cancer_type=raw_data["cancer_type"],
                organism=raw_data["organism"],
                platform=raw_data["platform"],
                expression_file=Path(raw_data["expression_file"]),
                metadata_file=(
                    Path(raw_data["metadata_file"])
                    if raw_data.get("metadata_file") is not None
                    else None
                ),
                number_of_samples=raw_data.get("number_of_samples"),
                number_of_genes=raw_data.get("number_of_genes"),
                date_added=datetime.fromisoformat(raw_data["date_added"]),
                version=raw_data.get("version", "1.0"),
                notes=raw_data.get("notes"),
                expression_file_sha256=raw_data.get("expression_file_sha256"),
                metadata_file_sha256=raw_data.get("metadata_file_sha256"),
                expression_file_size=raw_data.get("expression_file_size"),
                metadata_file_size=raw_data.get("metadata_file_size"),
            )
        except KeyError as error:
            raise KeyError(
                f"Metadata file for '{dataset_id}' is missing required field: {error}"
            ) from error

    def validate_registered_dataset(self, dataset_id: str) -> bool:
        """
        Verify a registered dataset's files exist, are readable, and
        its metadata record is complete.

        This checks structure and fingerprint consistency only -- file
        presence, readability, completeness of the catalog record,
        and (where a fingerprint was recorded) that the file on disk
        still matches the SHA-256 hash and size captured at
        registration time. It never opens expression_file to inspect
        its biological content; that is DataValidator's and
        QualityControl's responsibility.

        Args:
            dataset_id: The dataset identifier to validate.

        Returns:
            True if every check passes.

        Raises:
            FileNotFoundError: If no dataset is registered under this
                dataset_id, or if expression_file or metadata_file
                (when set) does not exist.
            ValueError: If the metadata record is missing a required
                field, a file exists but is not readable, or a file's
                current SHA-256 hash or size no longer matches the
                value recorded at registration time.
        """
        record = self.load_metadata(dataset_id)

        missing_fields = [
            field_name
            for field_name in REQUIRED_METADATA_FIELDS
            if not getattr(record, field_name, None)
        ]
        if missing_fields:
            raise ValueError(
                f"Dataset '{dataset_id}' has incomplete metadata. "
                f"Missing field(s): {missing_fields}"
            )

        if not record.expression_file.exists():
            raise FileNotFoundError(
                f"Expression file for dataset '{dataset_id}' does not exist: "
                f"{record.expression_file}"
            )

        if not _is_readable(record.expression_file):
            raise ValueError(
                f"Expression file for dataset '{dataset_id}' exists but is "
                f"not readable: {record.expression_file}"
            )

        self._verify_fingerprint(
            dataset_id=dataset_id,
            file_path=record.expression_file,
            file_label="expression file",
            expected_sha256=record.expression_file_sha256,
            expected_size=record.expression_file_size,
        )

        if record.metadata_file is not None:
            if not record.metadata_file.exists():
                raise FileNotFoundError(
                    f"Metadata file for dataset '{dataset_id}' does not exist: "
                    f"{record.metadata_file}"
                )
            if not _is_readable(record.metadata_file):
                raise ValueError(
                    f"Metadata file for dataset '{dataset_id}' exists but is "
                    f"not readable: {record.metadata_file}"
                )

            self._verify_fingerprint(
                dataset_id=dataset_id,
                file_path=record.metadata_file,
                file_label="metadata file",
                expected_sha256=record.metadata_file_sha256,
                expected_size=record.metadata_file_size,
            )

        logger.info("Dataset '%s' passed registration validation.", dataset_id)
        return True

    def _verify_fingerprint(
        self,
        dataset_id: str,
        file_path: Path,
        file_label: str,
        expected_sha256: str | None,
        expected_size: int | None,
    ) -> None:
        """
        Confirm a file's current SHA-256 hash and size still match the
        values recorded at registration time.

        Records created before fingerprinting existed have
        expected_sha256 and expected_size set to None -- those are
        skipped rather than failed, since there is nothing to compare
        against. This keeps validate_registered_dataset() backward
        compatible with catalogs registered under an earlier version
        of this module.

        Args:
            dataset_id: The dataset identifier, used only for error
                messages.
            file_path: The file to re-fingerprint.
            file_label: Human-readable label for the file (e.g.
                "expression file"), used only for error messages.
            expected_sha256: The SHA-256 hex digest recorded at
                registration time, or None if not recorded.
            expected_size: The file size in bytes recorded at
                registration time, or None if not recorded.

        Raises:
            ValueError: If the file's current hash or size does not
                match the recorded value.
        """
        if expected_size is not None:
            current_size = _get_file_size(file_path)
            if current_size != expected_size:
                raise ValueError(
                    f"Dataset '{dataset_id}' {file_label} size has changed "
                    f"since registration: expected {expected_size} bytes, "
                    f"found {current_size} bytes ({file_path}). The file may "
                    "have been corrupted, truncated, or replaced."
                )

        if expected_sha256 is not None:
            current_sha256 = _calculate_sha256(file_path)
            if current_sha256 != expected_sha256:
                raise ValueError(
                    f"Dataset '{dataset_id}' {file_label} SHA-256 no longer "
                    f"matches the value recorded at registration: expected "
                    f"{expected_sha256}, found {current_sha256} ({file_path}). "
                    "The file may have been modified since it was registered."
                )

    def update_dataset(self, dataset_id: str, **fields) -> DatasetMetadata:
        """
        Update fields on an existing dataset record and persist the
        result.

        Since DatasetMetadata is frozen, this loads the existing
        record, applies the given field updates via
        dataclasses.replace(), and saves the resulting new record --
        the on-disk file is fully overwritten, not patched in place.

        Args:
            dataset_id: The dataset identifier to update.
            **fields: Field names and new values to apply.

        Returns:
            The updated DatasetMetadata record.

        Raises:
            FileNotFoundError: If no dataset is registered under this
                dataset_id.
            TypeError: If a field name in **fields does not exist on
                DatasetMetadata.
        """
        existing_record = self.load_metadata(dataset_id)
        updated_record = replace(existing_record, **fields)
        self.save_metadata(updated_record)

        logger.info("Updated dataset '%s': %s", dataset_id, list(fields.keys()))
        return updated_record

    def get_dataset_summary(self, dataset_id: str) -> dict[str, Any]:
        """
        Produce a human-readable summary of a registered dataset's
        stored metadata.

        This reads only the catalog record -- it never opens
        expression_file to inspect its contents, so it stays cheap to
        call regardless of dataset size.

        Args:
            dataset_id: The dataset identifier to summarize.

        Returns:
            A dictionary with the dataset's name, source, cancer
            type, organism, platform, sample/gene counts, version,
            registration date, file paths, file sizes, and SHA-256
            fingerprints.

        Raises:
            FileNotFoundError: If no dataset is registered under this
                dataset_id.
        """
        record = self.load_metadata(dataset_id)

        return {
            "dataset_id": record.dataset_id,
            "dataset_name": record.dataset_name,
            "source": record.source,
            "cancer_type": record.cancer_type,
            "organism": record.organism,
            "platform": record.platform,
            "number_of_samples": record.number_of_samples,
            "number_of_genes": record.number_of_genes,
            "version": record.version,
            "date_added": record.date_added.isoformat(),
            "expression_file": str(record.expression_file),
            "metadata_file": (
                str(record.metadata_file) if record.metadata_file is not None else None
            ),
            "expression_file_size": record.expression_file_size,
            "metadata_file_size": record.metadata_file_size,
            "expression_file_sha256": record.expression_file_sha256,
            "metadata_file_sha256": record.metadata_file_sha256,
        }

    def get_catalog_statistics(self) -> dict[str, Any]:
        """
        Compute summary statistics across the entire dataset catalog.

        Pure aggregation over already-registered records -- no file
        access beyond the metadata JSON files themselves.

        Returns:
            A dictionary with:
                - "total_datasets": count of registered datasets
                - "datasets_per_source": dict of source -> count
                - "datasets_per_cancer_type": dict of cancer_type -> count
                - "datasets_with_metadata": count with metadata_file set
                - "datasets_without_metadata": count with metadata_file None
                - "oldest_registration": ISO timestamp of the earliest
                  date_added, or None if the catalog is empty
                - "newest_registration": ISO timestamp of the latest
                  date_added, or None if the catalog is empty
        """
        records = self.list_datasets()

        datasets_per_source: dict[str, int] = {}
        datasets_per_cancer_type: dict[str, int] = {}
        datasets_with_metadata = 0

        for record in records:
            datasets_per_source[record.source] = datasets_per_source.get(record.source, 0) + 1
            datasets_per_cancer_type[record.cancer_type] = (
                datasets_per_cancer_type.get(record.cancer_type, 0) + 1
            )
            if record.metadata_file is not None:
                datasets_with_metadata += 1

        oldest_registration = min((r.date_added for r in records), default=None)
        newest_registration = max((r.date_added for r in records), default=None)

        return {
            "total_datasets": len(records),
            "datasets_per_source": datasets_per_source,
            "datasets_per_cancer_type": datasets_per_cancer_type,
            "datasets_with_metadata": datasets_with_metadata,
            "datasets_without_metadata": len(records) - datasets_with_metadata,
            "oldest_registration": (
                oldest_registration.isoformat() if oldest_registration is not None else None
            ),
            "newest_registration": (
                newest_registration.isoformat() if newest_registration is not None else None
            ),
        }

    def find_datasets(
        self,
        *,
        source: str | None = None,
        cancer_type: str | None = None,
        organism: str | None = None,
        platform: str | None = None,
        version: str | None = None,
    ) -> list[DatasetMetadata]:
        """
        Search the catalog for datasets matching the given filters.

        Every provided filter must match exactly (case-sensitive) for
        a record to be included; filters left as None are ignored.
        Passing no filters returns every registered dataset, the same
        as list_datasets().

        Args:
            source: Exact source to match, e.g. "TCGA".
            cancer_type: Exact cancer type to match.
            organism: Exact organism to match.
            platform: Exact platform to match.
            version: Exact version label to match.

        Returns:
            A list of matching DatasetMetadata records.
        """
        filters = {
            "source": source,
            "cancer_type": cancer_type,
            "organism": organism,
            "platform": platform,
            "version": version,
        }
        active_filters = {name: value for name, value in filters.items() if value is not None}

        return [
            record
            for record in self.list_datasets()
            if all(getattr(record, name) == value for name, value in active_filters.items())
        ]

    def export_catalog(self, destination: Path, file_format: str = "json") -> None:
        """
        Export the entire catalog to a single file.

        Args:
            destination: Path to write the export to.
            file_format: Either "json" or "csv" (case-insensitive).

        Raises:
            ValueError: If file_format is not "json" or "csv".
        """
        normalized_format = file_format.strip().lower()

        if normalized_format == "json":
            self._export_catalog_json(destination)
        elif normalized_format == "csv":
            self._export_catalog_csv(destination)
        else:
            raise ValueError(
                f"Unsupported export format '{file_format}'. Expected 'json' or 'csv'."
            )

    def _export_catalog_json(self, destination: Path) -> None:
        """
        Export the catalog as a single JSON array of dataset records.

        Args:
            destination: Path to write the JSON export to.
        """
        records = self.list_datasets()
        serializable_records = [self.get_dataset_summary(r.dataset_id) for r in records]

        destination.parent.mkdir(parents=True, exist_ok=True)
        with open(destination, "w") as export_file:
            json.dump(serializable_records, export_file, indent=2)

        logger.info("Exported %d dataset(s) to JSON: %s", len(records), destination)

    def _export_catalog_csv(self, destination: Path) -> None:
        """
        Export the catalog as a CSV file, one row per dataset.

        Args:
            destination: Path to write the CSV export to.
        """
        records = self.list_datasets()
        summaries = [self.get_dataset_summary(r.dataset_id) for r in records]

        destination.parent.mkdir(parents=True, exist_ok=True)

        if not summaries:
            destination.write_text("")
            logger.info("Exported empty catalog to CSV: %s", destination)
            return

        fieldnames = list(summaries[0].keys())
        with open(destination, "w", newline="") as export_file:
            writer = csv.DictWriter(export_file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(summaries)

        logger.info("Exported %d dataset(s) to CSV: %s", len(summaries), destination)

    def list_versions(self, dataset_id: str) -> list[DatasetMetadata]:
        """
        List every registered dataset that is a version of the given
        base dataset_id.

        Versioning here is a naming convention, not a stored
        relationship: a dataset_id is considered a version of
        `dataset_id` if it equals `dataset_id` exactly, or extends it
        with a "_v<digits>" suffix (e.g. both "TCGA_BRCA" and
        "TCGA_BRCA_v2" are versions of base "TCGA_BRCA"). This keeps
        multi-version support additive -- existing single-version
        datasets (registered with a plain dataset_id) work exactly as
        before, and nothing about DatasetMetadata's stored fields
        needed to change to support it.

        Args:
            dataset_id: The base dataset identifier to find versions
                of.

        Returns:
            Matching DatasetMetadata records, sorted by date_added
            ascending (oldest first).
        """
        version_pattern = re.compile(rf"^{re.escape(dataset_id)}(_v\d+)?$")

        matching_records = [
            record for record in self.list_datasets() if version_pattern.match(record.dataset_id)
        ]

        return sorted(matching_records, key=lambda record: record.date_added)


def _calculate_sha256(file_path: Path) -> str:
    """
    Compute the SHA-256 checksum of a file, streamed in fixed-size
    chunks so file size never dictates memory usage.

    Args:
        file_path: Path to an existing, readable file.

    Returns:
        The SHA-256 checksum as a lowercase hex digest.
    """
    digest = hashlib.sha256()
    with open(file_path, "rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(CHECKSUM_CHUNK_BYTES), b""):
            digest.update(chunk)

    return digest.hexdigest()


def _get_file_size(file_path: Path) -> int:
    """
    Get a file's size in bytes.

    Args:
        file_path: Path to an existing file.

    Returns:
        File size in bytes.
    """
    return file_path.stat().st_size


def _is_readable(file_path: Path) -> bool:
    """
    Check whether a file can actually be opened for reading.

    Uses an open/close attempt rather than only checking permission
    bits, since permission bits alone don't catch every reason a
    file might be unreadable (e.g. a broken symlink reporting valid
    permissions).

    Args:
        file_path: Path to check.

    Returns:
        True if the file could be opened for reading.
    """
    try:
        with open(file_path, "rb"):
            return True
    except OSError:
        return False

    