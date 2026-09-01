"""
BLC Mark - Metadata Manager

Purpose:
    Manage metadata for datasets that have already been downloaded
    into BLC Mark's data/raw/ directory -- what a file is, where
    it came from, and whether its bytes are intact.

Responsibilities:
    - Inspect local dataset files.
    - Generate SHA-256 checksums.
    - Store immutable metadata records.
    - Summarize dataset collections.
    - Define the future persistence interface.

Scope:
    This module manages metadata only.

    It DOES NOT:
        - download datasets
        - parse biological files
        - normalize expression matrices
        - validate biology
        - perform biomarker discovery

    Those responsibilities belong to other modules.
"""

from __future__ import annotations

import hashlib

from dataclasses import dataclass, replace
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path

from dataset_registry import CancerType, DataType


# ==========================================================
# Version Information
# ==========================================================

METADATA_MANAGER_VERSION = "1.0"
METADATA_SCHEMA_VERSION = "1.0"


# ==========================================================
# Project Paths
# ==========================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_METADATA_DIR = PROJECT_ROOT / "data" / "metadata"


# ==========================================================
# Checksum Configuration
# ==========================================================

# Stream files in 1 MB chunks so even multi-gigabyte TCGA files
# can be processed without loading everything into memory.
CHECKSUM_CHUNK_BYTES = 1024 * 1024


# ==========================================================
# Public API
# ==========================================================

__all__ = [
    "METADATA_MANAGER_VERSION",
    "METADATA_SCHEMA_VERSION",
    "DEFAULT_METADATA_DIR",
    "MetadataStatus",
    "DatasetMetadata",
    "MetadataManager",
]


# ==========================================================
# Metadata Status
# ==========================================================

class MetadataStatus(str, Enum):
    """
    Lifecycle state of a dataset.

    This is different from DownloadStatus.

    DownloadStatus describes ONE download attempt.

    MetadataStatus describes the lifetime of the dataset inside
    BLC Mark.
    """

    REGISTERED = "registered"
    DOWNLOADED = "downloaded"
    VALIDATED = "validated"
    PROCESSED = "processed"
    FAILED = "failed"
    ARCHIVED = "archived"


# ==========================================================
# Dataset Metadata
# ==========================================================

@dataclass(frozen=True)
class DatasetMetadata:
    """
    Immutable provenance record for a dataset.

    Frozen intentionally.

    Metadata records should never be modified in-place because
    that destroys provenance history.

    Any update creates a NEW record using dataclasses.replace().
    """

    dataset_id: str

    cancer_type: CancerType
    data_type: DataType

    local_path: Path

    sha256_checksum: str
    file_size_bytes: int

    source: str | None = None
    xena_host: str | None = None
    download_url: str | None = None

    created_at: datetime | None = None
    modified_at: datetime | None = None

    # Actual download timestamp (independent of filesystem timestamps)
    download_timestamp: datetime | None = None

    status: MetadataStatus = MetadataStatus.DOWNLOADED

    metadata_version: str = METADATA_SCHEMA_VERSION

    # Version of the downloaded dataset itself
    dataset_version: str | None = None

    notes: str | None = None


# ==========================================================
# Metadata Manager
# ==========================================================

# MetadataManager intentionally never parses biological data.
#
# Expression parsing belongs in preprocessing/.
#
# Differential expression belongs in biomarker_engine/.
#
# This separation keeps the architecture modular,
# reproducible and easy to test.


class MetadataManager:
    """
    Produces DatasetMetadata records and summarizes them.

    Currently implemented:

        • generate_checksum()
        • collect_metadata()
        • build_summary()
        • with_status()
        • render_summary_text()

    Future versions will add:

        • JSON persistence
        • CSV export
        • Metadata lookup
        • Duplicate detection
        • Dashboard integration
    """

    def __init__(
        self,
        metadata_directory: Path = DEFAULT_METADATA_DIR,
    ) -> None:
        """
        Store where metadata will eventually be persisted.

        No filesystem operations occur here.
        """

        self.metadata_directory = metadata_directory

    # ==========================================================
    # Core Metadata Collection
    # ==========================================================

    def generate_checksum(self, file_path: Path) -> str:
        """
        Compute the SHA-256 checksum of a file.

        The file is streamed in fixed-size chunks so memory usage
        remains constant regardless of file size.
        """

        if not file_path.exists():
            raise FileNotFoundError(
                f"Cannot checksum missing file: {file_path}"
            )

        if file_path.is_dir():
            raise IsADirectoryError(
                f"Cannot checksum a directory: {file_path}"
            )

        digest = hashlib.sha256()

        with open(file_path, "rb") as handle:
            for chunk in iter(
                lambda: handle.read(CHECKSUM_CHUNK_BYTES),
                b"",
            ):
                digest.update(chunk)

        return digest.hexdigest()

    def collect_metadata(
        self,
        file_path: Path,
        dataset_id: str,
        cancer_type: CancerType,
        data_type: DataType,
        *,
        source: str | None = None,
        xena_host: str | None = None,
        download_url: str | None = None,
        status: MetadataStatus = MetadataStatus.DOWNLOADED,
        notes: str | None = None,
        dataset_version: str | None = None,
    ) -> DatasetMetadata:
        """
        Build a DatasetMetadata record from an existing file.
        """

        if not file_path.exists():
            raise FileNotFoundError(
                f"Cannot collect metadata for missing file: {file_path}"
            )

        stat = file_path.stat()

        checksum = self.generate_checksum(file_path)

        timestamp = datetime.now(timezone.utc)

        return DatasetMetadata(
            dataset_id=dataset_id,
            cancer_type=cancer_type,
            data_type=data_type,
            local_path=file_path,
            sha256_checksum=checksum,
            file_size_bytes=stat.st_size,
            source=source,
            xena_host=xena_host,
            download_url=download_url,
            created_at=datetime.fromtimestamp(
                stat.st_ctime,
                tz=timezone.utc,
            ),
            modified_at=datetime.fromtimestamp(
                stat.st_mtime,
                tz=timezone.utc,
            ),
            download_timestamp=timestamp,
            status=status,
            dataset_version=dataset_version,
            notes=notes,
        )

    # ==========================================================
    # Summary Utilities
    # ==========================================================

    def build_summary(
        self,
        records: list[DatasetMetadata],
    ) -> dict:
        """
        Produce summary statistics for metadata records.
        """

        counts_by_status: dict[MetadataStatus, int] = {}
        counts_by_cancer_type: dict[CancerType, int] = {}

        for record in records:

            counts_by_status[record.status] = (
                counts_by_status.get(record.status, 0) + 1
            )

            counts_by_cancer_type[record.cancer_type] = (
                counts_by_cancer_type.get(record.cancer_type, 0) + 1
            )

        return {
            "manager_version": METADATA_MANAGER_VERSION,
            "metadata_schema_version": METADATA_SCHEMA_VERSION,
            "total_datasets": len(records),
            "total_size_bytes": sum(
                record.file_size_bytes
                for record in records
            ),
            "counts_by_status": counts_by_status,
            "counts_by_cancer_type": counts_by_cancer_type,
            "dataset_ids": [
                record.dataset_id
                for record in records
            ],
        }

    def render_summary_text(
        self,
        summary: dict,
    ) -> str:
        """
        Return a nicely formatted summary.
        """

        lines = [
            "Metadata Summary",
            "-" * 40,
            f"Manager Version : {summary['manager_version']}",
            f"Schema Version  : {summary['metadata_schema_version']}",
            "",
            f"Datasets        : {summary['total_datasets']}",
            f"Total Size      : {summary['total_size_bytes']:,} bytes",
            "",
            "Datasets by Status",
        ]

        for status, count in summary["counts_by_status"].items():
            lines.append(
                f"  {status.value:<12} {count}"
            )

        lines.append("")
        lines.append("Datasets by Cancer")

        for cancer, count in summary["counts_by_cancer_type"].items():
            lines.append(
                f"  {cancer.value:<12} {count}"
            )

        return "\n".join(lines)

    # ==========================================================
    # Convenience Helpers
    # ==========================================================

    def get_file_size_mb(
        self,
        record: DatasetMetadata,
    ) -> float:
        """
        Return file size in megabytes.
        """

        return record.file_size_bytes / (1024 * 1024)

    def metadata_exists(
        self,
        dataset_id: str,
    ) -> bool:
        """
        Placeholder for future persistence support.
        """

        raise NotImplementedError(
            "MetadataManager.metadata_exists() will be "
            "implemented after JSON persistence is added."
        )

    def with_status(
        self,
        record: DatasetMetadata,
        status: MetadataStatus,
    ) -> DatasetMetadata:
        """
        Return a copy of a DatasetMetadata object with an
        updated lifecycle status.
        """

        return replace(
            record,
            status=status,
        )

    # ==========================================================
    # Future Persistence Layer
    # ==========================================================

    def save_metadata(self, record: DatasetMetadata) -> None:
        """
        Persist a DatasetMetadata record to storage.

        Future implementation will most likely serialize metadata as
        JSON files inside the metadata directory.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError(
            "MetadataManager.save_metadata() is not implemented yet. "
            f"It will persist metadata for dataset "
            f"'{record.dataset_id}' under "
            f"{self.metadata_directory}."
        )

    def load_metadata(
        self,
        dataset_id: str,
    ) -> DatasetMetadata:
        """
        Load metadata for a dataset.

        Future implementation will deserialize a previously saved
        DatasetMetadata object.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError(
            "MetadataManager.load_metadata() is not implemented yet. "
            f"It will load metadata for '{dataset_id}'."
        )

    def export_json(
        self,
        records: list[DatasetMetadata],
        destination: Path,
    ) -> None:
        """
        Export metadata records as a JSON file.

        Future implementation will serialize Path, Enum and datetime
        objects into JSON-safe values.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError(
            "MetadataManager.export_json() is not implemented yet. "
            f"It will export {len(records)} record(s) "
            f"to '{destination}'."
        )

    def export_csv(
        self,
        records: list[DatasetMetadata],
        destination: Path,
    ) -> None:
        """
        Export metadata records to CSV.

        Future implementation will flatten DatasetMetadata into rows
        suitable for spreadsheet software.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError(
            "MetadataManager.export_csv() is not implemented yet. "
            f"It will export {len(records)} record(s) "
            f"to '{destination}'."
        )


# ==========================================================
# Future Extensions
#
# Planned after Phase 2:
#
# • JSON persistence
# • Metadata lookup
# • Duplicate detection
# • Metadata validation
# • Registry synchronization
# • Dashboard integration
#
# Candidate layout:
#
# data/
#     metadata/
#         TCGA-BRCA.json
#         TCGA-LUAD.json
#         TCGA-COAD.json
#
# ==========================================================


if __name__ == "__main__":
    import tempfile

    print("=" * 60)
    print("BLC Mark Metadata Manager")
    print("=" * 60)
    print(f"Version          : {METADATA_MANAGER_VERSION}")
    print(f"Schema Version   : {METADATA_SCHEMA_VERSION}")
    print(f"Metadata Folder  : {DEFAULT_METADATA_DIR}")
    print()

    manager = MetadataManager()

    # ------------------------------------------------------
    # Synthetic sanity check
    #
    # No biological data is used here. We create a temporary
    # text file purely to verify checksum generation and
    # metadata collection.
    # ------------------------------------------------------

    with tempfile.TemporaryDirectory() as scratch:

        scratch_file = Path(scratch) / "synthetic_test_file.txt"

        scratch_file.write_text(
            "BLC Mark Metadata Manager Sanity Check\n"
        )

        record = manager.collect_metadata(
            file_path=scratch_file,
            dataset_id="SYNTHETIC.TEST.FILE",
            cancer_type=CancerType.BRCA,
            data_type=DataType.GENE_EXPRESSION_RNASEQ,
            source="Synthetic Test",
            notes="Temporary file used for MetadataManager validation.",
            dataset_version="1.0",
        )

        summary = manager.build_summary([record])

        print("Metadata Record")
        print("-" * 40)
        print(f"Dataset ID       : {record.dataset_id}")
        print(f"Checksum         : {record.sha256_checksum}")
        print(f"File Size        : {record.file_size_bytes} bytes")
        print(f"File Size (MB)   : {manager.get_file_size_mb(record):.4f}")
        print(f"Status           : {record.status.value}")
        print(f"Download Time    : {record.download_timestamp}")
        print()

        print(manager.render_summary_text(summary))
        print()

        print("Sanity Check Results")
        print("-" * 40)
        print("✓ Metadata collection")
        print("✓ SHA-256 checksum generation")
        print("✓ Dataset summary generation")
        print("✓ File size calculation")
        print("✓ MetadataManager initialized successfully")
        print()

        print("Status : READY FOR PHASE 3")
