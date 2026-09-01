"""Result writing for BLC Mark Phase 4 evidence integration."""

import csv
import json
from dataclasses import asdict
from pathlib import Path

from .aggregation import GeneEvidenceProfile


EVIDENCE_OUTPUT_COLUMNS = (
    "gene_id",
    "cancer_cohort",
    "evidence_type",
    "source",
    "source_version",
    "evidence_id",
    "description",
    "retrieved_at",
    "source_url",
)


def write_evidence_profiles(
    profiles: list[GeneEvidenceProfile],
    output_csv: str | Path,
) -> Path:
    """Write Phase 4 evidence records to a deterministic CSV file."""
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    rows: list[dict[str, str | None]] = []

    for profile in profiles:
        for record in profile.evidence_records:
            rows.append(
                {
                    "gene_id": profile.gene_id,
                    "cancer_cohort": profile.cancer_cohort,
                    "evidence_type": record.evidence_type.value,
                    "source": record.source,
                    "source_version": record.source_version,
                    "evidence_id": record.evidence_id,
                    "description": record.description,
                    "retrieved_at": record.retrieved_at,
                    "source_url": record.source_url,
                }
            )

    rows.sort(
        key=lambda row: (
            str(row["gene_id"]),
            str(row["evidence_type"]),
            str(row["source"]),
            str(row["evidence_id"] or ""),
            str(row["description"]),
        )
    )

    with output_csv.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=EVIDENCE_OUTPUT_COLUMNS,
        )

        writer.writeheader()
        writer.writerows(rows)

    return output_csv


def write_json(
    payload: dict,
    output_path: str | Path,
) -> Path:
    """Write deterministic UTF-8 JSON metadata/QC output."""
    output_path = Path(output_path)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            payload,
            handle,
            indent=2,
            sort_keys=True,
        )

        handle.write("\n")

    return output_path