from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import scipy


PHASE6_FILES = [
    "results/phase6/tables/brca_top_25_biomarkers.csv",
    "results/phase6/tables/luad_top_25_biomarkers.csv",
    "results/phase6/tables/coad_top_25_biomarkers.csv",
    "results/phase6/tables/combined_top_25_biomarkers.csv",
    "results/phase6/figures/brca_top_15_biomarkers.png",
    "results/phase6/figures/luad_top_15_biomarkers.png",
    "results/phase6/figures/coad_top_15_biomarkers.png",
    "results/phase6/figures/cross_cohort_top_15_score_profiles.png",
    "results/phase6/report/BLC_Mark_V1_scientific_report.md",
]

GENERATOR_SCRIPTS = [
    "scripts/generate_phase6_tables.py",
    "scripts/generate_phase6_figures.py",
    "scripts/generate_phase6_report.py",
    "scripts/generate_phase6_manifest.py",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def file_record(
    repo_root: Path,
    relative_path: str,
) -> dict:
    path = repo_root / relative_path

    if not path.exists():
        raise FileNotFoundError(
            f"Required reproducibility file is missing: {path}"
        )

    stat = path.stat()

    return {
        "path": relative_path.replace("\\", "/"),
        "size_bytes": stat.st_size,
        "sha256": sha256_file(path),
        "modified_time_utc": datetime.fromtimestamp(
            stat.st_mtime,
            timezone.utc,
        ).isoformat(),
    }


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]

    output_directory = (
        repo_root
        / "results"
        / "phase6"
        / "reproducibility"
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    manifest_path = (
        output_directory
        / "phase6_reproducibility_manifest.json"
    )

    manifest = {
        "project": "BLC Mark",
        "version": "1.0",
        "phase": 6,
        "phase_name": "Scientific Outputs & Reproducibility",
        "generated_at_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "executed_cohorts": [
            "TCGA-BRCA",
            "TCGA-LUAD",
            "TCGA-COAD",
        ],
        "data_scope": "RNA-seq gene expression",
        "environment": {
            "python_version": sys.version.split()[0],
            "platform": platform.platform(),
            "pandas_version": pd.__version__,
            "scipy_version": scipy.__version__,
        },
        "phase6_outputs": [
            file_record(
                repo_root,
                relative_path,
            )
            for relative_path in PHASE6_FILES
        ],
        "generator_scripts": [
            file_record(
                repo_root,
                relative_path,
            )
            for relative_path in GENERATOR_SCRIPTS
        ],
        "reproducibility_notes": [
            (
                "Phase 6 tables are generated from frozen "
                "Phase 5 prioritization outputs."
            ),
            (
                "Phase 6 figures are generated programmatically "
                "from prioritization results."
            ),
            (
                "The scientific report is generated from project "
                "outputs rather than manually reconstructed values."
            ),
            (
                "SHA-256 hashes permit integrity verification of "
                "the recorded release artifacts."
            ),
        ],
    }

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        f"Wrote reproducibility manifest to {manifest_path}"
    )

    print(
        f"Recorded {len(manifest['phase6_outputs'])} Phase 6 outputs."
    )

    print(
        f"Recorded {len(manifest['generator_scripts'])} generator scripts."
    )


if __name__ == "__main__":
    main()