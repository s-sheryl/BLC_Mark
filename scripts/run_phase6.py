from __future__ import annotations

import subprocess
import sys
from pathlib import Path


PHASE6_STEPS = [
    ("Generate final biomarker tables", "generate_phase6_tables.py"),
    ("Generate scientific figures", "generate_phase6_figures.py"),
    ("Generate scientific report", "generate_phase6_report.py"),
    ("Generate reproducibility manifest", "generate_phase6_manifest.py"),
]


def run_step(
    scripts_directory: Path,
    description: str,
    script_name: str,
) -> None:
    script_path = scripts_directory / script_name

    if not script_path.is_file():
        raise FileNotFoundError(
            f"Required Phase 6 script is missing: {script_path}"
        )

    print()
    print("=" * 72)
    print(description)
    print(f"Script: {script_name}")
    print("=" * 72)

    completed = subprocess.run(
        [sys.executable, str(script_path)],
        cwd=scripts_directory.parent,
        check=False,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            f"Phase 6 execution failed while running "
            f"{script_name} with exit code "
            f"{completed.returncode}."
        )


def main() -> None:
    scripts_directory = Path(__file__).resolve().parent

    print("=" * 72)
    print("BLC Mark V1 — Phase 6 Reproducible Execution")
    print("=" * 72)

    print(
        "This runner regenerates Phase 6 scientific outputs "
        "from existing frozen upstream results."
    )

    print(
        "It does not rerun Phases 1–5 or retrieve new "
        "external biological evidence."
    )

    for description, script_name in PHASE6_STEPS:
        run_step(
            scripts_directory,
            description,
            script_name,
        )

    print()
    print("=" * 72)
    print("PHASE 6 EXECUTION COMPLETE")
    print("=" * 72)

    print(
        "Generated tables, figures, scientific report, "
        "and reproducibility manifest successfully."
    )


if __name__ == "__main__":
    main()