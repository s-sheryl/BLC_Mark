"""Final BLC Mark Phase 4 evidence-integration execution."""

from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path

from src.evidence_integration.analysis import run_phase4_analysis
from src.evidence_integration.candidates import extract_significant_candidates
from src.evidence_integration.cross_cancer import build_cross_cancer_index
from src.evidence_integration.hpa_prognostic import (
    load_hpa_prognostic_data,
)
from src.evidence_integration.ncbi_local import (
    NCBIGeneRecord,
    collect_local_functional_evidence,
    load_human_gene_info,
)
from src.evidence_integration.open_targets_local import (
    OpenTargetsAssociation,
    load_open_targets_associations,
)
from src.evidence_integration.qc import build_phase4_qc_report
from src.evidence_integration.reactome import load_reactome_mapping
from src.evidence_integration.results import (
    write_evidence_profiles,
    write_json,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------
# Phase 3 inputs
# ---------------------------------------------------------------------

PHASE3_FILES = {
    "TCGA-BRCA": (
        PROJECT_ROOT
        / "results"
        / "phase3"
        / "TCGA-BRCA"
        / "brca_v1_de_differential_expression_results.csv"
    ),
    "TCGA-LUAD": (
        PROJECT_ROOT
        / "results"
        / "phase3"
        / "TCGA-LUAD"
        / "luad_v1_de_differential_expression_results.csv"
    ),
    "TCGA-COAD": (
        PROJECT_ROOT
        / "results"
        / "phase3"
        / "TCGA-COAD"
        / "coad_v1_de_differential_expression_results.csv"
    ),
}


# ---------------------------------------------------------------------
# Frozen/local Phase 4 evidence sources
# ---------------------------------------------------------------------

NCBI_GENE_INFO_FILE = (
    PROJECT_ROOT
    / "data"
    / "external"
    / "ncbi"
    / "Homo_sapiens.gene_info.gz"
)

REACTOME_FILE = (
    PROJECT_ROOT
    / "data"
    / "external"
    / "reactome"
    / "NCBI2Reactome_V97.txt"
)

HPA_FILE = (
    PROJECT_ROOT
    / "data"
    / "external"
    / "hpa"
    / "cancer_prognostic_data_HPA25.1.tsv.zip"
)

OPEN_TARGETS_FILE = (
    PROJECT_ROOT
    / "data"
    / "external"
    / "open_targets"
    / "blc_mark_cancer_associations_26.06.parquet"
)


# ---------------------------------------------------------------------
# Source versions
# ---------------------------------------------------------------------

NCBI_SNAPSHOT_VERSION = "2026-08-26-snapshot"
REACTOME_VERSION = "97"
HPA_VERSION = "25.1"
OPEN_TARGETS_VERSION = "26.06"


# ---------------------------------------------------------------------
# Verified input hashes already recorded during acquisition
# ---------------------------------------------------------------------

NCBI_SHA256 = (
    "6D7EDBFBE9BA56523895A96C165005B83470DA002087C590142DCBFF9290824B"
)

REACTOME_SHA256 = (
    "44AEC2F3D6C0D6F03268FBDC66F95A14673E86DA31682B7EDBD71C6C6A1ECC8D"
)

HPA_SHA256 = (
    "86624E44FB366DF444D07721996A93394BF0490230FF4FBF75934E6362927EE8"
)


# ---------------------------------------------------------------------
# Loaded local indexes
# ---------------------------------------------------------------------

NCBI_INDEX: dict[str, NCBIGeneRecord] = {}


def calculate_sha256(path: Path) -> str:
    """Calculate SHA-256 for a local evidence-source file."""
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest().upper()


# ---------------------------------------------------------------------
# Local identifier and functional evidence access
# ---------------------------------------------------------------------

def local_ncbi_gene_id(
    gene_symbol: str,
) -> str | None:
    """Resolve a symbol to NCBI Gene ID from the local snapshot."""
    record = NCBI_INDEX.get(gene_symbol)

    if record is None:
        return None

    return record.gene_id


def local_ensembl_id(
    gene_symbol: str,
) -> str | None:
    """Resolve a symbol to Ensembl Gene ID from the local snapshot."""
    record = NCBI_INDEX.get(gene_symbol)

    if record is None:
        return None

    return record.ensembl_id


def local_functional_collector(
    gene_symbol: str,
    retrieved_at: str,
):
    """Collect functional evidence from local NCBI Gene metadata."""
    return collect_local_functional_evidence(
        gene_symbol,
        ncbi_index=NCBI_INDEX,
        retrieved_at=retrieved_at,
        source_version=NCBI_SNAPSHOT_VERSION,
    )


# ---------------------------------------------------------------------
# Cross-cancer candidate recurrence
# ---------------------------------------------------------------------

def build_real_cross_cancer_index() -> dict[str, tuple[str, ...]]:
    """Build cross-cancer evidence from real Phase 3 candidates."""
    cohort_candidates: dict[str, set[str]] = {}

    for cohort, path in PHASE3_FILES.items():
        candidates = extract_significant_candidates(
            path
        )

        cohort_candidates[cohort] = set(
            candidates["gene_id"]
        )

    return build_cross_cancer_index(
        cohort_candidates
    )


# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------

def validate_inputs() -> None:
    """Fail explicitly if any required Phase 4 input is absent."""
    for cohort, path in PHASE3_FILES.items():
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing Phase 3 file for {cohort}: {path}"
            )

    if not NCBI_GENE_INFO_FILE.is_file():
        raise FileNotFoundError(
            "Missing NCBI Homo sapiens gene_info file: "
            f"{NCBI_GENE_INFO_FILE}"
        )

    if not REACTOME_FILE.is_file():
        raise FileNotFoundError(
            f"Missing Reactome V97 file: {REACTOME_FILE}"
        )

    if not HPA_FILE.is_file():
        raise FileNotFoundError(
            f"Missing HPA 25.1 file: {HPA_FILE}"
        )

    if not OPEN_TARGETS_FILE.is_file():
        raise FileNotFoundError(
            "Missing filtered Open Targets 26.06 file: "
            f"{OPEN_TARGETS_FILE}"
        )


# ---------------------------------------------------------------------
# Cohort execution
# ---------------------------------------------------------------------

def run_cohort(
    cohort: str,
    phase3_path: Path,
    *,
    cross_cancer_index: dict[str, tuple[str, ...]],
    reactome_mapping: dict[str, list[dict[str, str]]],
    hpa_index: dict[tuple[str, str], dict[str, str]],
    open_targets_index: dict[
        tuple[str, str],
        OpenTargetsAssociation,
    ],
    retrieved_at: str,
    open_targets_sha256: str,
) -> None:
    """Execute and write Phase 4 outputs for one cohort."""
    print()
    print("=" * 72)
    print(f"Running Phase 4: {cohort}")
    print("=" * 72)

    result = run_phase4_analysis(
        phase3_path,
        cohort,
        cross_cancer_index=cross_cancer_index,
        reactome_mapping=reactome_mapping,
        hpa_index=hpa_index,
        open_targets_index=open_targets_index,
        retrieved_at=retrieved_at,
        reactome_version=REACTOME_VERSION,
        hpa_version=HPA_VERSION,
        open_targets_version=OPEN_TARGETS_VERSION,
        resolve_ncbi_gene_id=local_ncbi_gene_id,
        resolve_ensembl_id=local_ensembl_id,
        collect_functional=local_functional_collector,
    )

    output_dir = (
        PROJECT_ROOT
        / "results"
        / "phase4"
        / cohort
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    evidence_path = output_dir / "evidence.csv"
    qc_path = output_dir / "qc_report.json"
    metadata_path = output_dir / "metadata.json"

    profiles = list(
        result.profiles
    )

    write_evidence_profiles(
        profiles,
        evidence_path,
    )

    qc = build_phase4_qc_report(
        profiles,
        candidate_count=result.candidate_count,
        unresolved_identifier_count=(
            result.unresolved_identifier_count
        ),
    )

    metadata = {
        "phase": 4,
        "phase_name": "Evidence Integration",
        "blc_mark_version": "1.0",
        "cancer_cohort": cohort,
        "retrieved_at": retrieved_at,
        "candidate_count": result.candidate_count,
        "unresolved_identifier_count": (
            result.unresolved_identifier_count
        ),
        "ranking_performed": False,
        "sources": {
            "ncbi_gene": {
                "name": "NCBI Gene",
                "access_mode": "local bulk snapshot",
                "version": NCBI_SNAPSHOT_VERSION,
                "file": "Homo_sapiens.gene_info.gz",
                "sha256": NCBI_SHA256,
                "uses": [
                    "NCBI Gene ID resolution",
                    "Ensembl Gene ID resolution",
                    "functional description",
                ],
            },
            "reactome": {
                "name": "Reactome",
                "access_mode": "local bulk snapshot",
                "version": REACTOME_VERSION,
                "file": "NCBI2Reactome_V97.txt",
                "sha256": REACTOME_SHA256,
                "uses": [
                    "pathway evidence",
                ],
            },
            "open_targets": {
                "name": "Open Targets Platform",
                "access_mode": (
                    "local filtered bulk snapshot"
                ),
                "version": OPEN_TARGETS_VERSION,
                "file": (
                    "blc_mark_cancer_associations_26.06.parquet"
                ),
                "sha256": open_targets_sha256,
                "disease_mapping": {
                    "TCGA-BRCA": {
                        "id": "MONDO_0007254",
                        "name": "breast cancer",
                    },
                    "TCGA-LUAD": {
                        "id": "MONDO_0005061",
                        "name": "lung adenocarcinoma",
                    },
                    "TCGA-COAD": {
                        "id": "MONDO_0002271",
                        "name": "colon adenocarcinoma",
                    },
                },
                "uses": [
                    "direct cancer-association evidence",
                ],
            },
            "human_protein_atlas": {
                "name": "Human Protein Atlas",
                "access_mode": "local bulk snapshot",
                "version": HPA_VERSION,
                "file": (
                    "cancer_prognostic_data_HPA25.1.tsv.zip"
                ),
                "sha256": HPA_SHA256,
                "uses": [
                    "clinical/prognostic evidence",
                ],
            },
            "cross_cancer": {
                "name": (
                    "BLC Mark Phase 3 cross-cohort evidence"
                ),
                "access_mode": (
                    "internal Phase 3 outputs"
                ),
                "version": "1.0",
                "uses": [
                    "cross-cancer candidate recurrence",
                ],
            },
        },
    }

    write_json(
        qc,
        qc_path,
    )

    write_json(
        metadata,
        metadata_path,
    )

    print(
        "Candidates:",
        result.candidate_count,
    )

    print(
        "Unresolved identifiers:",
        result.unresolved_identifier_count,
    )

    print(
        "Genes with evidence:",
        qc["genes_with_evidence"],
    )

    print(
        "Genes without evidence:",
        qc["genes_without_evidence"],
    )

    print(
        "Evidence records:",
        qc["total_evidence_records"],
    )

    print(
        "Evidence type counts:",
        qc["evidence_type_counts"],
    )

    print(
        "Source counts:",
        qc["source_counts"],
    )

    print(
        "Evidence output:",
        evidence_path,
    )

    print(
        "QC output:",
        qc_path,
    )

    print(
        "Metadata output:",
        metadata_path,
    )


# ---------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------

def main() -> None:
    """Execute BLC Mark Phase 4 for all V1 cohorts."""
    global NCBI_INDEX

    validate_inputs()

    print(
        "Loading NCBI Homo sapiens gene_info..."
    )

    NCBI_INDEX = load_human_gene_info(
        NCBI_GENE_INFO_FILE
    )

    print(
        "NCBI genes loaded:",
        len(NCBI_INDEX),
    )

    print(
        "Loading Reactome V97..."
    )

    reactome_mapping = load_reactome_mapping(
        REACTOME_FILE
    )

    print(
        "Reactome mapped NCBI Gene IDs:",
        len(reactome_mapping),
    )

    print(
        "Loading Human Protein Atlas 25.1..."
    )

    hpa_index = load_hpa_prognostic_data(
        HPA_FILE
    )

    print(
        "HPA prognostic records loaded:",
        len(hpa_index),
    )

    print(
        "Loading Open Targets 26.06 cancer associations..."
    )

    open_targets_index = (
        load_open_targets_associations(
            OPEN_TARGETS_FILE
        )
    )

    print(
        "Open Targets associations loaded:",
        len(open_targets_index),
    )

    print(
        "Calculating Open Targets snapshot SHA-256..."
    )

    open_targets_sha256 = calculate_sha256(
        OPEN_TARGETS_FILE
    )

    print(
        "Open Targets SHA-256:",
        open_targets_sha256,
    )

    print(
        "Building cross-cancer index..."
    )

    cross_cancer_index = (
        build_real_cross_cancer_index()
    )

    print(
        "Unique Phase 3 candidate genes:",
        len(cross_cancer_index),
    )

    retrieved_at = (
        datetime.now(
            timezone.utc
        ).isoformat()
    )

    for cohort, phase3_path in PHASE3_FILES.items():
        run_cohort(
            cohort,
            phase3_path,
            cross_cancer_index=cross_cancer_index,
            reactome_mapping=reactome_mapping,
            hpa_index=hpa_index,
            open_targets_index=open_targets_index,
            retrieved_at=retrieved_at,
            open_targets_sha256=(
                open_targets_sha256
            ),
        )

    print()
    print("=" * 72)
    print(
        "BLC Mark Phase 4 execution completed."
    )
    print("=" * 72)


if __name__ == "__main__":
    main()