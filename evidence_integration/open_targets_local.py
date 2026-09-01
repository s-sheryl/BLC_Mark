"""Local Open Targets cancer-association evidence for BLC Mark Phase 4."""

from dataclasses import dataclass
from pathlib import Path

import pyarrow.parquet as pq

from .models import EvidenceRecord, EvidenceType


BLC_MARK_DISEASE_IDS = {
    "TCGA-BRCA": "MONDO_0007254",
    "TCGA-LUAD": "MONDO_0005061",
    "TCGA-COAD": "MONDO_0002271",
}


BLC_MARK_DISEASE_NAMES = {
    "MONDO_0007254": "breast cancer",
    "MONDO_0005061": "lung adenocarcinoma",
    "MONDO_0002271": "colon adenocarcinoma",
}


@dataclass(frozen=True)
class OpenTargetsAssociation:
    """One frozen Open Targets target-disease association."""

    disease_id: str
    target_id: str
    association_score: float
    evidence_count: int


def load_open_targets_associations(
    path: str | Path,
) -> dict[
    tuple[str, str],
    OpenTargetsAssociation,
]:
    """
    Load the filtered Open Targets 26.06 association snapshot.

    Records are indexed by (disease_id, target_id).
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"Open Targets association file does not exist: {path}"
        )

    if not path.is_file():
        raise ValueError(
            f"Open Targets association path is not a file: {path}"
        )

    table = pq.read_table(
        path,
        columns=[
            "diseaseId",
            "targetId",
            "associationScore",
            "evidenceCount",
        ],
    )

    index: dict[
        tuple[str, str],
        OpenTargetsAssociation,
    ] = {}

    for row in table.to_pylist():
        disease_id = str(
            row["diseaseId"]
        ).strip()

        target_id = str(
            row["targetId"]
        ).strip()

        score = row["associationScore"]
        evidence_count = row["evidenceCount"]

        if not disease_id:
            raise ValueError(
                "Open Targets record has blank diseaseId."
            )

        if not target_id:
            raise ValueError(
                "Open Targets record has blank targetId."
            )

        if disease_id not in BLC_MARK_DISEASE_NAMES:
            raise ValueError(
                "Unexpected disease ID in filtered Open Targets "
                f"dataset: {disease_id}"
            )

        if score is None:
            raise ValueError(
                "Open Targets record has missing associationScore."
            )

        if evidence_count is None:
            raise ValueError(
                "Open Targets record has missing evidenceCount."
            )

        association = OpenTargetsAssociation(
            disease_id=disease_id,
            target_id=target_id,
            association_score=float(score),
            evidence_count=int(evidence_count),
        )

        key = (
            disease_id,
            target_id,
        )

        if key in index:
            raise ValueError(
                "Duplicate Open Targets association for "
                f"{disease_id} / {target_id}."
            )

        index[key] = association

    return index


def collect_local_cancer_association_evidence(
    gene_symbol: str,
    ensembl_id: str,
    cancer_cohort: str,
    *,
    association_index: dict[
        tuple[str, str],
        OpenTargetsAssociation,
    ],
    retrieved_at: str,
    source_version: str,
) -> list[EvidenceRecord]:
    """
    Create cohort-specific cancer-association evidence from the local
    Open Targets snapshot.

    No ranking or weighting is performed.
    """
    if not isinstance(gene_symbol, str) or not gene_symbol.strip():
        raise ValueError(
            "gene_symbol must be a non-empty string."
        )

    if not isinstance(ensembl_id, str) or not ensembl_id.strip():
        raise ValueError(
            "ensembl_id must be a non-empty string."
        )

    if cancer_cohort not in BLC_MARK_DISEASE_IDS:
        raise ValueError(
            f"Unsupported BLC Mark cancer cohort: {cancer_cohort}"
        )

    if not isinstance(source_version, str) or not source_version.strip():
        raise ValueError(
            "source_version must be a non-empty string."
        )

    disease_id = BLC_MARK_DISEASE_IDS[
        cancer_cohort
    ]

    association = association_index.get(
        (
            disease_id,
            ensembl_id.strip(),
        )
    )

    if association is None:
        return []

    disease_name = BLC_MARK_DISEASE_NAMES[
        disease_id
    ]

    record = EvidenceRecord(
        gene_id=gene_symbol.strip(),
        evidence_type=EvidenceType.CANCER_ASSOCIATION,
        source="Open Targets Platform",
        source_version=source_version.strip(),
        evidence_id=disease_id,
        description=(
            f"{disease_name}; "
            f"Open Targets association score="
            f"{association.association_score:.6f}; "
            f"evidence count={association.evidence_count}"
        ),
        cancer_cohort=cancer_cohort,
        retrieved_at=retrieved_at,
        source_url=(
            "https://platform.opentargets.org/"
            f"disease/{disease_id}"
        ),
    )

    return [record]