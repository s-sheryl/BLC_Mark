"""End-to-end orchestration for BLC Mark Phase 4 evidence integration."""

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from .aggregation import (
    GeneEvidenceProfile,
    aggregate_gene_evidence,
)
from .candidates import extract_significant_candidates
from .cross_cancer import collect_cross_cancer_evidence
from .hpa_prognostic import collect_clinical_evidence
from .identifiers import normalize_gene_symbol
from .models import EvidenceRecord
from .open_targets_local import (
    OpenTargetsAssociation,
    collect_local_cancer_association_evidence,
)
from .reactome import collect_pathway_evidence


@dataclass(frozen=True)
class Phase4AnalysisResult:
    """Result of one cohort-level Phase 4 evidence-integration run."""

    cancer_cohort: str
    candidate_count: int
    unresolved_identifier_count: int
    profiles: tuple[GeneEvidenceProfile, ...]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.cancer_cohort, str)
            or not self.cancer_cohort.strip()
        ):
            raise ValueError(
                "cancer_cohort must be a non-empty string."
            )

        if self.candidate_count < 0:
            raise ValueError(
                "candidate_count cannot be negative."
            )

        if self.unresolved_identifier_count < 0:
            raise ValueError(
                "unresolved_identifier_count cannot be negative."
            )

        if len(self.profiles) != self.candidate_count:
            raise ValueError(
                "Profile count must equal candidate_count."
            )


def run_phase4_analysis(
    phase3_results_path: str | Path,
    cancer_cohort: str,
    *,
    cross_cancer_index: dict[str, tuple[str, ...]],
    reactome_mapping: dict[str, list[dict[str, str]]],
    hpa_index: dict[tuple[str, str], dict[str, str]],
    open_targets_index: dict[
        tuple[str, str],
        OpenTargetsAssociation,
    ],
    retrieved_at: str,
    reactome_version: str,
    hpa_version: str,
    open_targets_version: str,
    resolve_ncbi_gene_id: Callable[[str], str | None],
    resolve_ensembl_id: Callable[[str], str | None],
    collect_functional: Callable[
        [str, str],
        EvidenceRecord | None,
    ],
) -> Phase4AnalysisResult:
    """
    Run Phase 4 evidence integration for one cancer cohort.

    Evidence is aggregated only. No scoring, weighting, ranking,
    or biomarker prioritization is performed in Phase 4.
    """
    if (
        not isinstance(cancer_cohort, str)
        or not cancer_cohort.strip()
    ):
        raise ValueError(
            "cancer_cohort must be a non-empty string."
        )

    if not isinstance(open_targets_index, dict):
        raise TypeError(
            "open_targets_index must be a dictionary."
        )

    candidates = extract_significant_candidates(
        phase3_results_path
    )

    profiles: list[GeneEvidenceProfile] = []
    unresolved_identifier_count = 0

    for gene_id in candidates["gene_id"]:
        identifier = normalize_gene_symbol(
            str(gene_id)
        )

        records: list[EvidenceRecord] = []

        if not identifier.resolvable:
            unresolved_identifier_count += 1

            profiles.append(
                aggregate_gene_evidence(
                    identifier.original_id,
                    cancer_cohort,
                    records,
                )
            )
            continue

        symbol = identifier.normalized_symbol

        if symbol is None:
            raise RuntimeError(
                "Resolvable identifier unexpectedly lacks "
                "normalized_symbol."
            )

        records.extend(
            collect_cross_cancer_evidence(
                symbol,
                cancer_cohort,
                cross_cancer_index=cross_cancer_index,
                retrieved_at=retrieved_at,
            )
        )

        ncbi_gene_id = resolve_ncbi_gene_id(
            symbol
        )

        if ncbi_gene_id is not None:
            functional_record = collect_functional(
                symbol,
                retrieved_at,
            )

            if functional_record is not None:
                records.append(
                    functional_record
                )

            records.extend(
                collect_pathway_evidence(
                    symbol,
                    ncbi_gene_id,
                    reactome_mapping=reactome_mapping,
                    retrieved_at=retrieved_at,
                    source_version=reactome_version,
                )
            )

        ensembl_id = resolve_ensembl_id(
            symbol
        )

        if ensembl_id is not None:
            records.extend(
                collect_local_cancer_association_evidence(
                    symbol,
                    ensembl_id,
                    cancer_cohort,
                    association_index=open_targets_index,
                    retrieved_at=retrieved_at,
                    source_version=open_targets_version,
                )
            )

            records.extend(
                collect_clinical_evidence(
                    symbol,
                    ensembl_id,
                    cancer_cohort,
                    hpa_index=hpa_index,
                    retrieved_at=retrieved_at,
                    source_version=hpa_version,
                )
            )

        profiles.append(
            aggregate_gene_evidence(
                symbol,
                cancer_cohort,
                records,
            )
        )

    profiles.sort(
        key=lambda profile: profile.gene_id
    )

    return Phase4AnalysisResult(
        cancer_cohort=cancer_cohort.strip(),
        candidate_count=len(candidates),
        unresolved_identifier_count=(
            unresolved_identifier_count
        ),
        profiles=tuple(profiles),
    )