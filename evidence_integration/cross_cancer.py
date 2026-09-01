"""Cross-cancer evidence integration for BLC Mark Phase 4."""

from collections.abc import Mapping

from .models import EvidenceRecord, EvidenceType


def build_cross_cancer_index(
    cohort_candidates: Mapping[str, set[str]],
) -> dict[str, tuple[str, ...]]:
    """
    Build a deterministic gene-to-cohort membership index.

    Parameters
    ----------
    cohort_candidates:
        Mapping from cohort name to significant candidate gene symbols.

    Returns
    -------
    dict[str, tuple[str, ...]]
        Gene symbol mapped to sorted cohort names in which the gene
        was a significant Phase 3 candidate.
    """
    if not isinstance(cohort_candidates, Mapping):
        raise TypeError("cohort_candidates must be a mapping.")

    index: dict[str, set[str]] = {}

    for cohort, genes in cohort_candidates.items():
        if not isinstance(cohort, str) or not cohort.strip():
            raise ValueError(
                "Cohort names must be non-empty strings."
            )

        if not isinstance(genes, set):
            raise TypeError(
                "Each cohort candidate collection must be a set."
            )

        cohort_name = cohort.strip()

        for gene in genes:
            if not isinstance(gene, str) or not gene.strip():
                raise ValueError(
                    "Candidate gene identifiers must be non-empty strings."
                )

            symbol = gene.strip()

            index.setdefault(symbol, set()).add(cohort_name)

    return {
        gene: tuple(sorted(cohorts))
        for gene, cohorts in sorted(index.items())
    }


def collect_cross_cancer_evidence(
    gene_symbol: str,
    cancer_cohort: str,
    *,
    cross_cancer_index: dict[str, tuple[str, ...]],
    retrieved_at: str,
) -> list[EvidenceRecord]:
    """
    Create cross-cancer evidence for one candidate gene.

    Evidence is created only when the gene is significant in more than
    one cancer cohort. No ranking or weighting is performed.
    """
    if not isinstance(gene_symbol, str) or not gene_symbol.strip():
        raise ValueError(
            "gene_symbol must be a non-empty string."
        )

    if not isinstance(cancer_cohort, str) or not cancer_cohort.strip():
        raise ValueError(
            "cancer_cohort must be a non-empty string."
        )

    symbol = gene_symbol.strip()
    cohort = cancer_cohort.strip()

    cohorts = cross_cancer_index.get(symbol, ())

    if cohort not in cohorts:
        return []

    if len(cohorts) <= 1:
        return []

    cohort_text = ", ".join(cohorts)

    record = EvidenceRecord(
        gene_id=symbol,
        evidence_type=EvidenceType.CROSS_CANCER,
        source="BLC Mark Phase 3",
        source_version="1.0",
        evidence_id=f"CROSS_CANCER:{symbol}",
        description=(
            f"Significant differential-expression candidate in "
            f"{len(cohorts)} cohorts: {cohort_text}"
        ),
        cancer_cohort=cohort,
        retrieved_at=retrieved_at,
        source_url=None,
    )

    return [record]