"""Quality control for BLC Mark Phase 4 evidence integration."""

from collections import Counter

from .aggregation import GeneEvidenceProfile


def build_phase4_qc_report(
    profiles: list[GeneEvidenceProfile],
    *,
    candidate_count: int,
    unresolved_identifier_count: int,
) -> dict:
    """Build deterministic Phase 4 QC summary statistics."""
    genes_with_evidence = 0
    genes_without_evidence = 0

    evidence_type_counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()

    total_evidence_records = 0

    for profile in profiles:
        if profile.evidence_records:
            genes_with_evidence += 1
        else:
            genes_without_evidence += 1

        for record in profile.evidence_records:
            total_evidence_records += 1
            evidence_type_counts[
                record.evidence_type.value
            ] += 1
            source_counts[
                record.source
            ] += 1

    if len(profiles) != candidate_count:
        raise ValueError(
            "Profile count does not match candidate_count."
        )

    return {
        "candidate_count": candidate_count,
        "profile_count": len(profiles),
        "genes_with_evidence": genes_with_evidence,
        "genes_without_evidence": genes_without_evidence,
        "unresolved_identifier_count": unresolved_identifier_count,
        "total_evidence_records": total_evidence_records,
        "evidence_type_counts": dict(
            sorted(evidence_type_counts.items())
        ),
        "source_counts": dict(
            sorted(source_counts.items())
        ),
    }