"""Evidence aggregation for BLC Mark Phase 4."""

from dataclasses import dataclass

from .models import EvidenceRecord


@dataclass(frozen=True)
class GeneEvidenceProfile:
    """
    Aggregated Phase 4 evidence for one candidate gene.

    This object groups evidence records without assigning scores,
    weights, ranks, or priority labels.
    """

    gene_id: str
    cancer_cohort: str
    evidence_records: tuple[EvidenceRecord, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.gene_id, str) or not self.gene_id.strip():
            raise ValueError("gene_id must be a non-empty string.")

        if (
            not isinstance(self.cancer_cohort, str)
            or not self.cancer_cohort.strip()
        ):
            raise ValueError(
                "cancer_cohort must be a non-empty string."
            )

        if not isinstance(self.evidence_records, tuple):
            raise TypeError(
                "evidence_records must be a tuple."
            )

        for record in self.evidence_records:
            if not isinstance(record, EvidenceRecord):
                raise TypeError(
                    "Every evidence record must be an EvidenceRecord."
                )

            if record.gene_id != self.gene_id.strip():
                raise ValueError(
                    "Evidence record gene_id does not match profile gene_id."
                )

        object.__setattr__(
            self,
            "gene_id",
            self.gene_id.strip(),
        )

        object.__setattr__(
            self,
            "cancer_cohort",
            self.cancer_cohort.strip(),
        )


def aggregate_gene_evidence(
    gene_id: str,
    cancer_cohort: str,
    records: list[EvidenceRecord],
) -> GeneEvidenceProfile:
    """
    Aggregate evidence records for one candidate gene.

    Duplicate records are removed using source, evidence type,
    evidence identifier, and description. Output ordering is
    deterministic.

    No evidence scoring or ranking is performed.
    """
    if not isinstance(records, list):
        raise TypeError("records must be a list.")

    gene = gene_id.strip() if isinstance(gene_id, str) else gene_id

    for record in records:
        if not isinstance(record, EvidenceRecord):
            raise TypeError(
                "Every item in records must be an EvidenceRecord."
            )

        if record.gene_id != gene:
            raise ValueError(
                "Cannot aggregate evidence belonging to another gene."
            )

    unique: dict[
        tuple[str, str, str, str],
        EvidenceRecord,
    ] = {}

    for record in records:
        key = (
            record.evidence_type.value,
            record.source,
            record.evidence_id or "",
            record.description,
        )

        unique[key] = record

    ordered = tuple(
        sorted(
            unique.values(),
            key=lambda record: (
                record.evidence_type.value,
                record.source,
                record.evidence_id or "",
                record.description,
            ),
        )
    )

    return GeneEvidenceProfile(
        gene_id=gene_id,
        cancer_cohort=cancer_cohort,
        evidence_records=ordered,
    )