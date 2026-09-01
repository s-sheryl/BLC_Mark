"""Core data models for BLC Mark Phase 4 evidence integration."""

from dataclasses import dataclass
from enum import Enum


class EvidenceType(str, Enum):
    """Supported Phase 4 evidence categories."""

    FUNCTIONAL = "functional"
    PATHWAY = "pathway"
    CANCER_ASSOCIATION = "cancer_association"
    CLINICAL = "clinical"
    CROSS_CANCER = "cross_cancer"


@dataclass(frozen=True)
class EvidenceRecord:
    """
    A single traceable piece of evidence associated with a candidate gene.

    Parameters
    ----------
    gene_id:
        Candidate gene symbol inherited from Phase 3.

    evidence_type:
        Category of evidence represented by this record.

    source:
        Name of the external or internal evidence source.

    source_version:
        Version or release identifier when available.

    evidence_id:
        Source-specific stable identifier when available.

    description:
        Human-readable description of the evidence.

    cancer_cohort:
        TCGA cohort associated with the evidence when applicable.

    retrieved_at:
        Retrieval timestamp recorded by the evidence collector.

    source_url:
        Source URL when available.
    """

    gene_id: str
    evidence_type: EvidenceType
    source: str
    description: str
    source_version: str | None = None
    evidence_id: str | None = None
    cancer_cohort: str | None = None
    retrieved_at: str | None = None
    source_url: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.gene_id, str) or not self.gene_id.strip():
            raise ValueError("gene_id must be a non-empty string.")

        if not isinstance(self.evidence_type, EvidenceType):
            raise TypeError("evidence_type must be an EvidenceType.")

        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("source must be a non-empty string.")

        if not isinstance(self.description, str) or not self.description.strip():
            raise ValueError("description must be a non-empty string.")

        object.__setattr__(self, "gene_id", self.gene_id.strip())
        object.__setattr__(self, "source", self.source.strip())
        object.__setattr__(self, "description", self.description.strip())