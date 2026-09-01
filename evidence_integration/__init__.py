"""BLC Mark Phase 4 evidence-integration package."""

from .aggregation import (
    GeneEvidenceProfile,
    aggregate_gene_evidence,
)
from .candidates import extract_significant_candidates
from .exceptions import (
    EvidenceInputError,
    EvidenceIntegrationError,
)
from .identifiers import (
    GeneIdentifier,
    normalize_gene_symbol,
)
from .models import (
    EvidenceRecord,
    EvidenceType,
)
from .open_targets import (
    collect_cancer_association_evidence,
    fetch_target_disease_associations,
)
from .open_targets_identifiers import (
    resolve_gene_symbol_to_ensembl,
)
from .provenance import EvidenceSourceMetadata
from .reactome import (
    collect_pathway_evidence,
    load_reactome_mapping,
)
from .validation import (
    REQUIRED_PHASE3_COLUMNS,
    validate_phase3_results,
)


__all__ = [
    "EvidenceInputError",
    "EvidenceIntegrationError",
    "EvidenceRecord",
    "EvidenceSourceMetadata",
    "EvidenceType",
    "GeneEvidenceProfile",
    "GeneIdentifier",
    "REQUIRED_PHASE3_COLUMNS",
    "aggregate_gene_evidence",
    "collect_cancer_association_evidence",
    "collect_pathway_evidence",
    "extract_significant_candidates",
    "fetch_target_disease_associations",
    "load_reactome_mapping",
    "normalize_gene_symbol",
    "resolve_gene_symbol_to_ensembl",
    "validate_phase3_results",
]
