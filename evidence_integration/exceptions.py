"""Exceptions for the BLC Mark evidence-integration package."""


class EvidenceIntegrationError(Exception):
    """Base exception for Phase 4 evidence-integration failures."""


class EvidenceInputError(EvidenceIntegrationError):
    """Raised when Phase 3 input data is invalid or incompatible."""