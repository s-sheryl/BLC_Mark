"""Exceptions for BLC Mark Phase 5 biomarker prioritization."""


class BiomarkerPrioritizationError(Exception):
    """Base exception for Phase 5 prioritization failures."""


class PrioritizationValidationError(BiomarkerPrioritizationError):
    """Raised when Phase 5 input validation fails."""