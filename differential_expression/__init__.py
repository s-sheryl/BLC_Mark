"""
BLC Mark Phase 3 -- Differential Expression Analysis.

Public API:
    - run_analysis(configuration) -> AnalysisResult
        The complete orchestrated pipeline (analysis.py).
    - build_configuration(...) -> DEAnalysisConfiguration
        Build and validate a DE analysis configuration
        (configuration.py).

Every other name in this package (validation, comparison, filtering,
methods, multiple_testing, qc, reproducibility, results submodules)
is available via its own submodule import
(e.g. `from src.differential_expression import methods`) for callers
that need finer-grained access, but is not re-exported here, per the
project's "keep the package interface minimal and explicit" rule.

Models and exceptions are intentionally not re-exported wholesale
either: import them from
`src.differential_expression.models` and
`src.differential_expression.exceptions` directly, so it is always
clear which submodule a given type or exception comes from.
"""

from src.differential_expression.analysis import run_analysis
from src.differential_expression.configuration import build_configuration

DIFFERENTIAL_EXPRESSION_PACKAGE_VERSION = "1.0"

__all__ = [
    "DIFFERENTIAL_EXPRESSION_PACKAGE_VERSION",
    "run_analysis",
    "build_configuration",
]