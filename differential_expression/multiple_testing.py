"""
Purpose:
    Apply Benjamini-Hochberg false-discovery-rate correction to the
    raw p-values produced by methods.py, using an established
    scientific library rather than a hand-rolled implementation.

Responsibilities:
    - Accept the raw, per-gene statistical results produced by
      methods.py (adjusted_p_value and significant still None at that
      point).
    - Apply Benjamini-Hochberg correction (scipy.stats.
      false_discovery_control) across every gene that has a valid raw
      p-value.
    - Preserve gene identity and original ordering.
    - Exclude genes with a missing raw p-value from the correction
      itself, while keeping them in the output with adjusted_p_value
      left as None (specification Section 6.4).
    - Annotate each gene's significance against the configured
      significance_threshold.
    - Fail explicitly (MultipleTestingError) if no gene has a valid
      raw p-value to correct, or if an unsupported correction method
      is configured.

Scope:
    This module performs no statistical testing itself -- it only
    corrects p-values already produced by methods.py. It does not
    decide which correction method to use; models.MultipleTestingMethod
    currently defines only Benjamini-Hochberg, per specification
    Section 6.1, and this module rejects anything else rather than
    silently substituting a different correction.
"""

import numpy as np
from scipy.stats import false_discovery_control

from src.differential_expression.exceptions import MultipleTestingError
from src.differential_expression.models import GeneResult, MultipleTestingMethod

MULTIPLE_TESTING_VERSION = "1.0"

__all__ = [
    "MULTIPLE_TESTING_VERSION",
    "apply_multiple_testing_correction",
]


def apply_multiple_testing_correction(
    gene_results: tuple[GeneResult, ...],
    multiple_testing_method: MultipleTestingMethod,
    significance_threshold: float,
) -> tuple[GeneResult, ...]:
    """Apply the configured multiple-testing correction.

    Args:
        gene_results: Per-gene results from methods.py, with
            adjusted_p_value and significant both None.
        multiple_testing_method: The configured correction method.
            Only MultipleTestingMethod.BENJAMINI_HOCHBERG is
            supported in Version 1.
        significance_threshold: The adjusted-p-value (FDR) threshold
            used to annotate `significant`.

    Returns:
        A new tuple of GeneResult, in the same order as
        `gene_results`, with adjusted_p_value and significant
        populated for every gene that had a valid raw p-value, and
        left as None for every gene that did not.

    Raises:
        MultipleTestingError: If `multiple_testing_method` is not
            Benjamini-Hochberg, or if no gene in `gene_results` has a
            valid raw p-value to correct (a dataset-level statistical
            failure, per specification Section 10.3).
    """
    if multiple_testing_method != MultipleTestingMethod.BENJAMINI_HOCHBERG:
        raise MultipleTestingError(
            "Only Benjamini-Hochberg multiple-testing correction is "
            f"supported in Version 1; got "
            f"'{multiple_testing_method.value}'."
        )

    testable_indices = [
        index
        for index, result in enumerate(gene_results)
        if result.raw_p_value is not None
    ]

    if not testable_indices:
        raise MultipleTestingError(
            "No gene has a valid raw p-value; multiple-testing "
            "correction cannot be performed. This indicates a "
            "dataset-level statistical failure upstream, not a "
            "multiple-testing configuration problem."
        )

    raw_p_values = np.array(
        [gene_results[index].raw_p_value for index in testable_indices],
        dtype=float,
    )

    adjusted_p_values = false_discovery_control(raw_p_values, method="bh")

    adjusted_by_index: dict[int, float] = {
        index: float(adjusted_p_value)
        for index, adjusted_p_value in zip(testable_indices, adjusted_p_values)
    }

    corrected_results: list[GeneResult] = []
    for index, result in enumerate(gene_results):
        adjusted_p_value = adjusted_by_index.get(index)
        significant = (
            adjusted_p_value < significance_threshold
            if adjusted_p_value is not None
            else None
        )
        corrected_results.append(
            GeneResult(
                gene_id=result.gene_id,
                tested=result.tested,
                effect_size=result.effect_size,
                effect_size_label=result.effect_size_label,
                raw_p_value=result.raw_p_value,
                adjusted_p_value=adjusted_p_value,
                significant=significant,
                missing_reason=result.missing_reason,
            )
        )

    return tuple(corrected_results)