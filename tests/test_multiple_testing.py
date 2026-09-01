"""Tests for src.differential_expression.multiple_testing."""

import pytest
from statistics import NormalDist

from src.differential_expression.exceptions import MultipleTestingError
from src.differential_expression.models import GeneResult, MultipleTestingMethod
from src.differential_expression.multiple_testing import (
    apply_multiple_testing_correction,
)


def _gene_result(gene_id, raw_p_value, missing_reason=None):
    return GeneResult(
        gene_id=gene_id,
        tested=True,
        effect_size=1.0 if raw_p_value is not None else None,
        effect_size_label="mean_difference" if raw_p_value is not None else None,
        raw_p_value=raw_p_value,
        adjusted_p_value=None,
        significant=None,
        missing_reason=missing_reason,
    )


def test_bh_correction_matches_known_values():
    # Known BH example: p-values 0.01, 0.02, 0.03, 0.04, 0.05 for 5
    # genes. BH adjusted p-values (monotone, non-decreasing from the
    # largest rank downward) are: 0.05, 0.05, 0.05, 0.05, 0.05
    # (each raw_p * n / rank, then cumulative min from the top).
    results = tuple(
        _gene_result(f"G{i}", p)
        for i, p in enumerate([0.01, 0.02, 0.03, 0.04, 0.05], start=1)
    )
    corrected = apply_multiple_testing_correction(
        results, MultipleTestingMethod.BENJAMINI_HOCHBERG, significance_threshold=0.05
    )
    adjusted = [gene.adjusted_p_value for gene in corrected]
    assert adjusted == pytest.approx([0.05, 0.05, 0.05, 0.05, 0.05])


def test_gene_identity_preserved():
    results = (
        _gene_result("GENE_X", 0.01),
        _gene_result("GENE_Y", 0.5),
    )
    corrected = apply_multiple_testing_correction(
        results, MultipleTestingMethod.BENJAMINI_HOCHBERG, significance_threshold=0.05
    )
    assert [gene.gene_id for gene in corrected] == ["GENE_X", "GENE_Y"]


def test_raw_p_values_preserved():
    results = (_gene_result("G1", 0.01), _gene_result("G2", 0.5))
    corrected = apply_multiple_testing_correction(
        results, MultipleTestingMethod.BENJAMINI_HOCHBERG, significance_threshold=0.05
    )
    assert corrected[0].raw_p_value == 0.01
    assert corrected[1].raw_p_value == 0.5


def test_missing_p_values_excluded_from_correction_but_retained():
    results = (
        _gene_result("G1", 0.01),
        _gene_result("G2", None, missing_reason="Zero-variance gene."),
        _gene_result("G3", 0.5),
    )
    corrected = apply_multiple_testing_correction(
        results, MultipleTestingMethod.BENJAMINI_HOCHBERG, significance_threshold=0.05
    )
    assert corrected[1].gene_id == "G2"
    assert corrected[1].adjusted_p_value is None
    assert corrected[1].significant is None
    assert corrected[1].missing_reason == "Zero-variance gene."
    assert corrected[0].adjusted_p_value is not None
    assert corrected[2].adjusted_p_value is not None


def test_significance_threshold_applied():
    results = (_gene_result("G1", 0.001), _gene_result("G2", 0.9))
    corrected = apply_multiple_testing_correction(
        results, MultipleTestingMethod.BENJAMINI_HOCHBERG, significance_threshold=0.05
    )
    significant_by_gene = {gene.gene_id: gene.significant for gene in corrected}
    assert significant_by_gene["G1"] is True
    assert significant_by_gene["G2"] is False


def test_deterministic_output():
    results = tuple(_gene_result(f"G{i}", 0.01 * i) for i in range(1, 11))
    first = apply_multiple_testing_correction(
        results, MultipleTestingMethod.BENJAMINI_HOCHBERG, significance_threshold=0.05
    )
    second = apply_multiple_testing_correction(
        results, MultipleTestingMethod.BENJAMINI_HOCHBERG, significance_threshold=0.05
    )
    assert [g.adjusted_p_value for g in first] == [g.adjusted_p_value for g in second]


def test_all_missing_p_values_raises():
    results = (
        _gene_result("G1", None, missing_reason="failed"),
        _gene_result("G2", None, missing_reason="failed"),
    )
    with pytest.raises(MultipleTestingError):
        apply_multiple_testing_correction(
            results, MultipleTestingMethod.BENJAMINI_HOCHBERG, significance_threshold=0.05
        )


def test_unsupported_correction_method_rejected():
    results = (_gene_result("G1", 0.01),)

    class _FakeMethod:
        value = "bonferroni"

    with pytest.raises(MultipleTestingError):
        apply_multiple_testing_correction(
            results, _FakeMethod(), significance_threshold=0.05  # type: ignore[arg-type]
        )