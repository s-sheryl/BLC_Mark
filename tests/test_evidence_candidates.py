from pathlib import Path

import pandas as pd

from src.evidence_integration.candidates import extract_significant_candidates


def write_results(tmp_path: Path, rows: list[dict]) -> Path:
    path = tmp_path / "phase3.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def make_row(
    gene_id: str,
    adjusted_p_value: float,
    significant: bool,
) -> dict:
    return {
        "gene_id": gene_id,
        "tested": True,
        "effect_size": 1.0,
        "effect_size_label": "difference_in_group_means",
        "raw_p_value": adjusted_p_value / 2,
        "adjusted_p_value": adjusted_p_value,
        "significant": significant,
        "missing_reason": None,
    }


def test_only_significant_candidates_are_retained(tmp_path):
    rows = [
        make_row("TP53", 0.001, True),
        make_row("GENE2", 0.50, False),
    ]

    path = write_results(tmp_path, rows)

    candidates = extract_significant_candidates(path)

    assert list(candidates["gene_id"]) == ["TP53"]


def test_candidates_are_sorted_deterministically(tmp_path):
    rows = [
        make_row("GENE_B", 0.01, True),
        make_row("GENE_C", 0.001, True),
        make_row("GENE_A", 0.01, True),
    ]

    path = write_results(tmp_path, rows)

    candidates = extract_significant_candidates(path)

    assert list(candidates["gene_id"]) == [
        "GENE_C",
        "GENE_A",
        "GENE_B",
    ]


def test_no_additional_effect_size_filter_is_applied(tmp_path):
    row = make_row("LOW_EFFECT", 0.001, True)
    row["effect_size"] = 0.0001

    path = write_results(tmp_path, [row])

    candidates = extract_significant_candidates(path)

    assert len(candidates) == 1
    assert candidates.loc[0, "gene_id"] == "LOW_EFFECT"


def test_empty_candidate_set_is_allowed(tmp_path):
    rows = [
        make_row("GENE1", 0.50, False),
        make_row("GENE2", 0.75, False),
    ]

    path = write_results(tmp_path, rows)

    candidates = extract_significant_candidates(path)

    assert candidates.empty