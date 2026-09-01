from pathlib import Path

import pandas as pd
import pytest

from src.evidence_integration.exceptions import EvidenceInputError
from src.evidence_integration.validation import validate_phase3_results


VALID_ROW = {
    "gene_id": "TP53",
    "tested": True,
    "effect_size": 1.25,
    "effect_size_label": "difference_in_group_means",
    "raw_p_value": 0.001,
    "adjusted_p_value": 0.01,
    "significant": True,
    "missing_reason": None,
}


def write_csv(tmp_path: Path, rows: list[dict]) -> Path:
    path = tmp_path / "phase3.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_valid_phase3_results_are_loaded(tmp_path):
    path = write_csv(tmp_path, [VALID_ROW])

    dataframe = validate_phase3_results(path)

    assert len(dataframe) == 1
    assert dataframe.loc[0, "gene_id"] == "TP53"


def test_missing_file_is_rejected(tmp_path):
    path = tmp_path / "missing.csv"

    with pytest.raises(EvidenceInputError):
        validate_phase3_results(path)


def test_missing_required_column_is_rejected(tmp_path):
    row = VALID_ROW.copy()
    row.pop("adjusted_p_value")
    path = write_csv(tmp_path, [row])

    with pytest.raises(EvidenceInputError, match="adjusted_p_value"):
        validate_phase3_results(path)


def test_blank_gene_id_is_rejected(tmp_path):
    row = VALID_ROW.copy()
    row["gene_id"] = " "
    path = write_csv(tmp_path, [row])

    with pytest.raises(EvidenceInputError, match="blank gene_id"):
        validate_phase3_results(path)


def test_duplicate_gene_ids_are_rejected(tmp_path):
    row_two = VALID_ROW.copy()
    path = write_csv(tmp_path, [VALID_ROW, row_two])

    with pytest.raises(EvidenceInputError, match="duplicate gene_id"):
        validate_phase3_results(path)


def test_untested_gene_is_rejected(tmp_path):
    row = VALID_ROW.copy()
    row["tested"] = False
    path = write_csv(tmp_path, [row])

    with pytest.raises(EvidenceInputError, match="untested"):
        validate_phase3_results(path)