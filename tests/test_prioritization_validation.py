import json
from pathlib import Path

import pandas as pd
import pytest

from src.biomarker_prioritization.configuration import (
    build_configuration,
)
from src.biomarker_prioritization.exceptions import (
    PrioritizationValidationError,
)
from src.biomarker_prioritization.validation import (
    validate_prioritization_inputs,
)


def _write_valid_inputs(
    tmp_path: Path,
):
    phase3_path = tmp_path / "phase3.csv"
    evidence_path = tmp_path / "evidence.csv"
    metadata_path = tmp_path / "metadata.json"

    pd.DataFrame(
        [
            {
                "gene_id": "TP53",
                "tested": True,
                "effect_size": 2.0,
                "effect_size_label": "log2_fold_change",
                "adjusted_p_value": 0.001,
                "significant": True,
            },
            {
                "gene_id": "EGFR",
                "tested": True,
                "effect_size": 0.5,
                "effect_size_label": "log2_fold_change",
                "adjusted_p_value": 0.01,
                "significant": True,
            },
            {
                "gene_id": "NOTSIG",
                "tested": True,
                "effect_size": 0.1,
                "effect_size_label": "log2_fold_change",
                "adjusted_p_value": 0.5,
                "significant": False,
            },
        ]
    ).to_csv(
        phase3_path,
        index=False,
    )

    pd.DataFrame(
        [
            {
                "gene_id": "TP53",
                "cancer_cohort": "TCGA-BRCA",
                "evidence_type": "cancer_association",
                "source": "Open Targets Platform",
                "source_version": "26.06",
                "evidence_id": "MONDO_TEST",
                "description": "test evidence",
                "retrieved_at": "2026-08-27T00:00:00+00:00",
                "source_url": "",
            }
        ]
    ).to_csv(
        evidence_path,
        index=False,
    )

    metadata_path.write_text(
        json.dumps(
            {
                "phase": 4,
                "phase_name": "Evidence Integration",
                "cancer_cohort": "TCGA-BRCA",
                "candidate_count": 2,
                "ranking_performed": False,
                "unresolved_identifier_count": 0,
            }
        ),
        encoding="utf-8",
    )

    return (
        phase3_path,
        evidence_path,
        metadata_path,
    )


def _configuration(
    tmp_path: Path,
):
    phase3, evidence, metadata = (
        _write_valid_inputs(tmp_path)
    )

    return build_configuration(
        analysis_id="phase5_test",
        cancer_cohort="TCGA-BRCA",
        phase3_results_path=phase3,
        phase4_evidence_path=evidence,
        phase4_metadata_path=metadata,
        output_dir=tmp_path / "out",
    )


def test_validation_accepts_consistent_inputs(
    tmp_path,
):
    config = _configuration(tmp_path)

    candidates, evidence, metadata = (
        validate_prioritization_inputs(
            config
        )
    )

    assert len(candidates) == 2
    assert set(candidates["gene_id"]) == {
        "TP53",
        "EGFR",
    }
    assert len(evidence) == 1
    assert metadata["candidate_count"] == 2


def test_validation_rejects_missing_phase3_file(
    tmp_path,
):
    _, evidence, metadata = (
        _write_valid_inputs(tmp_path)
    )

    config = build_configuration(
        analysis_id="phase5_test",
        cancer_cohort="TCGA-BRCA",
        phase3_results_path=tmp_path / "missing.csv",
        phase4_evidence_path=evidence,
        phase4_metadata_path=metadata,
        output_dir=tmp_path / "out",
    )

    with pytest.raises(
        PrioritizationValidationError
    ):
        validate_prioritization_inputs(
            config
        )


def test_validation_rejects_wrong_metadata_phase(
    tmp_path,
):
    config = _configuration(tmp_path)

    metadata = json.loads(
        config.phase4_metadata_path.read_text(
            encoding="utf-8"
        )
    )

    metadata["phase"] = 3

    config.phase4_metadata_path.write_text(
        json.dumps(metadata),
        encoding="utf-8",
    )

    with pytest.raises(
        PrioritizationValidationError
    ):
        validate_prioritization_inputs(
            config
        )


def test_validation_rejects_cohort_mismatch(
    tmp_path,
):
    config = _configuration(tmp_path)

    metadata = json.loads(
        config.phase4_metadata_path.read_text(
            encoding="utf-8"
        )
    )

    metadata["cancer_cohort"] = (
        "TCGA-LUAD"
    )

    config.phase4_metadata_path.write_text(
        json.dumps(metadata),
        encoding="utf-8",
    )

    with pytest.raises(
        PrioritizationValidationError
    ):
        validate_prioritization_inputs(
            config
        )


def test_validation_rejects_candidate_count_mismatch(
    tmp_path,
):
    config = _configuration(tmp_path)

    metadata = json.loads(
        config.phase4_metadata_path.read_text(
            encoding="utf-8"
        )
    )

    metadata["candidate_count"] = 999

    config.phase4_metadata_path.write_text(
        json.dumps(metadata),
        encoding="utf-8",
    )

    with pytest.raises(
        PrioritizationValidationError
    ):
        validate_prioritization_inputs(
            config
        )


def test_validation_rejects_duplicate_phase3_candidates(
    tmp_path,
):
    config = _configuration(tmp_path)

    dataframe = pd.read_csv(
        config.phase3_results_path
    )

    duplicate = dataframe.iloc[[0]]

    dataframe = pd.concat(
        [
            dataframe,
            duplicate,
        ],
        ignore_index=True,
    )

    dataframe.to_csv(
        config.phase3_results_path,
        index=False,
    )

    metadata = json.loads(
        config.phase4_metadata_path.read_text(
            encoding="utf-8"
        )
    )

    metadata["candidate_count"] = 3

    config.phase4_metadata_path.write_text(
        json.dumps(metadata),
        encoding="utf-8",
    )

    with pytest.raises(
        PrioritizationValidationError
    ):
        validate_prioritization_inputs(
            config
        )


def test_validation_rejects_unexpected_evidence_gene(
    tmp_path,
):
    config = _configuration(tmp_path)

    evidence = pd.read_csv(
        config.phase4_evidence_path
    )

    extra = evidence.iloc[[0]].copy()
    extra.loc[:, "gene_id"] = "NOT_A_CANDIDATE"
    extra.loc[:, "evidence_id"] = "OTHER"

    evidence = pd.concat(
        [
            evidence,
            extra,
        ],
        ignore_index=True,
    )

    evidence.to_csv(
        config.phase4_evidence_path,
        index=False,
    )

    with pytest.raises(
        PrioritizationValidationError
    ):
        validate_prioritization_inputs(
            config
        )


def test_validation_rejects_duplicate_evidence_record(
    tmp_path,
):
    config = _configuration(tmp_path)

    evidence = pd.read_csv(
        config.phase4_evidence_path
    )

    evidence = pd.concat(
        [
            evidence,
            evidence.iloc[[0]],
        ],
        ignore_index=True,
    )

    evidence.to_csv(
        config.phase4_evidence_path,
        index=False,
    )

    with pytest.raises(
        PrioritizationValidationError
    ):
        validate_prioritization_inputs(
            config
        )