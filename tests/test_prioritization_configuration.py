from pathlib import Path

import pytest

from src.biomarker_prioritization.configuration import (
    InvalidPrioritizationConfigurationError,
    build_configuration,
)
from src.biomarker_prioritization.models import (
    CANCER_ASSOCIATION_WEIGHT,
    CLINICAL_WEIGHT,
    CROSS_CANCER_WEIGHT,
    DE_WEIGHT,
    PrioritizationConfiguration,
)


def test_build_configuration_accepts_string_paths():
    config = build_configuration(
        analysis_id="brca_phase5_v1",
        cancer_cohort="TCGA-BRCA",
        phase3_results_path="results/phase3/brca.csv",
        phase4_evidence_path="results/phase4/evidence.csv",
        phase4_metadata_path="results/phase4/metadata.json",
        output_dir="results/phase5/TCGA-BRCA",
    )

    assert isinstance(config.phase3_results_path, Path)
    assert isinstance(config.phase4_evidence_path, Path)
    assert isinstance(config.phase4_metadata_path, Path)
    assert isinstance(config.output_dir, Path)


def test_configuration_uses_fixed_v1_weights():
    config = build_configuration(
        analysis_id="brca_phase5_v1",
        cancer_cohort="TCGA-BRCA",
        phase3_results_path="phase3.csv",
        phase4_evidence_path="evidence.csv",
        phase4_metadata_path="metadata.json",
        output_dir="out",
    )

    assert config.de_weight == DE_WEIGHT
    assert config.cancer_association_weight == CANCER_ASSOCIATION_WEIGHT
    assert config.clinical_weight == CLINICAL_WEIGHT
    assert config.cross_cancer_weight == CROSS_CANCER_WEIGHT

    assert (
        config.de_weight
        + config.cancer_association_weight
        + config.clinical_weight
        + config.cross_cancer_weight
    ) == pytest.approx(1.0)


def test_configuration_rejects_empty_analysis_id():
    with pytest.raises(
        InvalidPrioritizationConfigurationError
    ):
        build_configuration(
            analysis_id="",
            cancer_cohort="TCGA-BRCA",
            phase3_results_path="phase3.csv",
            phase4_evidence_path="evidence.csv",
            phase4_metadata_path="metadata.json",
            output_dir="out",
        )


def test_configuration_rejects_empty_cancer_cohort():
    with pytest.raises(
        InvalidPrioritizationConfigurationError
    ):
        build_configuration(
            analysis_id="phase5",
            cancer_cohort="",
            phase3_results_path="phase3.csv",
            phase4_evidence_path="evidence.csv",
            phase4_metadata_path="metadata.json",
            output_dir="out",
        )


def test_configuration_rejects_empty_path():
    with pytest.raises(
        InvalidPrioritizationConfigurationError
    ):
        build_configuration(
            analysis_id="phase5",
            cancer_cohort="TCGA-BRCA",
            phase3_results_path="",
            phase4_evidence_path="evidence.csv",
            phase4_metadata_path="metadata.json",
            output_dir="out",
        )


def test_configuration_rejects_non_path_value():
    with pytest.raises(
        InvalidPrioritizationConfigurationError
    ):
        build_configuration(
            analysis_id="phase5",
            cancer_cohort="TCGA-BRCA",
            phase3_results_path=123,
            phase4_evidence_path="evidence.csv",
            phase4_metadata_path="metadata.json",
            output_dir="out",
        )


def test_model_rejects_custom_v1_weight():
    with pytest.raises(ValueError):
        PrioritizationConfiguration(
            analysis_id="phase5",
            cancer_cohort="TCGA-BRCA",
            phase3_results_path=Path("phase3.csv"),
            phase4_evidence_path=Path("evidence.csv"),
            phase4_metadata_path=Path("metadata.json"),
            output_dir=Path("out"),
            de_weight=0.40,
        )