import json

import pandas as pd
import pytest

from src.biomarker_prioritization.analysis import (
    run_prioritization_analysis,
)
from src.biomarker_prioritization.configuration import (
    build_configuration,
)
from src.biomarker_prioritization.models import (
    EvidenceAvailability,
)


def _write_analysis_inputs(
    tmp_path,
    *,
    include_unresolved=False,
):
    phase3_path = tmp_path / "phase3.csv"
    evidence_path = tmp_path / "evidence.csv"
    metadata_path = tmp_path / "metadata.json"

    phase3_rows = [
        {
            "gene_id": "TP53",
            "tested": True,
            "effect_size": 4.0,
            "effect_size_label": "log2_fold_change",
            "adjusted_p_value": 0.001,
            "significant": True,
        },
        {
            "gene_id": "EGFR",
            "tested": True,
            "effect_size": 2.0,
            "effect_size_label": "log2_fold_change",
            "adjusted_p_value": 0.01,
            "significant": True,
        },
        {
            "gene_id": "GENE_NO_EVIDENCE",
            "tested": True,
            "effect_size": 1.0,
            "effect_size_label": "log2_fold_change",
            "adjusted_p_value": 0.02,
            "significant": True,
        },
    ]

    if include_unresolved:
        phase3_rows.append(
            {
                "gene_id": "?|90288",
                "tested": True,
                "effect_size": 3.0,
                "effect_size_label": "log2_fold_change",
                "adjusted_p_value": 0.0001,
                "significant": True,
            }
        )

    pd.DataFrame(
        phase3_rows
    ).to_csv(
        phase3_path,
        index=False,
    )

    evidence_rows = [
        {
            "gene_id": "TP53",
            "cancer_cohort": "TCGA-BRCA",
            "evidence_type": "cancer_association",
            "source": "Open Targets Platform",
            "source_version": "26.06",
            "evidence_id": "MONDO_TEST",
            "description": (
                "breast cancer; Open Targets association "
                "score=0.800000; evidence count=20"
            ),
            "retrieved_at": "2026-08-27T00:00:00+00:00",
            "source_url": "",
        },
        {
            "gene_id": "TP53",
            "cancer_cohort": "TCGA-BRCA",
            "evidence_type": "clinical",
            "source": "Human Protein Atlas",
            "source_version": "25.1",
            "evidence_id": "HPA_TEST",
            "description": (
                "validated prognostic - unfavorable; "
                "reported value=1e-4"
            ),
            "retrieved_at": "2026-08-27T00:00:00+00:00",
            "source_url": "",
        },
        {
            "gene_id": "TP53",
            "cancer_cohort": "TCGA-BRCA",
            "evidence_type": "cross_cancer",
            "source": "BLC Mark Phase 3",
            "source_version": "1.0",
            "evidence_id": "CROSS_CANCER:TP53",
            "description": (
                "Significant differential-expression candidate "
                "in 3 cohorts: TCGA-BRCA, TCGA-COAD, TCGA-LUAD"
            ),
            "retrieved_at": "2026-08-27T00:00:00+00:00",
            "source_url": "",
        },
        {
            "gene_id": "TP53",
            "cancer_cohort": "TCGA-BRCA",
            "evidence_type": "functional",
            "source": "NCBI Gene",
            "source_version": "snapshot",
            "evidence_id": "NCBI_GENE:7157",
            "description": "tumor protein p53",
            "retrieved_at": "2026-08-27T00:00:00+00:00",
            "source_url": "",
        },
        {
            "gene_id": "TP53",
            "cancer_cohort": "TCGA-BRCA",
            "evidence_type": "pathway",
            "source": "Reactome",
            "source_version": "97",
            "evidence_id": "R-HSA-1",
            "description": "Pathway one",
            "retrieved_at": "2026-08-27T00:00:00+00:00",
            "source_url": "",
        },
        {
            "gene_id": "TP53",
            "cancer_cohort": "TCGA-BRCA",
            "evidence_type": "pathway",
            "source": "Reactome",
            "source_version": "97",
            "evidence_id": "R-HSA-2",
            "description": "Pathway two",
            "retrieved_at": "2026-08-27T00:00:00+00:00",
            "source_url": "",
        },
        {
            "gene_id": "EGFR",
            "cancer_cohort": "TCGA-BRCA",
            "evidence_type": "cancer_association",
            "source": "Open Targets Platform",
            "source_version": "26.06",
            "evidence_id": "MONDO_TEST",
            "description": (
                "breast cancer; Open Targets association "
                "score=0.400000; evidence count=10"
            ),
            "retrieved_at": "2026-08-27T00:00:00+00:00",
            "source_url": "",
        },
        {
            "gene_id": "EGFR",
            "cancer_cohort": "TCGA-BRCA",
            "evidence_type": "clinical",
            "source": "Human Protein Atlas",
            "source_version": "25.1",
            "evidence_id": "HPA_EGFR",
            "description": (
                "potential prognostic - favorable; "
                "reported value=0.02"
            ),
            "retrieved_at": "2026-08-27T00:00:00+00:00",
            "source_url": "",
        },
        {
            "gene_id": "EGFR",
            "cancer_cohort": "TCGA-BRCA",
            "evidence_type": "cross_cancer",
            "source": "BLC Mark Phase 3",
            "source_version": "1.0",
            "evidence_id": "CROSS_CANCER:EGFR",
            "description": (
                "Significant differential-expression candidate "
                "in 2 cohorts: TCGA-BRCA, TCGA-LUAD"
            ),
            "retrieved_at": "2026-08-27T00:00:00+00:00",
            "source_url": "",
        },
    ]

    pd.DataFrame(
        evidence_rows
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
                "candidate_count": len(
                    phase3_rows
                ),
                "ranking_performed": False,
                "unresolved_identifier_count": (
                    1 if include_unresolved else 0
                ),
            }
        ),
        encoding="utf-8",
    )

    return build_configuration(
        analysis_id="phase5_test",
        cancer_cohort="TCGA-BRCA",
        phase3_results_path=phase3_path,
        phase4_evidence_path=evidence_path,
        phase4_metadata_path=metadata_path,
        output_dir=tmp_path / "out",
    )


def test_analysis_scores_and_ranks_candidates(
    tmp_path,
):
    configuration = _write_analysis_inputs(
        tmp_path
    )

    results = run_prioritization_analysis(
        configuration
    )

    assert len(results) == 3

    assert results[0].gene_id == "TP53"
    assert results[0].rank == 1

    assert results[1].gene_id == "EGFR"
    assert results[1].rank == 2

    assert (
        results[2].gene_id
        == "GENE_NO_EVIDENCE"
    )
    assert results[2].rank == 3


def test_analysis_preserves_context_only_evidence(
    tmp_path,
):
    configuration = _write_analysis_inputs(
        tmp_path
    )

    results = run_prioritization_analysis(
        configuration
    )

    tp53 = next(
        item
        for item in results
        if item.gene_id == "TP53"
    )

    assert (
        tp53.raw_input.functional_description
        == "tumor protein p53"
    )

    assert tp53.raw_input.pathway_count == 2


def test_no_evidence_means_no_positive_support(
    tmp_path,
):
    configuration = _write_analysis_inputs(
        tmp_path
    )

    results = run_prioritization_analysis(
        configuration
    )

    item = next(
        result
        for result in results
        if result.gene_id
        == "GENE_NO_EVIDENCE"
    )

    assert (
        item.raw_input.cancer_association_availability
        == EvidenceAvailability.NO_SUPPORT
    )

    assert (
        item.raw_input.clinical_availability
        == EvidenceAvailability.NO_SUPPORT
    )

    assert item.component_scores.cancer_association_score == 0.0
    assert item.component_scores.clinical_score == 0.0
    assert item.component_scores.cross_cancer_score == 0.0

    assert item.final_score is not None


def test_unresolved_candidate_is_retained_without_final_score(
    tmp_path,
):
    configuration = _write_analysis_inputs(
        tmp_path,
        include_unresolved=True,
    )

    results = run_prioritization_analysis(
        configuration
    )

    unresolved = next(
        item
        for item in results
        if item.gene_id == "?|90288"
    )

    assert (
        unresolved.raw_input.cancer_association_availability
        == EvidenceAvailability.UNAVAILABLE
    )

    assert (
        unresolved.raw_input.clinical_availability
        == EvidenceAvailability.UNAVAILABLE
    )

    assert unresolved.final_score is None
    assert unresolved.rank is None

    assert results[-1].gene_id == "?|90288"


def test_analysis_uses_absolute_effect_size_for_de_rank(
    tmp_path,
):
    configuration = _write_analysis_inputs(
        tmp_path
    )

    phase3 = pd.read_csv(
        configuration.phase3_results_path
    )

    phase3.loc[
        phase3["gene_id"] == "TP53",
        "effect_size",
    ] = -4.0

    phase3.to_csv(
        configuration.phase3_results_path,
        index=False,
    )

    results = run_prioritization_analysis(
        configuration
    )

    tp53 = next(
        item
        for item in results
        if item.gene_id == "TP53"
    )

    assert tp53.component_scores.de_score == pytest.approx(
        1.0
    )


def test_analysis_rejects_malformed_open_targets_description(
    tmp_path,
):
    configuration = _write_analysis_inputs(
        tmp_path
    )

    evidence = pd.read_csv(
        configuration.phase4_evidence_path
    )

    mask = (
        evidence["evidence_type"]
        == "cancer_association"
    )

    evidence.loc[
        mask,
        "description",
    ] = "malformed association record"

    evidence.to_csv(
        configuration.phase4_evidence_path,
        index=False,
    )

    with pytest.raises(ValueError):
        run_prioritization_analysis(
            configuration
        )


def test_analysis_rejects_multiple_scoring_records_of_same_type(
    tmp_path,
):
    configuration = _write_analysis_inputs(
        tmp_path
    )

    evidence = pd.read_csv(
        configuration.phase4_evidence_path
    )

    duplicate = evidence[
        (
            evidence["gene_id"] == "TP53"
        )
        & (
            evidence["evidence_type"]
            == "clinical"
        )
    ].copy()

    duplicate.loc[
        :,
        "evidence_id",
    ] = "SECOND_CLINICAL_RECORD"

    evidence = pd.concat(
        [
            evidence,
            duplicate,
        ],
        ignore_index=True,
    )

    evidence.to_csv(
        configuration.phase4_evidence_path,
        index=False,
    )

    with pytest.raises(ValueError):
        run_prioritization_analysis(
            configuration
        )