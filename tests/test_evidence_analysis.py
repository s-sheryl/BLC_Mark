import pandas as pd

from src.evidence_integration.analysis import (
    Phase4AnalysisResult,
    run_phase4_analysis,
)
from src.evidence_integration.models import (
    EvidenceRecord,
    EvidenceType,
)
from src.evidence_integration.open_targets_local import (
    OpenTargetsAssociation,
)


def write_phase3_results(tmp_path):
    rows = [
        {
            "gene_id": "TP53",
            "tested": True,
            "effect_size": 1.0,
            "effect_size_label": "difference_in_group_means",
            "raw_p_value": 0.001,
            "adjusted_p_value": 0.01,
            "significant": True,
            "missing_reason": None,
        },
        {
            "gene_id": "?|10431",
            "tested": True,
            "effect_size": 1.0,
            "effect_size_label": "difference_in_group_means",
            "raw_p_value": 0.001,
            "adjusted_p_value": 0.01,
            "significant": True,
            "missing_reason": None,
        },
        {
            "gene_id": "NOT_SIGNIFICANT",
            "tested": True,
            "effect_size": 0.0,
            "effect_size_label": "difference_in_group_means",
            "raw_p_value": 0.5,
            "adjusted_p_value": 0.8,
            "significant": False,
            "missing_reason": None,
        },
    ]

    path = tmp_path / "phase3.csv"

    pd.DataFrame(rows).to_csv(
        path,
        index=False,
    )

    return path


def fake_functional(
    gene_symbol,
    retrieved_at,
):
    return EvidenceRecord(
        gene_id=gene_symbol,
        evidence_type=EvidenceType.FUNCTIONAL,
        source="NCBI Gene",
        evidence_id="NCBI_GENE:7157",
        description="Functional evidence.",
        retrieved_at=retrieved_at,
    )


def test_phase4_analysis_integrates_available_evidence(
    tmp_path,
):
    phase3_path = write_phase3_results(
        tmp_path
    )

    reactome_mapping = {
        "7157": [
            {
                "pathway_id": "R-HSA-123",
                "pathway_url": "https://reactome.org/123",
                "pathway_name": "DNA repair",
                "evidence_code": "TAS",
                "species": "Homo sapiens",
            }
        ]
    }

    hpa_index = {
        (
            "ENSG00000141510",
            "Breast Invasive Carcinoma (TCGA)",
        ): {
            "Gene": "ENSG00000141510",
            "Gene name": "TP53",
            "Cancer": "Breast Invasive Carcinoma (TCGA)",
            "potential prognostic - favorable": "",
            "unprognostic - favorable": "",
            "potential prognostic - unfavorable": "",
            "unprognostic - unfavorable": "0.01",
            "validated prognostic - favorable": "",
            "validated prognostic - unfavorable": "",
        }
    }

    open_targets_index = {
        (
            "MONDO_0007254",
            "ENSG00000141510",
        ): OpenTargetsAssociation(
            disease_id="MONDO_0007254",
            target_id="ENSG00000141510",
            association_score=0.75,
            evidence_count=10,
        )
    }

    cross_cancer_index = {
        "TP53": (
            "TCGA-BRCA",
            "TCGA-LUAD",
        ),
    }

    result = run_phase4_analysis(
        phase3_path,
        "TCGA-BRCA",
        cross_cancer_index=cross_cancer_index,
        reactome_mapping=reactome_mapping,
        hpa_index=hpa_index,
        open_targets_index=open_targets_index,
        retrieved_at="2026-08-26T13:30:00+05:30",
        reactome_version="97",
        hpa_version="25.1",
        open_targets_version="26.06",
        resolve_ncbi_gene_id=(
            lambda symbol: (
                "7157"
                if symbol == "TP53"
                else None
            )
        ),
        resolve_ensembl_id=(
            lambda symbol: (
                "ENSG00000141510"
                if symbol == "TP53"
                else None
            )
        ),
        collect_functional=fake_functional,
    )

    assert isinstance(
        result,
        Phase4AnalysisResult,
    )

    assert result.candidate_count == 2
    assert result.unresolved_identifier_count == 1
    assert len(result.profiles) == 2

    tp53 = next(
        profile
        for profile in result.profiles
        if profile.gene_id == "TP53"
    )

    evidence_types = {
        record.evidence_type
        for record in tp53.evidence_records
    }

    assert EvidenceType.FUNCTIONAL in evidence_types
    assert EvidenceType.PATHWAY in evidence_types
    assert (
        EvidenceType.CANCER_ASSOCIATION
        in evidence_types
    )
    assert EvidenceType.CLINICAL in evidence_types
    assert EvidenceType.CROSS_CANCER in evidence_types


def test_unresolved_identifier_is_preserved(
    tmp_path,
):
    phase3_path = write_phase3_results(
        tmp_path
    )

    result = run_phase4_analysis(
        phase3_path,
        "TCGA-BRCA",
        cross_cancer_index={},
        reactome_mapping={},
        hpa_index={},
        open_targets_index={},
        retrieved_at="2026-08-26T13:30:00+05:30",
        reactome_version="97",
        hpa_version="25.1",
        open_targets_version="26.06",
        resolve_ncbi_gene_id=lambda symbol: None,
        resolve_ensembl_id=lambda symbol: None,
        collect_functional=lambda *args: None,
    )

    unresolved = next(
        profile
        for profile in result.profiles
        if profile.gene_id == "?|10431"
    )

    assert unresolved.evidence_records == ()