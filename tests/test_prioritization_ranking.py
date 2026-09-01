from src.biomarker_prioritization.models import (
    ComponentScores,
    EvidenceAvailability,
    PrioritizationInput,
    PrioritizedBiomarker,
)
from src.biomarker_prioritization.ranking import rank_biomarkers


def _make_biomarker(
    gene_id: str,
    final_score: float | None,
) -> PrioritizedBiomarker:
    raw_input = PrioritizationInput(
        gene_id=gene_id,
        cancer_cohort="TCGA-BRCA",
        effect_size=1.0,
        effect_size_label="log2_fold_change",
        adjusted_p_value=0.01,
        cancer_association_score=0.5,
        clinical_category="potential prognostic - favorable",
        clinical_direction="favorable",
        cross_cancer_cohort_count=2,
        cancer_association_availability=EvidenceAvailability.AVAILABLE,
        clinical_availability=EvidenceAvailability.AVAILABLE,
        cross_cancer_availability=EvidenceAvailability.AVAILABLE,
        functional_description=None,
        pathway_count=0,
    )

    component_scores = ComponentScores(
        de_score=0.5,
        cancer_association_score=0.5,
        clinical_score=0.5,
        cross_cancer_score=0.5,
    )

    return PrioritizedBiomarker(
        gene_id=gene_id,
        cancer_cohort="TCGA-BRCA",
        raw_input=raw_input,
        component_scores=component_scores,
        final_score=final_score,
        rank=None,
    )


def test_rank_biomarkers_orders_scores_descending():
    biomarkers = (
        _make_biomarker("GENE_A", 0.20),
        _make_biomarker("GENE_B", 0.90),
        _make_biomarker("GENE_C", 0.50),
    )

    ranked = rank_biomarkers(biomarkers)

    assert [item.gene_id for item in ranked] == [
        "GENE_B",
        "GENE_C",
        "GENE_A",
    ]

    assert [item.rank for item in ranked] == [
        1,
        2,
        3,
    ]


def test_rank_biomarkers_breaks_ties_by_gene_id():
    biomarkers = (
        _make_biomarker("TP53", 0.80),
        _make_biomarker("BRCA1", 0.80),
        _make_biomarker("EGFR", 0.80),
    )

    ranked = rank_biomarkers(biomarkers)

    assert [item.gene_id for item in ranked] == [
        "BRCA1",
        "EGFR",
        "TP53",
    ]

    assert [item.rank for item in ranked] == [
        1,
        2,
        3,
    ]


def test_unavailable_scores_are_placed_last():
    biomarkers = (
        _make_biomarker("GENE_A", None),
        _make_biomarker("GENE_B", 0.40),
        _make_biomarker("GENE_C", None),
        _make_biomarker("GENE_D", 0.70),
    )

    ranked = rank_biomarkers(biomarkers)

    assert [item.gene_id for item in ranked] == [
        "GENE_D",
        "GENE_B",
        "GENE_A",
        "GENE_C",
    ]

    assert [item.rank for item in ranked] == [
        1,
        2,
        None,
        None,
    ]


def test_unavailable_candidates_are_sorted_deterministically():
    biomarkers = (
        _make_biomarker("ZZZ", None),
        _make_biomarker("AAA", None),
        _make_biomarker("MMM", None),
    )

    ranked = rank_biomarkers(biomarkers)

    assert [item.gene_id for item in ranked] == [
        "AAA",
        "MMM",
        "ZZZ",
    ]

    assert all(item.rank is None for item in ranked)


def test_empty_input_returns_empty_tuple():
    assert rank_biomarkers(()) == ()


def test_original_objects_are_not_mutated():
    biomarker = _make_biomarker("TP53", 0.90)

    ranked = rank_biomarkers((biomarker,))

    assert biomarker.rank is None
    assert ranked[0].rank == 1


def test_rank_biomarkers_rejects_invalid_item():
    try:
        rank_biomarkers(("not-a-biomarker",))
    except TypeError:
        pass
    else:
        raise AssertionError(
            "Expected TypeError for invalid ranking input."
        )