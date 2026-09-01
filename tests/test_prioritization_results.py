import csv
import json

from src.biomarker_prioritization.models import (
    ComponentScores,
    EvidenceAvailability,
    PrioritizationInput,
    PrioritizedBiomarker,
)
from src.biomarker_prioritization.results import (
    PRIORITIZATION_OUTPUT_COLUMNS,
    write_json,
    write_prioritization_results,
)


def _biomarker(
    gene_id: str,
    final_score: float | None,
    rank: int | None,
) -> PrioritizedBiomarker:
    raw_input = PrioritizationInput(
        gene_id=gene_id,
        cancer_cohort="TCGA-BRCA",
        effect_size=2.0,
        effect_size_label="log2_fold_change",
        adjusted_p_value=0.01,
        cancer_association_score=0.4,
        clinical_category="potential prognostic - favorable",
        clinical_direction="favorable",
        cross_cancer_cohort_count=2,
        cancer_association_availability=EvidenceAvailability.AVAILABLE,
        clinical_availability=EvidenceAvailability.AVAILABLE,
        cross_cancer_availability=EvidenceAvailability.AVAILABLE,
        functional_description="test gene",
        pathway_count=3,
    )

    scores = ComponentScores(
        de_score=0.8,
        cancer_association_score=0.4,
        clinical_score=0.5,
        cross_cancer_score=0.5,
    )

    return PrioritizedBiomarker(
        gene_id=gene_id,
        cancer_cohort="TCGA-BRCA",
        raw_input=raw_input,
        component_scores=scores,
        final_score=final_score,
        rank=rank,
    )


def test_write_prioritization_results(tmp_path):
    output_path = tmp_path / "prioritization_results.csv"

    biomarkers = [
        _biomarker("TP53", 0.8, 1),
        _biomarker("EGFR", 0.6, 2),
    ]

    returned = write_prioritization_results(
        biomarkers,
        output_path,
    )

    assert returned == output_path
    assert output_path.exists()

    with output_path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 2
    assert tuple(rows[0].keys()) == PRIORITIZATION_OUTPUT_COLUMNS
    assert rows[0]["gene_id"] == "TP53"
    assert rows[0]["rank"] == "1"
    assert rows[0]["functional_description"] == "test gene"


def test_results_are_written_in_rank_order(tmp_path):
    output_path = tmp_path / "results.csv"

    biomarkers = [
        _biomarker("EGFR", 0.6, 2),
        _biomarker("TP53", 0.8, 1),
    ]

    write_prioritization_results(
        biomarkers,
        output_path,
    )

    with output_path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))

    assert [row["gene_id"] for row in rows] == [
        "TP53",
        "EGFR",
    ]


def test_unavailable_score_written_last(tmp_path):
    output_path = tmp_path / "results.csv"

    biomarkers = [
        _biomarker("UNRESOLVED", None, None),
        _biomarker("TP53", 0.8, 1),
    ]

    write_prioritization_results(
        biomarkers,
        output_path,
    )

    with output_path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))

    assert rows[-1]["gene_id"] == "UNRESOLVED"
    assert rows[-1]["final_score"] == ""
    assert rows[-1]["rank"] == ""


def test_write_json_is_deterministic(tmp_path):
    output_path = tmp_path / "metadata.json"

    payload = {
        "z": 1,
        "a": 2,
    }

    write_json(
        payload,
        output_path,
    )

    text = output_path.read_text(
        encoding="utf-8"
    )

    assert text.endswith("\n")

    loaded = json.loads(text)

    assert loaded == payload
    assert text.index('"a"') < text.index('"z"')


def test_results_reject_invalid_item(tmp_path):
    output_path = tmp_path / "results.csv"

    try:
        write_prioritization_results(
            ["bad"],
            output_path,
        )
    except TypeError:
        pass
    else:
        raise AssertionError(
            "Expected TypeError."
        )
