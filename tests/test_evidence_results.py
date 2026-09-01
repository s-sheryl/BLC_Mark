import csv
import json

import pytest

from src.evidence_integration.aggregation import GeneEvidenceProfile
from src.evidence_integration.models import EvidenceRecord, EvidenceType
from src.evidence_integration.qc import build_phase4_qc_report
from src.evidence_integration.results import (
    write_evidence_profiles,
    write_json,
)


def make_profile(
    gene_id: str,
    records: tuple[EvidenceRecord, ...],
) -> GeneEvidenceProfile:
    return GeneEvidenceProfile(
        gene_id=gene_id,
        cancer_cohort="TCGA-BRCA",
        evidence_records=records,
    )


def test_write_evidence_profiles(tmp_path):
    record = EvidenceRecord(
        gene_id="TP53",
        evidence_type=EvidenceType.FUNCTIONAL,
        source="NCBI Gene",
        source_version=None,
        evidence_id="NCBI_GENE:7157",
        description="Functional evidence.",
    )

    profile = make_profile(
        "TP53",
        (record,),
    )

    output = tmp_path / "evidence.csv"

    write_evidence_profiles(
        [profile],
        output,
    )

    with output.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 1
    assert rows[0]["gene_id"] == "TP53"
    assert rows[0]["evidence_type"] == "functional"
    assert rows[0]["source"] == "NCBI Gene"


def test_empty_profile_writes_header_only(tmp_path):
    profile = make_profile(
        "TP53",
        (),
    )

    output = tmp_path / "evidence.csv"

    write_evidence_profiles(
        [profile],
        output,
    )

    with output.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))

    assert rows == []


def test_write_json(tmp_path):
    output = tmp_path / "metadata.json"

    write_json(
        {"b": 2, "a": 1},
        output,
    )

    with output.open(
        encoding="utf-8",
    ) as handle:
        payload = json.load(handle)

    assert payload == {
        "a": 1,
        "b": 2,
    }


def test_qc_report_counts():
    functional = EvidenceRecord(
        gene_id="TP53",
        evidence_type=EvidenceType.FUNCTIONAL,
        source="NCBI Gene",
        description="Functional.",
    )

    pathway = EvidenceRecord(
        gene_id="TP53",
        evidence_type=EvidenceType.PATHWAY,
        source="Reactome",
        description="Pathway.",
    )

    profiles = [
        make_profile(
            "TP53",
            (
                functional,
                pathway,
            ),
        ),
        make_profile(
            "BRCA1",
            (),
        ),
    ]

    qc = build_phase4_qc_report(
        profiles,
        candidate_count=2,
        unresolved_identifier_count=1,
    )

    assert qc["candidate_count"] == 2
    assert qc["genes_with_evidence"] == 1
    assert qc["genes_without_evidence"] == 1
    assert qc["unresolved_identifier_count"] == 1
    assert qc["total_evidence_records"] == 2

    assert qc["evidence_type_counts"] == {
        "functional": 1,
        "pathway": 1,
    }


def test_qc_profile_count_mismatch_is_rejected():
    with pytest.raises(
        ValueError,
        match="Profile count",
    ):
        build_phase4_qc_report(
            [],
            candidate_count=1,
            unresolved_identifier_count=0,
        )