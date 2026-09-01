from unittest.mock import patch

from src.evidence_integration.models import EvidenceType
from src.evidence_integration.ncbi_gene import (
    collect_functional_evidence,
    fetch_gene_summary,
    find_human_gene_id,
)


def test_find_human_gene_id_returns_first_id():
    payload = {
        "esearchresult": {
            "idlist": ["7157"],
        }
    }

    with patch(
        "src.evidence_integration.ncbi_gene._get_json",
        return_value=payload,
    ):
        result = find_human_gene_id("TP53")

    assert result == "7157"


def test_find_human_gene_id_returns_none_when_missing():
    payload = {
        "esearchresult": {
            "idlist": [],
        }
    }

    with patch(
        "src.evidence_integration.ncbi_gene._get_json",
        return_value=payload,
    ):
        result = find_human_gene_id("UNKNOWN_GENE")

    assert result is None


def test_fetch_gene_summary_returns_record():
    payload = {
        "result": {
            "7157": {
                "uid": "7157",
                "name": "TP53",
                "description": "tumor protein p53",
                "summary": "TP53 participates in cellular stress responses.",
                "organism": {
                    "scientificname": "Homo sapiens",
                },
            }
        }
    }

    with patch(
        "src.evidence_integration.ncbi_gene._get_json",
        return_value=payload,
    ):
        result = fetch_gene_summary("7157")

    assert result["name"] == "TP53"


def test_collect_functional_evidence_constructs_record():
    summary = {
        "uid": "7157",
        "name": "TP53",
        "description": "tumor protein p53",
        "summary": "TP53 participates in cellular stress responses.",
        "organism": {
            "scientificname": "Homo sapiens",
        },
    }

    with (
        patch(
            "src.evidence_integration.ncbi_gene.find_human_gene_id",
            return_value="7157",
        ),
        patch(
            "src.evidence_integration.ncbi_gene.fetch_gene_summary",
            return_value=summary,
        ),
    ):
        record = collect_functional_evidence(
            "TP53",
            retrieved_at="2026-08-26T12:00:00+05:30",
        )

    assert record is not None
    assert record.gene_id == "TP53"
    assert record.evidence_type is EvidenceType.FUNCTIONAL
    assert record.source == "NCBI Gene"
    assert record.evidence_id == "NCBI_GENE:7157"


def test_missing_gene_returns_no_evidence():
    with patch(
        "src.evidence_integration.ncbi_gene.find_human_gene_id",
        return_value=None,
    ):
        result = collect_functional_evidence(
            "UNKNOWN",
            retrieved_at="2026-08-26T12:00:00+05:30",
        )

    assert result is None


def test_missing_summary_text_returns_no_evidence():
    summary = {
        "uid": "7157",
        "name": "TP53",
        "description": "",
        "summary": "",
        "organism": {
            "scientificname": "Homo sapiens",
        },
    }

    with (
        patch(
            "src.evidence_integration.ncbi_gene.find_human_gene_id",
            return_value="7157",
        ),
        patch(
            "src.evidence_integration.ncbi_gene.fetch_gene_summary",
            return_value=summary,
        ),
    ):
        result = collect_functional_evidence(
            "TP53",
            retrieved_at="2026-08-26T12:00:00+05:30",
        )

    assert result is None


def test_non_human_result_is_rejected():
    summary = {
        "uid": "123",
        "name": "TP53",
        "description": "example",
        "summary": "Example summary.",
        "organism": {
            "scientificname": "Mus musculus",
        },
    }

    with (
        patch(
            "src.evidence_integration.ncbi_gene.find_human_gene_id",
            return_value="123",
        ),
        patch(
            "src.evidence_integration.ncbi_gene.fetch_gene_summary",
            return_value=summary,
        ),
    ):
        result = collect_functional_evidence(
            "TP53",
            retrieved_at="2026-08-26T12:00:00+05:30",
        )

    assert result is None


def test_symbol_mismatch_is_rejected():
    summary = {
        "uid": "7157",
        "name": "OTHER",
        "description": "example",
        "summary": "Example summary.",
        "organism": {
            "scientificname": "Homo sapiens",
        },
    }

    with (
        patch(
            "src.evidence_integration.ncbi_gene.find_human_gene_id",
            return_value="7157",
        ),
        patch(
            "src.evidence_integration.ncbi_gene.fetch_gene_summary",
            return_value=summary,
        ),
    ):
        result = collect_functional_evidence(
            "TP53",
            retrieved_at="2026-08-26T12:00:00+05:30",
        )

    assert result is None