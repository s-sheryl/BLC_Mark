from unittest.mock import patch

import pytest

from src.evidence_integration.open_targets_identifiers import (
    resolve_gene_symbol_to_ensembl,
)


def test_exact_target_match_returns_ensembl_id():
    payload = {
        "data": {
            "search": {
                "hits": [
                    {
                        "id": "ENSG00000141510",
                        "entity": "target",
                        "name": "TP53",
                    }
                ]
            }
        }
    }

    with patch(
        "src.evidence_integration."
        "open_targets_identifiers._post_graphql",
        return_value=payload,
    ):
        result = resolve_gene_symbol_to_ensembl("TP53")

    assert result == "ENSG00000141510"


def test_non_exact_symbol_is_rejected():
    payload = {
        "data": {
            "search": {
                "hits": [
                    {
                        "id": "ENSG00000141510",
                        "entity": "target",
                        "name": "TP53BP1",
                    }
                ]
            }
        }
    }

    with patch(
        "src.evidence_integration."
        "open_targets_identifiers._post_graphql",
        return_value=payload,
    ):
        result = resolve_gene_symbol_to_ensembl("TP53")

    assert result is None


def test_non_target_entity_is_rejected():
    payload = {
        "data": {
            "search": {
                "hits": [
                    {
                        "id": "MONDO_123",
                        "entity": "disease",
                        "name": "TP53",
                    }
                ]
            }
        }
    }

    with patch(
        "src.evidence_integration."
        "open_targets_identifiers._post_graphql",
        return_value=payload,
    ):
        result = resolve_gene_symbol_to_ensembl("TP53")

    assert result is None


def test_multiple_exact_matches_are_rejected():
    payload = {
        "data": {
            "search": {
                "hits": [
                    {
                        "id": "ENSG00000111111",
                        "entity": "target",
                        "name": "GENE1",
                    },
                    {
                        "id": "ENSG00000222222",
                        "entity": "target",
                        "name": "GENE1",
                    },
                ]
            }
        }
    }

    with patch(
        "src.evidence_integration."
        "open_targets_identifiers._post_graphql",
        return_value=payload,
    ):
        result = resolve_gene_symbol_to_ensembl("GENE1")

    assert result is None


def test_no_hits_returns_none():
    payload = {
        "data": {
            "search": {
                "hits": [],
            }
        }
    }

    with patch(
        "src.evidence_integration."
        "open_targets_identifiers._post_graphql",
        return_value=payload,
    ):
        result = resolve_gene_symbol_to_ensembl("UNKNOWN")

    assert result is None


@pytest.mark.parametrize(
    "gene_symbol",
    ["", " ", "\t"],
)
def test_blank_gene_symbol_is_rejected(gene_symbol):
    with pytest.raises(
        ValueError,
        match="gene_symbol",
    ):
        resolve_gene_symbol_to_ensembl(
            gene_symbol
        )