import pytest

from src.evidence_integration.identifiers import (
    GeneIdentifier,
    normalize_gene_symbol,
)


def test_valid_gene_identifier_constructs():
    identifier = GeneIdentifier(
        original_id="TP53",
        normalized_symbol="TP53",
        resolvable=True,
    )

    assert identifier.original_id == "TP53"
    assert identifier.normalized_symbol == "TP53"
    assert identifier.resolvable is True


def test_whitespace_is_removed():
    identifier = GeneIdentifier(
        original_id=" TP53 ",
        normalized_symbol=" TP53 ",
        resolvable=True,
    )

    assert identifier.original_id == "TP53"
    assert identifier.normalized_symbol == "TP53"


def test_normalize_gene_symbol_preserves_gene_symbol():
    identifier = normalize_gene_symbol("TP53")

    assert identifier.original_id == "TP53"
    assert identifier.normalized_symbol == "TP53"
    assert identifier.resolvable is True
    assert identifier.unresolved_reason is None


def test_normalize_gene_symbol_does_not_change_case():
    identifier = normalize_gene_symbol("Tp53")

    assert identifier.normalized_symbol == "Tp53"


def test_normalize_gene_symbol_does_not_replace_aliases():
    identifier = normalize_gene_symbol("P53")

    assert identifier.normalized_symbol == "P53"


def test_legacy_unresolved_identifier_is_flagged():
    identifier = normalize_gene_symbol("?|10431")

    assert identifier.original_id == "?|10431"
    assert identifier.normalized_symbol is None
    assert identifier.resolvable is False
    assert (
        identifier.unresolved_reason
        == "Legacy unresolved TCGA/Xena identifier."
    )


def test_non_numeric_pipe_identifier_is_not_silently_rewritten():
    identifier = normalize_gene_symbol("?|ABC")

    assert identifier.original_id == "?|ABC"
    assert identifier.normalized_symbol == "?|ABC"
    assert identifier.resolvable is True


def test_unresolvable_identifier_requires_reason():
    with pytest.raises(
        ValueError,
        match="unresolved_reason",
    ):
        GeneIdentifier(
            original_id="?|10431",
            normalized_symbol=None,
            resolvable=False,
            unresolved_reason=None,
        )


def test_unresolvable_identifier_cannot_have_normalized_symbol():
    with pytest.raises(
        ValueError,
        match="must not have",
    ):
        GeneIdentifier(
            original_id="?|10431",
            normalized_symbol="10431",
            resolvable=False,
            unresolved_reason="Legacy unresolved identifier.",
        )


def test_resolvable_identifier_cannot_have_unresolved_reason():
    with pytest.raises(
        ValueError,
        match="cannot have",
    ):
        GeneIdentifier(
            original_id="TP53",
            normalized_symbol="TP53",
            resolvable=True,
            unresolved_reason="Should not exist.",
        )


@pytest.mark.parametrize(
    "value",
    ["", " ", "\t"],
)
def test_blank_gene_id_is_rejected(value):
    with pytest.raises(
        ValueError,
        match="gene_id",
    ):
        normalize_gene_symbol(value)


def test_invalid_gene_symbol_input_type_is_rejected():
    with pytest.raises(
        ValueError,
        match="gene_id",
    ):
        normalize_gene_symbol(None)


def test_gene_identifier_is_immutable():
    identifier = normalize_gene_symbol("TP53")

    with pytest.raises(Exception):
        identifier.normalized_symbol = "BRCA1"