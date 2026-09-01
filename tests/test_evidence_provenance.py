import pytest

from src.evidence_integration.provenance import EvidenceSourceMetadata


def test_valid_source_metadata_constructs():
    metadata = EvidenceSourceMetadata(
        source_name="ExampleSource",
        source_version="2026.1",
        source_url="https://example.org",
        retrieved_at="2026-08-26T12:00:00+05:30",
        query_description="Cancer-related evidence for Phase 4 candidates.",
    )

    assert metadata.source_name == "ExampleSource"
    assert metadata.source_version == "2026.1"
    assert metadata.source_url == "https://example.org"


def test_source_version_is_optional():
    metadata = EvidenceSourceMetadata(
        source_name="ExampleSource",
        source_url="https://example.org",
        retrieved_at="2026-08-26T12:00:00+05:30",
        query_description="Evidence retrieval.",
    )

    assert metadata.source_version is None


@pytest.mark.parametrize(
    "field_name",
    [
        "source_name",
        "source_url",
        "retrieved_at",
        "query_description",
    ],
)
def test_required_blank_fields_are_rejected(field_name):
    values = {
        "source_name": "ExampleSource",
        "source_url": "https://example.org",
        "retrieved_at": "2026-08-26T12:00:00+05:30",
        "query_description": "Evidence retrieval.",
    }

    values[field_name] = " "

    with pytest.raises(ValueError, match=field_name):
        EvidenceSourceMetadata(**values)


def test_blank_source_version_is_rejected_when_provided():
    with pytest.raises(ValueError, match="source_version"):
        EvidenceSourceMetadata(
            source_name="ExampleSource",
            source_version=" ",
            source_url="https://example.org",
            retrieved_at="2026-08-26T12:00:00+05:30",
            query_description="Evidence retrieval.",
        )


def test_whitespace_is_removed():
    metadata = EvidenceSourceMetadata(
        source_name=" ExampleSource ",
        source_version=" 2026.1 ",
        source_url=" https://example.org ",
        retrieved_at=" 2026-08-26T12:00:00+05:30 ",
        query_description=" Evidence retrieval. ",
    )

    assert metadata.source_name == "ExampleSource"
    assert metadata.source_version == "2026.1"
    assert metadata.source_url == "https://example.org"
    assert metadata.query_description == "Evidence retrieval."


def test_source_metadata_is_immutable():
    metadata = EvidenceSourceMetadata(
        source_name="ExampleSource",
        source_url="https://example.org",
        retrieved_at="2026-08-26T12:00:00+05:30",
        query_description="Evidence retrieval.",
    )

    with pytest.raises(Exception):
        metadata.source_name = "OtherSource"