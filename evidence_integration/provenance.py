"""Provenance metadata for BLC Mark Phase 4 evidence integration."""

from dataclasses import dataclass


@dataclass(frozen=True)
class EvidenceSourceMetadata:
    """
    Provenance information for a Phase 4 evidence source.

    Parameters
    ----------
    source_name:
        Human-readable source name.

    source_version:
        Source release or version when available.

    source_url:
        Canonical URL for the evidence source.

    retrieved_at:
        ISO 8601 retrieval timestamp.

    query_description:
        Human-readable description of what was retrieved.
    """

    source_name: str
    source_url: str
    retrieved_at: str
    query_description: str
    source_version: str | None = None

    def __post_init__(self) -> None:
        required_fields = {
            "source_name": self.source_name,
            "source_url": self.source_url,
            "retrieved_at": self.retrieved_at,
            "query_description": self.query_description,
        }

        for field_name, value in required_fields.items():
            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    f"{field_name} must be a non-empty string."
                )

        object.__setattr__(
            self,
            "source_name",
            self.source_name.strip(),
        )
        object.__setattr__(
            self,
            "source_url",
            self.source_url.strip(),
        )
        object.__setattr__(
            self,
            "retrieved_at",
            self.retrieved_at.strip(),
        )
        object.__setattr__(
            self,
            "query_description",
            self.query_description.strip(),
        )

        if self.source_version is not None:
            if (
                not isinstance(self.source_version, str)
                or not self.source_version.strip()
            ):
                raise ValueError(
                    "source_version must be a non-empty string when provided."
                )

            object.__setattr__(
                self,
                "source_version",
                self.source_version.strip(),
            )