"""Gene identifier handling for BLC Mark Phase 4 evidence integration."""

from dataclasses import dataclass


@dataclass(frozen=True)
class GeneIdentifier:
    """
    A normalized candidate-gene identifier.

    The original identifier is preserved for traceability while the
    normalized symbol is used internally by Phase 4 when resolvable.
    """

    original_id: str
    normalized_symbol: str | None
    resolvable: bool
    unresolved_reason: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.original_id, str) or not self.original_id.strip():
            raise ValueError("original_id must be a non-empty string.")

        object.__setattr__(
            self,
            "original_id",
            self.original_id.strip(),
        )

        if self.resolvable:
            if (
                not isinstance(self.normalized_symbol, str)
                or not self.normalized_symbol.strip()
            ):
                raise ValueError(
                    "Resolvable identifiers require a non-empty "
                    "normalized_symbol."
                )

            object.__setattr__(
                self,
                "normalized_symbol",
                self.normalized_symbol.strip(),
            )

            if self.unresolved_reason is not None:
                raise ValueError(
                    "Resolvable identifiers cannot have "
                    "an unresolved_reason."
                )

        else:
            if self.normalized_symbol is not None:
                raise ValueError(
                    "Unresolvable identifiers must not have "
                    "a normalized_symbol."
                )

            if (
                not isinstance(self.unresolved_reason, str)
                or not self.unresolved_reason.strip()
            ):
                raise ValueError(
                    "Unresolvable identifiers require "
                    "an unresolved_reason."
                )

            object.__setattr__(
                self,
                "unresolved_reason",
                self.unresolved_reason.strip(),
            )


def normalize_gene_symbol(gene_id: str) -> GeneIdentifier:
    """
    Conservatively normalize a Phase 3 gene identifier.

    Legacy unresolved labels of the form '?|<numeric_id>' are explicitly
    retained as unresolved and are not sent to external evidence sources.

    No alias replacement or case conversion is performed.
    """
    if not isinstance(gene_id, str) or not gene_id.strip():
        raise ValueError("gene_id must be a non-empty string.")

    stripped = gene_id.strip()

    if stripped.startswith("?|"):
        numeric_part = stripped[2:].strip()

        if numeric_part.isdigit():
            return GeneIdentifier(
                original_id=stripped,
                normalized_symbol=None,
                resolvable=False,
                unresolved_reason=(
                    "Legacy unresolved TCGA/Xena identifier."
                ),
            )

    return GeneIdentifier(
        original_id=stripped,
        normalized_symbol=stripped,
        resolvable=True,
        unresolved_reason=None,
    )