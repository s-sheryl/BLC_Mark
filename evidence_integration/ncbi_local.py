"""Local NCBI Gene metadata integration for BLC Mark Phase 4."""

import gzip
from dataclasses import dataclass
from pathlib import Path

from .models import EvidenceRecord, EvidenceType


@dataclass(frozen=True)
class NCBIGeneRecord:
    """Normalized human NCBI Gene metadata."""

    gene_id: str
    symbol: str
    ensembl_id: str | None
    description: str
    aliases: tuple[str, ...]


def _extract_ensembl_id(
    db_xrefs: str,
) -> str | None:
    """
    Extract a unique Ensembl gene identifier from NCBI dbXrefs.

    Returns None when no Ensembl identifier exists or when more than
    one distinct Ensembl gene identifier is present.
    """
    identifiers: list[str] = []

    for item in db_xrefs.split("|"):
        item = item.strip()

        if not item.startswith("Ensembl:"):
            continue

        identifier = item.split(
            ":",
            1,
        )[1].strip()

        if identifier.startswith("ENSG"):
            identifiers.append(identifier)

    unique_identifiers = sorted(
        set(identifiers)
    )

    if len(unique_identifiers) == 1:
        return unique_identifiers[0]

    return None


def load_human_gene_info(
    path: str | Path,
) -> dict[str, NCBIGeneRecord]:
    """
    Load the NCBI Homo sapiens gene_info file.

    The returned dictionary contains only gene symbols that resolve
    uniquely to one NCBI Gene record.

    If the same current symbol occurs for multiple NCBI Gene records,
    that symbol is treated as ambiguous and excluded from the index
    rather than silently selecting one record.
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"NCBI gene_info file does not exist: {path}"
        )

    if not path.is_file():
        raise ValueError(
            f"NCBI gene_info path is not a file: {path}"
        )

    records: dict[str, NCBIGeneRecord] = {}

    ambiguous_symbols: set[str] = set()

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
    ) as handle:
        header = handle.readline()

        if not header.startswith("#tax_id"):
            raise ValueError(
                "NCBI gene_info file has an unexpected header."
            )

        for line_number, raw_line in enumerate(
            handle,
            start=2,
        ):
            line = raw_line.rstrip(
                "\r\n"
            )

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 9:
                raise ValueError(
                    "Invalid NCBI gene_info row at "
                    f"line {line_number}."
                )

            tax_id = fields[0].strip()
            gene_id = fields[1].strip()
            symbol = fields[2].strip()
            synonyms = fields[4].strip()
            db_xrefs = fields[5].strip()
            description = fields[8].strip()

            if tax_id != "9606":
                continue

            if not gene_id:
                raise ValueError(
                    f"Missing GeneID at line {line_number}."
                )

            if not symbol or symbol == "-":
                continue

            aliases = tuple(
                sorted(
                    {
                        alias.strip()
                        for alias in synonyms.split("|")
                        if alias.strip()
                        and alias.strip() != "-"
                    }
                )
            )

            record = NCBIGeneRecord(
                gene_id=gene_id,
                symbol=symbol,
                ensembl_id=_extract_ensembl_id(
                    db_xrefs
                ),
                description=(
                    ""
                    if description == "-"
                    else description
                ),
                aliases=aliases,
            )

            if symbol in ambiguous_symbols:
                continue

            existing = records.get(
                symbol
            )

            if existing is None:
                records[symbol] = record
                continue

            if existing.gene_id == record.gene_id:
                continue

            # Multiple distinct NCBI Gene records use the same
            # current symbol. Preserve scientific ambiguity rather
            # than selecting one record silently.
            records.pop(
                symbol,
                None,
            )

            ambiguous_symbols.add(
                symbol
            )

    return records


def collect_local_functional_evidence(
    gene_symbol: str,
    *,
    ncbi_index: dict[str, NCBIGeneRecord],
    retrieved_at: str,
    source_version: str,
) -> EvidenceRecord | None:
    """
    Create functional evidence from local NCBI gene_info metadata.

    Returns None when the symbol is absent, ambiguous, or lacks a
    description.
    """
    if (
        not isinstance(gene_symbol, str)
        or not gene_symbol.strip()
    ):
        raise ValueError(
            "gene_symbol must be a non-empty string."
        )

    if (
        not isinstance(source_version, str)
        or not source_version.strip()
    ):
        raise ValueError(
            "source_version must be a non-empty string."
        )

    symbol = gene_symbol.strip()

    record = ncbi_index.get(
        symbol
    )

    if record is None:
        return None

    if not record.description:
        return None

    return EvidenceRecord(
        gene_id=symbol,
        evidence_type=EvidenceType.FUNCTIONAL,
        source="NCBI Gene",
        source_version=source_version.strip(),
        evidence_id=(
            f"NCBI_GENE:{record.gene_id}"
        ),
        description=record.description,
        retrieved_at=retrieved_at,
        source_url=(
            "https://www.ncbi.nlm.nih.gov/"
            f"gene/{record.gene_id}"
        ),
    )