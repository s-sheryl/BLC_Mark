"""Local Reactome pathway evidence integration for BLC Mark Phase 4."""

from pathlib import Path

from .models import EvidenceRecord, EvidenceType


EXPECTED_COLUMN_COUNT = 6


def load_reactome_mapping(
    path: str | Path,
) -> dict[str, list[dict[str, str]]]:
    """
    Load an NCBI Gene ID to Reactome pathway mapping file.

    Expected tab-separated columns:
    1. NCBI Gene ID
    2. Reactome stable pathway ID
    3. Reactome pathway URL
    4. Pathway name
    5. Evidence code
    6. Species

    Only Homo sapiens records are retained.
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"Reactome mapping file does not exist: {path}"
        )

    if not path.is_file():
        raise ValueError(
            f"Reactome mapping path is not a file: {path}"
        )

    mapping: dict[str, list[dict[str, str]]] = {}

    with path.open(
        "r",
        encoding="utf-8",
        errors="strict",
    ) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) != EXPECTED_COLUMN_COUNT:
                raise ValueError(
                    "Invalid Reactome mapping row at "
                    f"line {line_number}: expected "
                    f"{EXPECTED_COLUMN_COUNT} columns, "
                    f"found {len(fields)}."
                )

            (
                ncbi_gene_id,
                pathway_id,
                pathway_url,
                pathway_name,
                evidence_code,
                species,
            ) = (field.strip() for field in fields)

            if species != "Homo sapiens":
                continue

            if not ncbi_gene_id:
                raise ValueError(
                    f"Missing NCBI Gene ID at line {line_number}."
                )

            if not pathway_id:
                raise ValueError(
                    f"Missing Reactome pathway ID at line {line_number}."
                )

            if not pathway_name:
                raise ValueError(
                    f"Missing pathway name at line {line_number}."
                )

            mapping.setdefault(
                ncbi_gene_id,
                [],
            ).append(
                {
                    "pathway_id": pathway_id,
                    "pathway_url": pathway_url,
                    "pathway_name": pathway_name,
                    "evidence_code": evidence_code,
                    "species": species,
                }
            )

    for pathways in mapping.values():
        pathways.sort(
            key=lambda item: (
                item["pathway_id"],
                item["pathway_name"],
            )
        )

    return mapping


def collect_pathway_evidence(
    gene_symbol: str,
    ncbi_gene_id: str,
    *,
    reactome_mapping: dict[str, list[dict[str, str]]],
    retrieved_at: str,
    source_version: str,
) -> list[EvidenceRecord]:
    """
    Convert local Reactome pathway mappings into evidence records.

    No pathway ranking or weighting is performed.
    """
    if not isinstance(gene_symbol, str) or not gene_symbol.strip():
        raise ValueError("gene_symbol must be a non-empty string.")

    if not isinstance(ncbi_gene_id, str) or not ncbi_gene_id.strip():
        raise ValueError("ncbi_gene_id must be a non-empty string.")

    if not isinstance(source_version, str) or not source_version.strip():
        raise ValueError("source_version must be a non-empty string.")

    pathways = reactome_mapping.get(
        ncbi_gene_id.strip(),
        [],
    )

    records: list[EvidenceRecord] = []

    for pathway in pathways:
        records.append(
            EvidenceRecord(
                gene_id=gene_symbol.strip(),
                evidence_type=EvidenceType.PATHWAY,
                source="Reactome",
                source_version=source_version.strip(),
                evidence_id=pathway["pathway_id"],
                description=pathway["pathway_name"],
                retrieved_at=retrieved_at,
                source_url=pathway["pathway_url"],
            )
        )

    records.sort(
        key=lambda record: (
            record.evidence_id or "",
            record.description,
        )
    )

    return records