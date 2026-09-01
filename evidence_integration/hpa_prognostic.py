"""Human Protein Atlas prognostic evidence for BLC Mark Phase 4."""

import csv
import io
import zipfile
from pathlib import Path

from .models import EvidenceRecord, EvidenceType


HPA_REQUIRED_COLUMNS = (
    "Gene",
    "Gene name",
    "Cancer",
    "potential prognostic - favorable",
    "unprognostic - favorable",
    "potential prognostic - unfavorable",
    "unprognostic - unfavorable",
    "validated prognostic - favorable",
    "validated prognostic - unfavorable",
)


BLC_MARK_TO_HPA_CANCER = {
    "TCGA-BRCA": "Breast Invasive Carcinoma (TCGA)",
    "TCGA-LUAD": "Lung Adenocarcinoma (TCGA)",
    "TCGA-COAD": "Colon Adenocarcinoma (TCGA)",
}


PROGNOSTIC_COLUMNS = (
    "potential prognostic - favorable",
    "unprognostic - favorable",
    "potential prognostic - unfavorable",
    "unprognostic - unfavorable",
    "validated prognostic - favorable",
    "validated prognostic - unfavorable",
)


def load_hpa_prognostic_data(
    zip_path: str | Path,
) -> dict[tuple[str, str], dict[str, str]]:
    """
    Load HPA cancer prognostic data from its versioned ZIP archive.

    Records are indexed by (Ensembl gene ID, HPA cancer label).
    """
    zip_path = Path(zip_path)

    if not zip_path.exists():
        raise FileNotFoundError(
            f"HPA prognostic archive does not exist: {zip_path}"
        )

    if not zip_path.is_file():
        raise ValueError(
            f"HPA prognostic path is not a file: {zip_path}"
        )

    if not zipfile.is_zipfile(zip_path):
        raise ValueError(
            f"HPA prognostic archive is not a valid ZIP file: {zip_path}"
        )

    with zipfile.ZipFile(zip_path) as archive:
        members = archive.namelist()

        if "cancer_prognostic_data.tsv" not in members:
            raise ValueError(
                "HPA archive does not contain "
                "cancer_prognostic_data.tsv."
            )

        with archive.open("cancer_prognostic_data.tsv") as raw_handle:
            text_handle = io.TextIOWrapper(
                raw_handle,
                encoding="utf-8",
            )

            reader = csv.DictReader(
                text_handle,
                delimiter="\t",
            )

            if reader.fieldnames is None:
                raise ValueError(
                    "HPA prognostic file has no header."
                )

            missing_columns = (
                set(HPA_REQUIRED_COLUMNS)
                - set(reader.fieldnames)
            )

            if missing_columns:
                missing = ", ".join(
                    sorted(missing_columns)
                )

                raise ValueError(
                    "HPA prognostic file is missing "
                    f"required columns: {missing}"
                )

            index: dict[
                tuple[str, str],
                dict[str, str],
            ] = {}

            for row_number, row in enumerate(
                reader,
                start=2,
            ):
                ensembl_id = str(
                    row.get("Gene", "")
                ).strip()

                gene_name = str(
                    row.get("Gene name", "")
                ).strip()

                cancer = str(
                    row.get("Cancer", "")
                ).strip()

                if not ensembl_id:
                    raise ValueError(
                        "Missing Ensembl gene ID at "
                        f"HPA row {row_number}."
                    )

                if not gene_name:
                    raise ValueError(
                        "Missing gene name at "
                        f"HPA row {row_number}."
                    )

                if not cancer:
                    raise ValueError(
                        "Missing cancer label at "
                        f"HPA row {row_number}."
                    )

                key = (
                    ensembl_id,
                    cancer,
                )

                if key in index:
                    raise ValueError(
                        "Duplicate HPA prognostic record for "
                        f"{ensembl_id} / {cancer}."
                    )

                index[key] = {
                    column: str(
                        row.get(column, "")
                    ).strip()
                    for column in HPA_REQUIRED_COLUMNS
                }

    return index


def collect_clinical_evidence(
    gene_symbol: str,
    ensembl_id: str,
    cancer_cohort: str,
    *,
    hpa_index: dict[
        tuple[str, str],
        dict[str, str],
    ],
    retrieved_at: str,
    source_version: str,
) -> list[EvidenceRecord]:
    """
    Convert HPA prognostic annotations into Phase 4 clinical evidence.

    HPA prognostic categories and reported values are preserved exactly.
    BLC Mark does not reinterpret the values or apply a new threshold.
    """
    if not isinstance(gene_symbol, str) or not gene_symbol.strip():
        raise ValueError(
            "gene_symbol must be a non-empty string."
        )

    if not isinstance(ensembl_id, str) or not ensembl_id.strip():
        raise ValueError(
            "ensembl_id must be a non-empty string."
        )

    if cancer_cohort not in BLC_MARK_TO_HPA_CANCER:
        raise ValueError(
            f"Unsupported BLC Mark cancer cohort: {cancer_cohort}"
        )

    if not isinstance(source_version, str) or not source_version.strip():
        raise ValueError(
            "source_version must be a non-empty string."
        )

    symbol = gene_symbol.strip()
    identifier = ensembl_id.strip()

    hpa_cancer = BLC_MARK_TO_HPA_CANCER[
        cancer_cohort
    ]

    row = hpa_index.get(
        (
            identifier,
            hpa_cancer,
        )
    )

    if row is None:
        return []

    hpa_gene_name = str(
        row.get("Gene name", "")
    ).strip()

    if hpa_gene_name and hpa_gene_name != symbol:
        return []

    records: list[EvidenceRecord] = []

    for prognostic_category in PROGNOSTIC_COLUMNS:
        reported_value = str(
            row.get(
                prognostic_category,
                "",
            )
        ).strip()

        if not reported_value:
            continue

        records.append(
            EvidenceRecord(
                gene_id=symbol,
                evidence_type=EvidenceType.CLINICAL,
                source="Human Protein Atlas",
                source_version=source_version.strip(),
                evidence_id=(
                    f"HPA_PROGNOSTIC:"
                    f"{identifier}:"
                    f"{cancer_cohort}:"
                    f"{prognostic_category}"
                ),
                description=(
                    f"{prognostic_category}; "
                    f"reported value={reported_value}; "
                    f"cancer={hpa_cancer}"
                ),
                cancer_cohort=cancer_cohort,
                retrieved_at=retrieved_at,
                source_url=(
                    "https://www.proteinatlas.org/"
                    f"{identifier}-{symbol}/cancer"
                ),
            )
        )

    records.sort(
        key=lambda record: (
            record.evidence_id or "",
            record.description,
        )
    )

    return records