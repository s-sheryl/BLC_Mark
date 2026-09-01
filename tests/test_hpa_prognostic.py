import csv
import zipfile
from pathlib import Path

import pytest

from src.evidence_integration.hpa_prognostic import (
    BLC_MARK_TO_HPA_CANCER,
    HPA_REQUIRED_COLUMNS,
    collect_clinical_evidence,
    load_hpa_prognostic_data,
)
from src.evidence_integration.models import EvidenceType


def write_hpa_zip(
    tmp_path: Path,
    rows: list[dict[str, str]],
) -> Path:
    tsv_path = tmp_path / "cancer_prognostic_data.tsv"

    with tsv_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=HPA_REQUIRED_COLUMNS,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)

    zip_path = tmp_path / "hpa.zip"

    with zipfile.ZipFile(
        zip_path,
        "w",
    ) as archive:
        archive.write(
            tsv_path,
            arcname="cancer_prognostic_data.tsv",
        )

    return zip_path


def make_row(
    *,
    gene: str = "ENSG00000141510",
    gene_name: str = "TP53",
    cancer: str = "Breast Invasive Carcinoma (TCGA)",
) -> dict[str, str]:
    return {
        "Gene": gene,
        "Gene name": gene_name,
        "Cancer": cancer,
        "potential prognostic - favorable": "",
        "unprognostic - favorable": "",
        "potential prognostic - unfavorable": "",
        "unprognostic - unfavorable": "",
        "validated prognostic - favorable": "",
        "validated prognostic - unfavorable": "",
    }


def test_cohort_mapping_is_explicit():
    assert (
        BLC_MARK_TO_HPA_CANCER["TCGA-BRCA"]
        == "Breast Invasive Carcinoma (TCGA)"
    )

    assert (
        BLC_MARK_TO_HPA_CANCER["TCGA-LUAD"]
        == "Lung Adenocarcinoma (TCGA)"
    )

    assert (
        BLC_MARK_TO_HPA_CANCER["TCGA-COAD"]
        == "Colon Adenocarcinoma (TCGA)"
    )


def test_load_hpa_prognostic_data(tmp_path):
    row = make_row()

    path = write_hpa_zip(
        tmp_path,
        [row],
    )

    index = load_hpa_prognostic_data(path)

    key = (
        "ENSG00000141510",
        "Breast Invasive Carcinoma (TCGA)",
    )

    assert key in index
    assert index[key]["Gene name"] == "TP53"


def test_invalid_zip_is_rejected(tmp_path):
    path = tmp_path / "bad.zip"
    path.write_text(
        "not a zip",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="valid ZIP",
    ):
        load_hpa_prognostic_data(path)


def test_missing_file_is_rejected(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_hpa_prognostic_data(
            tmp_path / "missing.zip"
        )


def test_collect_clinical_evidence():
    row = make_row()

    row[
        "potential prognostic - unfavorable"
    ] = "1.250e-03"

    index = {
        (
            "ENSG00000141510",
            "Breast Invasive Carcinoma (TCGA)",
        ): row
    }

    records = collect_clinical_evidence(
        "TP53",
        "ENSG00000141510",
        "TCGA-BRCA",
        hpa_index=index,
        retrieved_at="2026-08-26T13:00:00+05:30",
        source_version="25.1",
    )

    assert len(records) == 1

    record = records[0]

    assert record.evidence_type is EvidenceType.CLINICAL
    assert record.source == "Human Protein Atlas"
    assert record.source_version == "25.1"
    assert record.cancer_cohort == "TCGA-BRCA"

    assert (
        "potential prognostic - unfavorable"
        in record.description
    )

    assert "1.250e-03" in record.description


def test_multiple_hpa_categories_are_preserved():
    row = make_row()

    row[
        "potential prognostic - favorable"
    ] = "2.000e-02"

    row[
        "validated prognostic - favorable"
    ] = "5.000e-03"

    index = {
        (
            "ENSG00000141510",
            "Breast Invasive Carcinoma (TCGA)",
        ): row
    }

    records = collect_clinical_evidence(
        "TP53",
        "ENSG00000141510",
        "TCGA-BRCA",
        hpa_index=index,
        retrieved_at="2026-08-26T13:00:00+05:30",
        source_version="25.1",
    )

    assert len(records) == 2


def test_no_prognostic_annotation_returns_empty():
    row = make_row()

    index = {
        (
            "ENSG00000141510",
            "Breast Invasive Carcinoma (TCGA)",
        ): row
    }

    records = collect_clinical_evidence(
        "TP53",
        "ENSG00000141510",
        "TCGA-BRCA",
        hpa_index=index,
        retrieved_at="2026-08-26T13:00:00+05:30",
        source_version="25.1",
    )

    assert records == []


def test_symbol_mismatch_returns_empty():
    row = make_row(
        gene_name="OTHER",
    )

    index = {
        (
            "ENSG00000141510",
            "Breast Invasive Carcinoma (TCGA)",
        ): row
    }

    records = collect_clinical_evidence(
        "TP53",
        "ENSG00000141510",
        "TCGA-BRCA",
        hpa_index=index,
        retrieved_at="2026-08-26T13:00:00+05:30",
        source_version="25.1",
    )

    assert records == []


def test_unsupported_cohort_is_rejected():
    with pytest.raises(
        ValueError,
        match="Unsupported",
    ):
        collect_clinical_evidence(
            "TP53",
            "ENSG00000141510",
            "TCGA-OV",
            hpa_index={},
            retrieved_at="2026-08-26T13:00:00+05:30",
            source_version="25.1",
        )


@pytest.mark.parametrize(
    "gene_symbol",
    ["", " ", "\t"],
)
def test_blank_gene_symbol_is_rejected(gene_symbol):
    with pytest.raises(
        ValueError,
        match="gene_symbol",
    ):
        collect_clinical_evidence(
            gene_symbol,
            "ENSG00000141510",
            "TCGA-BRCA",
            hpa_index={},
            retrieved_at="2026-08-26T13:00:00+05:30",
            source_version="25.1",
        )