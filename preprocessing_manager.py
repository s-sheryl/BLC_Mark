"""
Purpose:
    Convert structurally validated datasets (files that have already
    passed ValidationManager) into standardized tabular datasets fit
    for downstream analysis -- consistent column names, no duplicate
    or empty rows/columns, numeric expression columns that are
    actually numeric.

Responsibilities:
    - Load a validated tabular file into a pandas DataFrame.
    - Enforce structural expectations: required columns present,
      gene identifier column populated, sample columns identifiable.
    - Clean the table: standardize column names, drop duplicate and
      empty rows/columns, coerce expression columns to numeric,
      handle missing values with a documented, conservative default.
    - Record what changed during preprocessing and report it.

Scope:
    This module prepares files. It does not analyze biology in any
    sense: no differential expression, no statistical testing, no
    pathway analysis, no machine learning, no clinical
    interpretation. Where a real analytical step belongs (log2
    transformation, quantile normalization, batch effect correction,
    gene ID cross-referencing, outlier detection), it is deliberately
    left as a documented stub -- those require scientific judgment
    this module does not make on its own.
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from configuration_manager import ProjectConfig

PREPROCESSING_MANAGER_VERSION = "1.0"

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# The gene identifier column is expected to name genes, not measure
# them -- everything else in a wide expression table is treated as a
# sample column. This default matches the column name convention used
# elsewhere in the project's placeholder datasets; a real dataset may
# need to pass a different name explicitly.
DEFAULT_GENE_ID_COLUMN = "gene_id"

# Placeholder missing-value strategy. Filling with 0.0 is a
# deliberately conservative, easily-reversible default for the
# architecture phase -- it is NOT a scientifically validated
# imputation strategy. Real missing-data handling for expression
# data (e.g. k-NN imputation, dropping low-coverage genes) needs a
# domain decision, not a hardcoded default, and should replace this
# before any real analysis runs on the output.
DEFAULT_MISSING_VALUE_FILL: float | None = None

SUPPORTED_EXTENSIONS = (".txt", ".tsv", ".csv", ".gz")

__all__ = [
    "PREPROCESSING_MANAGER_VERSION",
    "PROJECT_ROOT",
    "DEFAULT_GENE_ID_COLUMN",
    "DEFAULT_MISSING_VALUE_FILL",
    "SUPPORTED_EXTENSIONS",
    "PreprocessingStatus",
    "PreprocessingResult",
    "PreprocessingManager",
]


class PreprocessingStatus(str, Enum):
    """
    Overall outcome of preprocessing one dataset.

    WARNING mirrors the same distinction ValidationManager makes:
    a dataset in WARNING state completed preprocessing and is usable,
    but something about it (e.g. a large fraction of missing values
    that had to be filled) is worth a human's attention before it
    feeds downstream analysis.
    """

    PENDING = "pending"
    PREPROCESSED = "preprocessed"
    WARNING = "warning"
    FAILED = "failed"


@dataclass(frozen=True)
class PreprocessingResult:
    """
    Immutable record of what happened when one dataset was
    preprocessed.

    Frozen for the same reason ValidationResult and DatasetMetadata
    are frozen elsewhere in the project: this is a record of what was
    true at the moment preprocessing ran, not a live object anything
    should mutate afterward. before/after counts are kept explicit
    rather than folded into a single "rows changed" figure, since a
    researcher reviewing this wants to see the shape of the table at
    both ends, not just the delta.
    """

    dataset_id: str
    status: PreprocessingStatus
    passed: bool
    preprocessing_time: datetime

    rows_before: int
    rows_after: int
    columns_before: int
    columns_after: int

    removed_duplicate_rows: int
    removed_empty_rows: int
    removed_empty_columns: int
    missing_values_fixed: int

    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    notes: str | None = None


class PreprocessingManager:
    """
    Prepares structurally validated datasets for downstream analysis.

    Every method here operates on table *structure* -- shape,
    column names, data types, missing values -- and never on what the
    numbers in the table biologically mean. That boundary is what
    keeps this class distinct from the differential expression and
    biomarker discovery stages that come after it in the pipeline.
    """

    def __init__(
    self,
    config: ProjectConfig,
    gene_id_column: str = DEFAULT_GENE_ID_COLUMN,
    missing_value_fill: float | None = DEFAULT_MISSING_VALUE_FILL,
) -> None:
        """
        Store preprocessing configuration. No file or DataFrame is
        touched here -- preprocessing only happens when
        preprocess_dataset() or one of the individual step methods is
        called explicitly.

        Args:
            gene_id_column: Name of the column identifying genes.
                Every other column in the table is treated as a
                sample's expression values.
            missing_value_fill: Placeholder fill value used by
                handle_missing_values(). See DEFAULT_MISSING_VALUE_FILL
                for why this is a conservative default rather than a
                scientific choice.
        """
        self.config = config
        self.processed_dir = config.processed_dir
        if self.processed_dir is None:
         raise ValueError(
        "ProjectConfig.processed_dir must be configured before "
        "preprocessing datasets."
    )
        self.gene_id_column = gene_id_column
        self.missing_value_fill = missing_value_fill

    def load_dataset(self, file_path: Path) -> pd.DataFrame:
        """
        Load a validated tabular file into a DataFrame.

        Delimiter is inferred from the file's extension rather than
        guessed from content: `.csv` is comma-delimited, `.tsv`/`.txt`
        are tab-delimited, and `.gz` is treated as a compressed
        version of whichever delimiter its inner suffix implies
        (e.g. `expression.tsv.gz`). pandas handles gzip decompression
        automatically via `compression="infer"`.

        Args:
            file_path: Path to a file that has already passed
                ValidationManager checks.

        Returns:
            The loaded table as a DataFrame, exactly as read -- no
            cleaning has happened yet at this point.

        Raises:
            ValueError: If the file's extension isn't one this
                method knows how to delimit.
            pandas.errors.EmptyDataError: If the file has no parseable
                content.
        """
        suffixes = file_path.suffixes
        inner_suffix = suffixes[-2] if file_path.suffix == ".gz" and len(suffixes) > 1 else file_path.suffix

        if inner_suffix == ".csv":
            separator = ","
        elif inner_suffix in (".tsv", ".txt"):
            separator = "\t"
        else:
            raise ValueError(
                f"Cannot infer delimiter for {file_path} "
                f"(recognized extensions: {SUPPORTED_EXTENSIONS})."
            )

        return pd.read_csv(file_path, sep=separator, compression="infer")

    def validate_required_columns(
        self, df: pd.DataFrame, required_columns: list[str]
    ) -> bool:
        """
        Confirm every required column is present in the table.

        Args:
            df: The table to check.
            required_columns: Column names that must be present.

        Returns:
            True if every required column exists.

        Raises:
            ValueError: If one or more required columns are missing.
        """
        missing = [column for column in required_columns if column not in df.columns]
        if missing:
            raise ValueError(f"Missing required column(s): {missing}")
        return True

    def standardize_column_names(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalize column names to a consistent convention: stripped
        whitespace, lowercase, spaces replaced with underscores.

        Datasets pulled from different public sources rarely agree on
        column naming conventions (" Gene ID" vs "gene_id" vs
        "GeneID"). Standardizing here means every later step in the
        pipeline can rely on one naming convention instead of each
        one re-normalizing independently.

        Args:
            df: The table whose columns should be standardized.

        Returns:
            A new DataFrame with standardized column names. The
            input DataFrame is not modified in place.
        """
        renamed = {
            column: str(column).strip().lower().replace(" ", "_")
            for column in df.columns
        }
        return df.rename(columns=renamed)

    def remove_duplicate_rows(self, df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
        """
        Drop rows that are exact duplicates of another row.

        Args:
            df: The table to deduplicate.

        Returns:
            A tuple of (deduplicated DataFrame, number of rows removed).
        """
        rows_before = len(df)
        deduplicated = df.drop_duplicates()
        removed = rows_before - len(deduplicated)
        return deduplicated, removed

    def remove_empty_rows(self, df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
        """
        Drop rows where every value is missing.

        A row with no values at all carries no information and would
        only distort later missing-value handling and statistics.

        Args:
            df: The table to clean.

        Returns:
            A tuple of (cleaned DataFrame, number of rows removed).
        """
        rows_before = len(df)
        cleaned = df.dropna(axis="index", how="all")
        removed = rows_before - len(cleaned)
        return cleaned, removed

    def remove_empty_columns(self, df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
        """
        Drop columns where every value is missing.

        Args:
            df: The table to clean.

        Returns:
            A tuple of (cleaned DataFrame, number of columns removed).
        """
        columns_before = df.shape[1]
        cleaned = df.dropna(axis="columns", how="all")
        removed = columns_before - cleaned.shape[1]
        return cleaned, removed

    def validate_gene_identifier_column(self, df: pd.DataFrame) -> bool:
        """
        Confirm the configured gene identifier column exists and is
        populated.

        Args:
            df: The table to check.

        Returns:
            True if the column exists and has at least one non-null
            value.

        Raises:
            ValueError: If the column is missing entirely, or every
                value in it is null.
        """
        if self.gene_id_column not in df.columns:
            raise ValueError(
                f"Gene identifier column '{self.gene_id_column}' not found. "
                f"Available columns: {list(df.columns)}"
            )

        if df[self.gene_id_column].isna().all():
            raise ValueError(
                f"Gene identifier column '{self.gene_id_column}' is "
                "entirely empty."
            )

        return True

    def validate_sample_identifier_columns(self, df: pd.DataFrame) -> list[str]:
        """
        Identify and confirm the presence of sample columns.

        Every column other than the gene identifier column is treated
        as a sample's expression values, consistent with the wide
        format TCGA expression matrices typically use (one row per
        gene, one column per sample).

        Args:
            df: The table to check.

        Returns:
            The list of column names identified as sample columns.

        Raises:
            ValueError: If no sample columns remain besides the gene
                identifier column.
        """
        sample_columns = [column for column in df.columns if column != self.gene_id_column]

        if not sample_columns:
            raise ValueError(
                "No sample columns found -- table only contains the "
                f"gene identifier column '{self.gene_id_column}'."
            )

        return sample_columns

    def validate_numeric_columns(
        self, df: pd.DataFrame, sample_columns: list[str]
    ) -> tuple[pd.DataFrame, list[str]]:
        """
        Coerce sample columns to numeric dtype and flag problems.

        Values that cannot be parsed as numbers become NaN (handled
        later by handle_missing_values()) rather than silently
        dropped rows, so a single malformed cell doesn't destroy an
        otherwise good row. Infinite values are flagged explicitly,
        since they usually indicate a division-by-zero or parsing
        error upstream rather than a real biological measurement.

        Args:
            df: The table containing the sample columns.
            sample_columns: Columns expected to hold numeric
                expression values.

        Returns:
            A tuple of (DataFrame with coerced numeric columns, list
            of warning messages describing any coercion problems or
            infinite values found).
        """
        coerced = df.copy()
        warnings: list[str] = []

        for column in sample_columns:
            original_non_null = coerced[column].notna().sum()
            coerced[column] = pd.to_numeric(coerced[column], errors="coerce")
            new_non_null = coerced[column].notna().sum()

            newly_invalid = original_non_null - new_non_null
            if newly_invalid > 0:
                warnings.append(
                    f"Column '{column}': {newly_invalid} value(s) could not "
                    "be parsed as numeric and were treated as missing."
                )

            infinite_count = int(np.isinf(coerced[column]).sum())
            if infinite_count > 0:
                warnings.append(
                    f"Column '{column}': {infinite_count} infinite value(s) "
                    "detected."
                )
                coerced[column] = coerced[column].replace(
                    [np.inf, -np.inf], np.nan
                )

        return coerced, warnings

    
    def handle_missing_values(
        self, df: pd.DataFrame, sample_columns: list[str]
    ) -> tuple[pd.DataFrame, int]:
        """
        Detect missing values in sample columns without applying
        scientifically unjustified imputation.

        Missing expression values are preserved as NaN. The method only
        records how many missing values are present so that a scientifically
        justified downstream analysis can decide how they should be handled.

        Args:
            df: The table to inspect.
            sample_columns: Columns containing expression values.

        Returns:
            A tuple of:
                - unchanged DataFrame with missing values preserved as NaN;
                - number of missing values detected.
        """
        checked = df.copy()

        missing_values = int(
            checked[sample_columns].isna().sum().sum()
        )

        return checked, missing_values

    def standardize_gene_identifiers(
        self, df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Apply light, purely textual standardization to gene
        identifiers: strip whitespace, uppercase.

        This is placeholder architecture only. It does NOT cross-
        reference identifiers between naming systems (e.g. Ensembl
        gene IDs vs. HGNC symbols vs. Entrez IDs) -- that requires an
        authoritative mapping table and is explicitly out of scope
        here.

        Args:
            df: The table containing the gene identifier column.

        Returns:
            A new DataFrame with gene identifiers stripped and
            uppercased. The input DataFrame is not modified in place.
        """
        standardized = df.copy()

        standardized[self.gene_id_column] = (
            standardized[self.gene_id_column]
            .astype(str)
            .str.strip()
            .str.upper()
        )

        return standardized   

    
        

    def collect_preprocessing_statistics(
        self,
        *,
        rows_before: int,
        rows_after: int,
        columns_before: int,
        columns_after: int,
        removed_duplicate_rows: int,
        removed_empty_rows: int,
        removed_empty_columns: int,
        missing_values_detected: int,
    ) -> dict[str, int]:
        """
        Package raw preprocessing counters into a single statistics
        dictionary.

        Kept as its own method, separate from preprocess_dataset(),
        so the shape of "what counts as a preprocessing statistic" is
        defined in exactly one place and can be reused (e.g. by
        build_preprocessing_report()) without recomputation.

        Args:
            rows_before: Row count before any cleaning.
            rows_after: Row count after all cleaning steps.
            columns_before: Column count before any cleaning.
            columns_after: Column count after all cleaning steps.
            removed_duplicate_rows: Count from remove_duplicate_rows().
            removed_empty_rows: Count from remove_empty_rows().
            removed_empty_columns: Count from remove_empty_columns().
            missing_values_detected: Count from handle_missing_values().

        Returns:
            A dictionary of the same values, keyed by name -- mainly
            useful for logging or serialization contexts that want a
            dict rather than positional fields.
        """
        return {
            "rows_before": rows_before,
            "rows_after": rows_after,
            "columns_before": columns_before,
            "columns_after": columns_after,
            "removed_duplicate_rows": removed_duplicate_rows,
            "removed_empty_rows": removed_empty_rows,
            "removed_empty_columns": removed_empty_columns,
            "missing_values_detected": missing_values_detected,
        }
    
    def preprocess_dataset(
        self,
        file_path: Path,
        dataset_id: str,
    ) -> PreprocessingResult:
        """
        Run the full preprocessing pipeline against one dataset file.

        Steps run in a fixed order: load, standardize column names,
        validate structure (required columns, gene ID column, sample
        columns), remove duplicate/empty rows and columns, coerce
        sample columns to numeric, fill missing values, and lightly
        standardize gene identifiers. Any step that raises is caught
        and recorded as an error, and preprocessing stops at that
        point -- unlike ValidationManager's checks, these steps
        depend on each other's output (e.g. numeric coercion needs a
        confirmed sample column list), so continuing after a
        structural failure here would only produce confusing
        downstream errors.

        Args:
            file_path: Path to a dataset file that has already
                passed ValidationManager.
            dataset_id: Identifier for this dataset, carried through
                to the returned result.

        Returns:
            A PreprocessingResult describing what happened. Never
            raises for expected failure modes -- pipeline failures
            are captured in the result's status/errors instead.
        """
        warnings: list[str] = []
        errors: list[str] = []
        removed_duplicate_rows = 0
        removed_empty_rows = 0
        removed_empty_columns = 0
        missing_values_fixed = 0
        rows_before = 0
        columns_before = 0
        rows_after = 0
        columns_after = 0

        try:
            df = self.load_dataset(file_path)
            rows_before, columns_before = df.shape

            df = self.standardize_column_names(df)
            self.validate_gene_identifier_column(df)
            sample_columns = self.validate_sample_identifier_columns(df)

            df, removed_duplicate_rows = self.remove_duplicate_rows(df)
            df, removed_empty_rows = self.remove_empty_rows(df)
            df, removed_empty_columns = self.remove_empty_columns(df)

            # Sample columns may have shifted if remove_empty_columns
            # dropped one -- recompute rather than trust the earlier list.
            sample_columns = [c for c in sample_columns if c in df.columns]
            if not sample_columns:raise ValueError(
        "No sample columns remain after removing empty columns."
    )
            df, numeric_warnings = self.validate_numeric_columns(df, sample_columns)
            warnings.extend(numeric_warnings)

            df, missing_values_detected = self.handle_missing_values(
    df, sample_columns
)

            rows_after, columns_after = df.shape

        except (ValueError, FileNotFoundError, pd.errors.EmptyDataError) as error:
            errors.append(str(error))

        if errors:
            status = PreprocessingStatus.FAILED
        elif warnings:
            status = PreprocessingStatus.WARNING
        else:
            status = PreprocessingStatus.PREPROCESSED

        return PreprocessingResult(
            dataset_id=dataset_id,
            status=status,
            passed=not errors,
            preprocessing_time=datetime.now(timezone.utc),
            rows_before=rows_before,
            rows_after=rows_after,
            columns_before=columns_before,
            columns_after=columns_after,
            removed_duplicate_rows=removed_duplicate_rows,
            removed_empty_rows=removed_empty_rows,
            removed_empty_columns=removed_empty_columns,
            missing_values_detected=missing_values_detected,
            warnings=warnings,
            errors=errors,
        )

    def build_preprocessing_report(
        self, results: list[PreprocessingResult]
    ) -> dict[str, Any]:
        """
        Summarize a batch of PreprocessingResult objects.

        Pure aggregation over results already computed -- no file
        access happens here.

        Args:
            results: Preprocessing results to summarize.

        Returns:
            A dictionary with:
                - "total": count of results
                - "preprocessed": count with status == PREPROCESSED
                - "warnings": count with status == WARNING
                - "failed": count with status == FAILED
                - "datasets": list of {"dataset_id", "status"} entries
        """
        return {
            "total": len(results),
            "preprocessed": sum(
                1 for r in results if r.status == PreprocessingStatus.PREPROCESSED
            ),
            "warnings": sum(
                1 for r in results if r.status == PreprocessingStatus.WARNING
            ),
            "failed": sum(1 for r in results if r.status == PreprocessingStatus.FAILED),
            "datasets": [
                {"dataset_id": r.dataset_id, "status": r.status.value} for r in results
            ],
        }

    def render_preprocessing_summary(self, report: dict[str, Any]) -> str:
        """
        Render a preprocessing report as a human-readable, multi-line
        string suitable for console output.

        Args:
            report: The dictionary produced by build_preprocessing_report().

        Returns:
            A formatted multi-line summary string.
        """
        lines = [
            "Preprocessing Summary",
            "----------------------",
            f"Total datasets:  {report['total']}",
            f"Preprocessed:    {report['preprocessed']}",
            f"Warnings:        {report['warnings']}",
            f"Failed:          {report['failed']}",
            "",
            "Per-dataset status:",
        ]

        for dataset in report["datasets"]:
            lines.append(f"  [{dataset['status'].upper():>12}] {dataset['dataset_id']}")

        return "\n".join(lines)

    # -- Not implemented yet: real analytical/statistical steps --

    def export_report_json(self, report: dict[str, Any], destination: Path) -> None:
        """
        Export a preprocessing report to a JSON file.

        TODO(implementation): serialize `report` and write it to
        `destination`.

        Raises:
            NotImplementedError: Always, until implemented.
        """
        raise NotImplementedError(
            "PreprocessingManager.export_report_json() is not "
            f"implemented yet. It will write the report to {destination}."
        )

    def export_report_csv(
        self, results: list[PreprocessingResult], destination: Path
    ) -> None:
        """
        Export per-dataset preprocessing results to a CSV file.

        TODO(implementation): flatten each PreprocessingResult to a
        row and write to `destination` using the standard library
        csv module.

        Raises:
            NotImplementedError: Always, until implemented.
        """
        raise NotImplementedError(
            "PreprocessingManager.export_report_csv() is not "
            f"implemented yet. It will write {len(results)} result(s) "
            f"to {destination}."
        )

    def normalize_expression(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalize expression values across samples (e.g. TMM, RLE, or
        another established normalization method).

        TODO(implementation): requires a scientifically justified
        choice of normalization method appropriate to the assay type
        (RNA-seq counts vs. already-normalized RSEM values, etc.) --
        not something to default silently.

        Raises:
            NotImplementedError: Always, until implemented.
        """
        raise NotImplementedError(
            "PreprocessingManager.normalize_expression() is not "
            "implemented yet. Requires a chosen, justified "
            "normalization method before implementation."
        )

    def batch_effect_correction(self, df: pd.DataFrame, batch_labels: pd.Series) -> pd.DataFrame:
        """
        Correct for batch effects across samples (e.g. ComBat).

        TODO(implementation): requires known batch labels per sample
        (e.g. sequencing center, processing date) which are not part
        of the raw expression matrix and must be sourced from
        clinical/sample metadata first.

        Raises:
            NotImplementedError: Always, until implemented.
        """
        raise NotImplementedError(
            "PreprocessingManager.batch_effect_correction() is not "
            "implemented yet. Requires verified batch labels per "
            "sample before implementation."
        )

    def gene_id_mapping(self, df: pd.DataFrame, target_id_system: str) -> pd.DataFrame:
        """
        Map gene identifiers between naming systems (e.g. Ensembl ID
        to HGNC symbol).

        TODO(implementation): requires an authoritative, versioned
        mapping table (e.g. from Ensembl BioMart or HGNC) rather than
        a hand-built dictionary, so mappings stay traceable to a
        specific reference release.

        Raises:
            NotImplementedError: Always, until implemented.
        """
        raise NotImplementedError(
            "PreprocessingManager.gene_id_mapping() is not implemented "
            f"yet. Requires a versioned mapping table to "
            f"'{target_id_system}' before implementation."
        )

    def log2_transformation(self, df: pd.DataFrame, sample_columns: list[str]) -> pd.DataFrame:
        """
        Apply a log2(x + 1) transformation to expression values.

        TODO(implementation): requires confirming whether the input
        values are already log-transformed (RSEM outputs sometimes
        are) before applying this -- doing it twice silently distorts
        the data.

        Raises:
            NotImplementedError: Always, until implemented.
        """
        raise NotImplementedError(
            "PreprocessingManager.log2_transformation() is not "
            "implemented yet. Requires confirming input values are "
            "not already log-transformed before implementation."
        )

    def quantile_normalization(self, df: pd.DataFrame, sample_columns: list[str]) -> pd.DataFrame:
        """
        Apply quantile normalization across sample columns.

        TODO(implementation): a real implementation choice (e.g.
        via a maintained library) rather than a hand-rolled quantile
        algorithm, to avoid subtle correctness bugs in a step that
        directly affects every downstream statistical result.

        Raises:
            NotImplementedError: Always, until implemented.
        """
        raise NotImplementedError(
            "PreprocessingManager.quantile_normalization() is not "
            "implemented yet."
        )

    def outlier_detection(self, df: pd.DataFrame, sample_columns: list[str]) -> list[str]:
        """
        Identify sample columns that are statistical outliers
        relative to the rest of the cohort.

        TODO(implementation): requires a chosen outlier criterion
        (e.g. PCA-based distance, correlation-based) validated against
        this project's actual data before being trusted to flag or
        exclude real samples.

        Raises:
            NotImplementedError: Always, until implemented.
        """
        raise NotImplementedError(
            "PreprocessingManager.outlier_detection() is not "
            "implemented yet. Requires a validated outlier criterion "
            "before implementation."
        )


if __name__ == "__main__":
    print("BLC Mark Preprocessing Manager")
    print(f"Version: {PREPROCESSING_MANAGER_VERSION}\n")

    # Synthetic, non-biological demonstration data: a small numeric
    # table with the shape of an expression matrix (gene identifier
    # column + sample columns), deliberately including a duplicate
    # row, an empty row, a missing value, and a non-numeric cell so
    # the pipeline has something real to clean.
    synthetic_data = pd.DataFrame(
        {
            "Gene ID": ["gene_a", "gene_b", "gene_c", "gene_a", None],
            "Sample 1": [1.0, 2.5, "not_a_number", 1.0, None],
            "Sample 2": [3.1, None, 4.4, 3.1, None],
        }
    )

    manager = PreprocessingManager(gene_id_column="gene_id")

    import tempfile

    with tempfile.TemporaryDirectory() as scratch_dir:
        scratch_file = Path(scratch_dir) / "synthetic_expression.csv"
        synthetic_data.to_csv(scratch_file, index=False)

        result = manager.preprocess_dataset(
            file_path=scratch_file, dataset_id="SYNTHETIC.PREPROCESS.TEST"
        )

    print(f"Dataset: {result.dataset_id}")
    print(f"Status: {result.status.value.upper()}")
    print(f"Rows before -> after: {result.rows_before} -> {result.rows_after}")
    print(f"Columns before -> after: {result.columns_before} -> {result.columns_after}")
    print(f"Duplicate rows removed: {result.removed_duplicate_rows}")
    print(f"Empty rows removed: {result.removed_empty_rows}")
    print(f"Empty columns removed: {result.removed_empty_columns}")
    print(f"Missing values detected: {result.missing_values_detected}")
    if result.warnings:
        print(f"Warnings: {result.warnings}")
    if result.errors:
        print(f"Errors: {result.errors}")

    report = manager.build_preprocessing_report([result])
    print()
    print(manager.render_preprocessing_summary(report))