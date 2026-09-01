"""
Purpose:
    Define the typed, frozen data structures shared across the Phase 3
    (Differential Expression) package: analysis configuration,
    comparison/group definitions, validated inputs, gene-level
    statistical results, QC reporting, reproducibility metadata, and
    the final analysis result.

Responsibilities:
    - Represent every scientifically relevant DE parameter as an
      explicit, typed field rather than leaving it implicit in code.
    - Stay pure data: no file I/O, no statistical computation, no
      validation logic beyond structural/type checks performed in
      __post_init__.

Scope:
    This module deliberately does not decide values for any field --
    it only defines what shape a value must have once decided
    elsewhere (validation.py, configuration.py, comparison.py,
    filtering.py, methods.py, multiple_testing.py, qc.py,
    reproducibility.py, results.py).

    Every field included here is traceable to a specific requirement
    in docs/differential_expression_specification.md. No field is
    included "for completeness" without that basis.

Version:
    Phase 3 Version 1.0 DE models.
"""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path

import pandas as pd

from src.dataset_registry import CancerType

DE_MODELS_VERSION = "1.0"

# Minimum per-group sample count below which no two-sample statistical
# comparison can even be computed (a sample variance requires at least
# two observations). This is a mathematical floor, not a scientifically
# meaningful minimum -- specification Section 4.3 explicitly defers the
# scientifically meaningful minimum to the statistical method chosen at
# implementation time, and requires it to be recorded rather than
# silently assumed. DEAnalysisConfiguration.minimum_replicates_per_group
# must therefore be supplied explicitly by the caller and is validated
# to be no lower than this floor.
ABSOLUTE_MINIMUM_REPLICATES_PER_GROUP = 2

DEFAULT_SIGNIFICANCE_THRESHOLD = 0.05

__all__ = [
    "DE_MODELS_VERSION",
    "ABSOLUTE_MINIMUM_REPLICATES_PER_GROUP",
    "DEFAULT_SIGNIFICANCE_THRESHOLD",
    "ExpressionRepresentation",
    "StatisticalMethod",
    "MultipleTestingMethod",
    "AnalysisStatus",
    "GeneFilterConfiguration",
    "DEAnalysisConfiguration",
    "ComparisonDefinition",
    "ExcludedSample",
    "GroupAssignment",
    "ValidatedExpressionMatrix",
    "ValidatedMetadata",
    "SampleMatchResult",
    "GeneFilterResult",
    "GeneResult",
    "QCReport",
    "AnalysisMetadata",
    "FailureInfo",
    "AnalysisResult",
]


class ExpressionRepresentation(str, Enum):
    """The numeric representation of values in an expression matrix.

    Differential expression specification Section 5.2/8.1 requires
    the pipeline to know this before selecting a statistical method,
    and explicitly forbids inferring it from filenames, magnitudes, or
    column names -- it must be supplied explicitly as configuration.

    RAW_COUNTS:
        Non-negative integer (or effectively integer) RNA-seq read
        counts, suitable in principle for a count-based method such
        as DESeq2 (specification Section 5.1).

    NORMALIZED_LOG2:
        Already-normalized expression values that have also been
        log2-transformed (e.g. RSEM values as distributed by some
        public TCGA resources; see src/dataset_registry.py platform
        notes). A log2 fold-change effect size is meaningful for this
        representation.

    NORMALIZED_LINEAR:
        Already-normalized expression values on a linear (non-log)
        scale (e.g. linear TPM/FPKM/RSEM). A "log2 fold change" label
        is not meaningful for this representation without an explicit
        log transform, which this package does not silently apply
        (specification Section 7 of the safety rules: never
        log-transform data as a side effect of analysis).
    """

    RAW_COUNTS = "raw_counts"
    NORMALIZED_LOG2 = "normalized_log2"
    NORMALIZED_LINEAR = "normalized_linear"


class StatisticalMethod(str, Enum):
    """Statistical methods a DE analysis may be configured to use.

    DESEQ2:
        Named explicitly by specification Section 5.1 as the
        preferred method for raw-count expression data. No
        DESeq2-capable backend (R/Bioconductor, or a Python
        equivalent such as pydeseq2) is installed in this
        environment. Configuring DESEQ2 therefore always fails
        explicitly at execution time with UnsupportedMethodError --
        see methods.py for the exact failure behavior. This is a
        genuine environment limitation, not a simulated one: DESeq2
        is never faked or silently replaced.

    WELCH_T_TEST:
        Specification Section 5.2 requires, for already-normalized
        expression data, "a statistically appropriate method for
        that representation" to be used, explicitly justified and
        documented (framework, model design, normalization
        assumptions, and the statistical test applied, per Section
        5.1's documentation requirement applied to whichever method
        is chosen). Welch's two-sample t-test (unequal-variance,
        two-sided) is the statistically appropriate, well-established
        method for comparing two independent groups' means without
        assuming raw-count structure or equal variances, and is
        implemented here as that justified method for the
        NORMALIZED_LOG2 and NORMALIZED_LINEAR expression
        representations only. It is not authorized, and methods.py
        explicitly rejects it, for RAW_COUNTS data, since a t-test is
        not the statistically appropriate method for count data with
        its own dispersion structure -- that is DESeq2's domain per
        Section 5.1, and DESeq2 is unavailable rather than
        substituted.
    """

    DESEQ2 = "deseq2"
    WELCH_T_TEST = "welch_t_test"


class MultipleTestingMethod(str, Enum):
    """Multiple-testing correction methods a DE analysis may use.

    Version 1 supports only Benjamini-Hochberg FDR control, per
    specification Section 6.1. No other correction method is
    authorized.
    """

    BENJAMINI_HOCHBERG = "benjamini_hochberg"


class AnalysisStatus(str, Enum):
    """Lifecycle status of a differential expression analysis."""

    PENDING = "pending"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"


@dataclass(frozen=True)
class GeneFilterConfiguration:
    """Explicit, caller-supplied gene-filtering configuration.

    Specification Section 8.2 requires any gene filtering to be
    "explicitly defined, reproducible, recorded ... and appropriate
    for the selected statistical method," but the specification does
    not itself define a filtering criterion or threshold. Per the
    project's implementation rules, no default biological filtering
    threshold is invented here.

    The Version 1 default is therefore no filtering at all
    (apply_filter=False). A caller may explicitly opt into a single,
    simple, fully documented filtering rule -- excluding genes whose
    mean expression across all included samples falls below an
    explicit threshold -- by setting apply_filter=True and providing
    minimum_mean_expression. No other filtering criterion is
    implemented in Version 1; if a different criterion is needed, it
    must be added as an explicit specification change, not invented
    here.

    Attributes:
        apply_filter:
            Whether gene filtering is applied at all. Defaults to
            False (no filtering), which is the only behavior directly
            supported by the specification without further scientific
            input.

        minimum_mean_expression:
            Required when apply_filter is True. A gene is filtered
            out (excluded from statistical testing) if its mean
            expression across all samples included in the comparison
            is strictly below this value. Must be None when
            apply_filter is False.

        criterion_description:
            Human-readable description of the filtering criterion
            actually applied, recorded verbatim in QC and
            reproducibility metadata for traceability.
    """

    apply_filter: bool = False
    minimum_mean_expression: float | None = None
    criterion_description: str = (
        "No gene filtering configured for this analysis "
        "(Version 1 default: no filter applied)."
    )

    def __post_init__(self) -> None:
        if not isinstance(self.apply_filter, bool):
            raise TypeError(
                "'apply_filter' must be a bool, "
                f"got {type(self.apply_filter).__name__}."
            )

        if self.apply_filter:
            if self.minimum_mean_expression is None:
                raise ValueError(
                    "'minimum_mean_expression' must be provided when "
                    "'apply_filter' is True: Version 1 does not invent "
                    "a default filtering threshold."
                )
            if not isinstance(self.minimum_mean_expression, (int, float)) or isinstance(
                self.minimum_mean_expression, bool
            ):
                raise TypeError(
                    "'minimum_mean_expression' must be an int or float, "
                    f"got {type(self.minimum_mean_expression).__name__}."
                )
        else:
            if self.minimum_mean_expression is not None:
                raise ValueError(
                    "'minimum_mean_expression' must be None when "
                    "'apply_filter' is False."
                )


@dataclass(frozen=True)
class DEAnalysisConfiguration:
    """Complete, explicit configuration for one differential
    expression analysis.

    Every scientifically relevant setting the specification requires
    to be explicit (Sections 3.4, 4.1, 4.3, 5.2, 5.3, 6.2, 6.3, 8.2,
    9.1) is a named field here with no silently-inferred value.

    Attributes:
        analysis_id:
            Caller-supplied unique identifier for this analysis
            (specification Section 9.5). Not generated implicitly,
            so the same identifier can be supplied again by a caller
            reproducing a prior analysis.

        cancer_cohort:
            One of the four Version 1 in-scope cancer types, reusing
            the existing src.dataset_registry.CancerType enum rather
            than duplicating that list.

        expression_matrix_path:
            Path to the processed expression matrix.

        metadata_path:
            Path to the sample metadata file.

        gene_id_column:
            Name of the gene-identifier column in the expression
            matrix.

        sample_id_column:
            Name of the sample-identifier column in the metadata
            file.

        group_column:
            Name of the group/condition column in the metadata file.

        reference_group:
            The explicitly configured reference group label
            (specification Section 4.1/4.2). Never inferred.

        comparison_group:
            The explicitly configured comparison group label. Must
            differ from reference_group.

        expression_representation:
            The declared numeric representation of the expression
            matrix. Required; never inferred (specification Section
            5.2 and the project's expression-representation rule).

        statistical_method:
            The configured statistical method. See StatisticalMethod
            for Version 1's single authorized value.

        multiple_testing_method:
            The configured multiple-testing correction method.
            Defaults to, and for Version 1 can only be,
            Benjamini-Hochberg (specification Section 6.1).

        significance_threshold:
            The adjusted-p-value (FDR) threshold used to annotate
            significance in the result table (specification Section
            6.2). Defaults to 0.05, the specification's stated
            default, but is always recorded explicitly and never
            silently changed during execution.

        minimum_replicates_per_group:
            The caller-supplied minimum number of samples required in
            each comparison group before the analysis may proceed.
            Required, with no default: specification Section 4.3
            explicitly declines to state a scientifically meaningful
            minimum and defers it to the statistical method chosen at
            implementation time, so this package requires the caller
            to state the value explicitly rather than inventing one.
            Must be at least ABSOLUTE_MINIMUM_REPLICATES_PER_GROUP.

        effect_size_threshold:
            Optional effect-size threshold. Per specification Section
            6.3, "an effect-size threshold must not be applied unless
            it has been explicitly configured" -- so this defaults to
            None (no threshold applied) and is only used if the
            caller explicitly sets it.

        gene_filter:
            Explicit gene-filtering configuration. See
            GeneFilterConfiguration.

        output_dir:
            Directory into which DE outputs (results, metadata, QC
            report) will be written.

        random_seed:
            Optional random seed. Recorded for reproducibility
            (specification Section 9.4) if any configured operation
            involves randomness. Both of Version 1's statistical
            methods (Welch's t-test via scipy.stats.ttest_ind, and
            DESeq2 when/if available) and Benjamini-Hochberg
            correction are deterministic, so this is expected to
            remain None in Version 1, but the field exists so a seed
            is never silently omitted if a future method requires
            one.
    """

    analysis_id: str
    cancer_cohort: CancerType
    expression_matrix_path: Path
    metadata_path: Path
    gene_id_column: str
    sample_id_column: str
    group_column: str
    reference_group: str
    comparison_group: str
    expression_representation: ExpressionRepresentation
    statistical_method: StatisticalMethod
    minimum_replicates_per_group: int
    output_dir: Path
    multiple_testing_method: MultipleTestingMethod = (
        MultipleTestingMethod.BENJAMINI_HOCHBERG
    )
    significance_threshold: float = DEFAULT_SIGNIFICANCE_THRESHOLD
    effect_size_threshold: float | None = None
    gene_filter: GeneFilterConfiguration = field(
        default_factory=GeneFilterConfiguration
    )
    random_seed: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.analysis_id, str) or not self.analysis_id.strip():
            raise ValueError("'analysis_id' must be a non-empty str.")

        if not isinstance(self.cancer_cohort, CancerType):
            raise TypeError(
                "'cancer_cohort' must be a CancerType, "
                f"got {type(self.cancer_cohort).__name__}."
            )

        for path_field_name in ("expression_matrix_path", "metadata_path", "output_dir"):
            value = getattr(self, path_field_name)
            if not isinstance(value, Path):
                raise TypeError(
                    f"'{path_field_name}' must be a pathlib.Path, "
                    f"got {type(value).__name__}."
                )

        for column_field_name in (
            "gene_id_column",
            "sample_id_column",
            "group_column",
        ):
            value = getattr(self, column_field_name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    f"'{column_field_name}' must be a non-empty str."
                )

        if not isinstance(self.reference_group, str) or not self.reference_group.strip():
            raise ValueError("'reference_group' must be a non-empty str.")

        if not isinstance(self.comparison_group, str) or not self.comparison_group.strip():
            raise ValueError("'comparison_group' must be a non-empty str.")

        if self.reference_group == self.comparison_group:
            raise ValueError(
                "'reference_group' and 'comparison_group' must differ; "
                f"both were {self.reference_group!r}."
            )

        if not isinstance(self.expression_representation, ExpressionRepresentation):
            raise TypeError(
                "'expression_representation' must be an "
                "ExpressionRepresentation, got "
                f"{type(self.expression_representation).__name__}."
            )

        if not isinstance(self.statistical_method, StatisticalMethod):
            raise TypeError(
                "'statistical_method' must be a StatisticalMethod, "
                f"got {type(self.statistical_method).__name__}."
            )

        if not isinstance(self.multiple_testing_method, MultipleTestingMethod):
            raise TypeError(
                "'multiple_testing_method' must be a "
                "MultipleTestingMethod, got "
                f"{type(self.multiple_testing_method).__name__}."
            )

        if isinstance(self.significance_threshold, bool) or not isinstance(
            self.significance_threshold, (int, float)
        ):
            raise TypeError(
                "'significance_threshold' must be an int or float, "
                f"got {type(self.significance_threshold).__name__}."
            )

        if not (0.0 < self.significance_threshold < 1.0):
            raise ValueError(
                "'significance_threshold' must be strictly between 0 "
                f"and 1, got {self.significance_threshold!r}."
            )

        if self.effect_size_threshold is not None:
            if isinstance(self.effect_size_threshold, bool) or not isinstance(
                self.effect_size_threshold, (int, float)
            ):
                raise TypeError(
                    "'effect_size_threshold' must be an int, float, "
                    f"or None, got {type(self.effect_size_threshold).__name__}."
                )
            if self.effect_size_threshold < 0:
                raise ValueError(
                    "'effect_size_threshold' must be non-negative, "
                    f"got {self.effect_size_threshold!r}."
                )

        if isinstance(self.minimum_replicates_per_group, bool) or not isinstance(
            self.minimum_replicates_per_group, int
        ):
            raise TypeError(
                "'minimum_replicates_per_group' must be an int, "
                f"got {type(self.minimum_replicates_per_group).__name__}."
            )

        if self.minimum_replicates_per_group < ABSOLUTE_MINIMUM_REPLICATES_PER_GROUP:
            raise ValueError(
                "'minimum_replicates_per_group' must be at least "
                f"{ABSOLUTE_MINIMUM_REPLICATES_PER_GROUP} (the minimum "
                "required to compute a two-sample statistic), got "
                f"{self.minimum_replicates_per_group!r}."
            )

        if not isinstance(self.gene_filter, GeneFilterConfiguration):
            raise TypeError(
                "'gene_filter' must be a GeneFilterConfiguration, "
                f"got {type(self.gene_filter).__name__}."
            )

        if self.random_seed is not None and (
            isinstance(self.random_seed, bool) or not isinstance(self.random_seed, int)
        ):
            raise TypeError(
                "'random_seed' must be an int or None, "
                f"got {type(self.random_seed).__name__}."
            )


@dataclass(frozen=True)
class ComparisonDefinition:
    """The explicit two-group comparison being tested.

    Constructed from DEAnalysisConfiguration; kept as its own type so
    comparison.py's output has a name distinct from the full
    configuration, and so downstream modules (methods.py, qc.py,
    reproducibility.py) can depend on just the comparison rather than
    the entire configuration object.
    """

    group_column: str
    reference_group: str
    comparison_group: str

    def __post_init__(self) -> None:
        if self.reference_group == self.comparison_group:
            raise ValueError(
                "'reference_group' and 'comparison_group' must differ; "
                f"both were {self.reference_group!r}."
            )


@dataclass(frozen=True)
class ExcludedSample:
    """Record of one sample excluded from a comparison, with reason
    and stage, per specification Section 8.3.
    """

    sample_id: str
    reason: str
    stage: str


@dataclass(frozen=True)
class GroupAssignment:
    """The resolved outcome of applying a ComparisonDefinition to a
    validated, matched sample set.

    Sample ID tuples preserve deterministic ordering (sorted) so the
    same input always produces the same GroupAssignment.
    """

    reference_group: str
    comparison_group: str
    reference_sample_ids: tuple[str, ...]
    comparison_sample_ids: tuple[str, ...]
    excluded_samples: tuple[ExcludedSample, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.reference_group == self.comparison_group:
            raise ValueError(
                "'reference_group' and 'comparison_group' must differ."
            )

        overlap = set(self.reference_sample_ids) & set(self.comparison_sample_ids)
        if overlap:
            raise ValueError(
                "A sample cannot belong to both groups simultaneously; "
                f"overlapping sample IDs: {sorted(overlap)}."
            )


@dataclass(frozen=True)
class ValidatedExpressionMatrix:
    """A structurally validated expression matrix, produced only by
    validation.py.

    The dataframe field is excluded from dataclass-generated equality
    comparisons because pandas DataFrame.__eq__ returns an
    element-wise result, which is not usable as a boolean and would
    break the default dataclass __eq__.
    """

    file_path: Path
    gene_id_column: str
    sample_columns: tuple[str, ...]
    gene_ids: tuple[str, ...]
    dataframe: pd.DataFrame = field(compare=False, repr=False)


@dataclass(frozen=True)
class ValidatedMetadata:
    """Structurally validated sample metadata, produced only by
    validation.py.
    """

    file_path: Path
    sample_id_column: str
    group_column: str
    dataframe: pd.DataFrame = field(compare=False, repr=False)


@dataclass(frozen=True)
class SampleMatchResult:
    """The outcome of cross-referencing expression-matrix sample
    columns against metadata sample IDs (specification Section 3.4).

    Every tuple is sorted for deterministic ordering.
    """

    matched_samples: tuple[str, ...]
    expression_only_samples: tuple[str, ...]
    metadata_only_samples: tuple[str, ...]


@dataclass(frozen=True)
class GeneFilterResult:
    """The outcome of applying a GeneFilterConfiguration to a
    validated expression matrix, restricted to the samples included
    in the resolved comparison (specification Section 8.2).

    Attributes:
        tested_gene_ids:
            Gene identifiers that will proceed to statistical testing,
            in their original expression-matrix order.

        filtered_gene_ids:
            Gene identifiers removed by the filter, in their original
            expression-matrix order. Never silently dropped: every
            identifier here is also recorded, with the shared
            criterion_description, in QCReport.

        criterion_description:
            The exact filtering criterion applied (or a statement that
            no filtering was applied), copied from
            GeneFilterConfiguration.criterion_description for
            traceability into QC and reproducibility metadata.

        input_gene_count:
            Total number of genes in the expression matrix before
            filtering.
    """

    tested_gene_ids: tuple[str, ...]
    filtered_gene_ids: tuple[str, ...]
    criterion_description: str
    input_gene_count: int

    def __post_init__(self) -> None:
        overlap = set(self.tested_gene_ids) & set(self.filtered_gene_ids)
        if overlap:
            raise ValueError(
                "A gene cannot be both tested and filtered; overlapping "
                f"gene IDs: {sorted(overlap)}."
            )
        if len(self.tested_gene_ids) + len(self.filtered_gene_ids) != self.input_gene_count:
            raise ValueError(
                "tested_gene_ids and filtered_gene_ids must together "
                f"account for all {self.input_gene_count} input gene(s); "
                f"got {len(self.tested_gene_ids)} tested + "
                f"{len(self.filtered_gene_ids)} filtered."
            )


@dataclass(frozen=True)
class GeneResult:
    """The statistical result for one gene.

    Per specification Section 6.4, a gene for which the statistical
    method cannot produce a valid result is retained with missing
    fields recorded as None rather than removed from the output.

    Attributes:
        gene_id:
            The gene identifier, preserved exactly as it appeared in
            the expression matrix.

        tested:
            Whether this gene was actually submitted to the
            statistical method (False for genes removed by explicit
            filtering before testing, per specification Section 8.2).

        effect_size:
            Signed effect-size estimate. Positive means higher
            expression in the comparison group relative to the
            reference group. None if the gene was not tested or the
            method could not produce a value.

        effect_size_label:
            What effect_size actually represents (e.g.
            "log2_fold_change" only when justified by the declared
            expression representation, otherwise a non-fold-change
            label such as "mean_difference"). See methods.py for the
            exact labeling rule.

        raw_p_value:
            Unadjusted p-value from the statistical test. None if not
            tested or not produced.

        adjusted_p_value:
            Multiple-testing-corrected p-value. None if not tested,
            not produced, or excluded from correction (e.g. because
            raw_p_value was missing).

        significant:
            Whether adjusted_p_value is below the analysis's
            configured significance_threshold. None if
            adjusted_p_value is None (significance cannot be
            evaluated without an adjusted p-value).

        missing_reason:
            Human-readable reason a statistical value is missing, when
            known (specification Section 6.4). None if the gene was
            fully tested and produced valid values, or if the
            statistical framework did not report a reason.
    """

    gene_id: str
    tested: bool
    effect_size: float | None
    effect_size_label: str | None
    raw_p_value: float | None
    adjusted_p_value: float | None
    significant: bool | None
    missing_reason: str | None = None


@dataclass(frozen=True)
class QCReport:
    """Structured quality-control report for one analysis, covering
    specification Section 8's pre- and post-analysis requirements.

    Fields are grouped into observed data facts (what was found in
    the input) and configuration decisions (what the analysis was
    configured to do), per specification Section 8.5's requirement
    that QC reporting distinguish the two.
    """

    # --- Observed data facts ---
    input_gene_count: int
    tested_gene_count: int
    filtered_gene_count: int
    initial_sample_count: int
    included_sample_count: int
    excluded_sample_count: int
    excluded_samples: tuple[ExcludedSample, ...]
    reference_group_size: int
    comparison_group_size: int

    # --- Configuration decisions, recorded for traceability ---
    gene_filter_criterion: str
    statistical_method: StatisticalMethod
    multiple_testing_method: MultipleTestingMethod
    significance_threshold: float

    notes: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class AnalysisMetadata:
    """Complete reproducibility metadata for one analysis, per
    specification Section 9 and Section 12.2.
    """

    analysis_id: str
    timestamp_utc: datetime
    cancer_cohort: CancerType

    expression_matrix_path: Path
    metadata_path: Path
    expression_matrix_sha256: str
    metadata_sha256: str

    expression_representation: ExpressionRepresentation
    statistical_method: StatisticalMethod
    statistical_method_version: str | None

    design: str
    reference_group: str
    comparison_group: str
    included_sample_ids: tuple[str, ...]
    reference_group_size: int
    comparison_group_size: int

    gene_filter_criterion: str
    multiple_testing_method: MultipleTestingMethod
    significance_threshold: float
    effect_size_threshold: float | None

    python_version: str
    package_versions: dict[str, str]
    de_package_version: str
    blc_mark_version: str | None

    analysis_status: AnalysisStatus


@dataclass(frozen=True)
class FailureInfo:
    """Structured description of why an analysis failed, per
    specification Sections 10.6 and 15.5.
    """

    stage: str
    category: str
    message: str
    cause: str | None = None


@dataclass(frozen=True)
class AnalysisResult:
    """The complete outcome of running analysis.py's orchestration
    pipeline, whether the analysis succeeded or failed.

    Exactly one of (gene_results, qc_report, metadata) being fully
    populated versus failure being populated distinguishes success
    from failure; analysis.py never returns a result that looks
    successful (populated gene_results/metadata) alongside a
    populated failure, and never returns status=SUCCEEDED without
    all three of gene_results, qc_report, and metadata populated.
    """

    analysis_id: str
    status: AnalysisStatus
    gene_results: tuple[GeneResult, ...] | None = None
    qc_report: QCReport | None = None
    metadata: AnalysisMetadata | None = None
    failure: FailureInfo | None = None
    output_paths: dict[str, Path] | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.status, AnalysisStatus):
            raise TypeError(
                "'status' must be an AnalysisStatus, "
                f"got {type(self.status).__name__}."
            )

        if self.status == AnalysisStatus.SUCCEEDED:
            if self.gene_results is None or self.qc_report is None or self.metadata is None:
                raise ValueError(
                    "A SUCCEEDED AnalysisResult must have gene_results, "
                    "qc_report, and metadata all populated; a failed or "
                    "incomplete analysis must not be represented as "
                    "successful (specification Section 15.5)."
                )
            if self.failure is not None:
                raise ValueError(
                    "A SUCCEEDED AnalysisResult must not also carry "
                    "failure information."
                )

        if self.status == AnalysisStatus.FAILED:
            if self.failure is None:
                raise ValueError(
                    "A FAILED AnalysisResult must carry FailureInfo "
                    "explaining what went wrong."
                )