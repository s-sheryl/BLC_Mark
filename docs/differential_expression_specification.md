# BLC Mark — Differential Expression Analysis Specification

**Document:** Differential Expression Analysis Specification
**Project:** BLC Mark
**Version:** 1.0
**Status:** Implemented and frozen for Version 1
**Scope:** Version 1

---

## 1. Purpose

This document defines the scientific and engineering specification for differential-expression analysis in BLC Mark Version 1.

The differential-expression module is responsible for identifying genes whose expression differs between explicitly defined biological groups while preserving statistical transparency, reproducibility, validation, and traceability.

The module is designed as a candidate-discovery stage.

A statistically significant gene identified by this module is considered a **candidate for downstream evidence integration**. Differential expression alone is not interpreted as clinical biomarker validation.

---

## 2. Version 1 Scientific Scope

BLC Mark Version 1 performs differential-expression analysis on RNA-seq gene-expression data.

The executed Version 1 analysis contains three TCGA cohorts:

| Cancer type | Cohort |
| --- | --- |
| Breast Cancer | TCGA-BRCA |
| Lung Adenocarcinoma | TCGA-LUAD |
| Colorectal Cancer | TCGA-COAD |

The implemented comparison for each executed cohort is:

```text
Primary Tumor
vs
Solid Tissue Normal
```

Version 1 is intentionally limited to explicitly defined two-group comparisons.

The following designs are outside the Version 1 scope:

- multi-group differential-expression models;
- paired analyses;
- repeated-measures designs;
- interaction models;
- complex multifactor designs;
- survival models;
- longitudinal expression models.

These analyses may be considered in future versions but are not silently approximated by the Version 1 implementation.

---

## 3. Scientific Question

The differential-expression stage addresses the following question:

> Which genes show statistically significant expression differences between Primary Tumor and Solid Tissue Normal samples within each analyzed cancer cohort?

The resulting significant genes are passed to the evidence-integration stage for further biological and clinical evaluation.

---

## 4. Input Requirements

Differential-expression analysis requires two validated inputs:

1. a processed gene-expression matrix;
2. matching sample metadata.

The expression matrix must contain gene-level expression measurements for the samples included in the analysis.

The metadata must provide sufficient information to assign samples to the explicitly configured biological comparison groups.

Input data must be validated before statistical analysis begins.

The differential-expression module must not silently infer missing scientific metadata or reconstruct comparison groups from assumptions.

---

## 5. Gene Identifiers

Gene identifiers are preserved through the differential-expression stage.

The module does not silently remap gene identifiers during statistical analysis.

Identifier resolution required by external biological evidence sources is handled separately during downstream evidence integration.

This separation prevents identifier-conversion decisions from silently changing the statistical candidate set.

---

## 6. Expression Representation

The statistical method used for differential expression depends on the representation of the expression data.

BLC Mark explicitly tracks expression representation rather than attempting to infer it from numerical values.

Representations recognized by the project include:

```text
RAW_COUNTS
NORMALIZED_LOG2
NORMALIZED_LINEAR
```

The executed Version 1 analyses use:

```text
NORMALIZED_LOG2
```

The representation supplied to the analysis must be compatible with the selected statistical method.

An incompatible representation-method combination must fail explicitly.

---

## 7. Statistical Method

For the executed Version 1 analyses, differential expression is evaluated using:

**Welch's unequal-variance two-sample t-test**

The test is two-sided.

Welch's t-test was selected for the executed workflow because the Version 1 datasets are analyzed as normalized log2 expression values rather than raw RNA-seq counts.

The method does not assume equal variance between the two comparison groups.

---

## 8. Raw-Count Data

Raw RNA-seq counts represent a different statistical analysis problem from normalized log2 expression.

Count-specific approaches such as DESeq2 are more appropriate for raw RNA-seq count data.

BLC Mark Version 1 does not silently apply Welch's t-test to data declared as raw counts.

DESeq2 is not implemented as an available Version 1 analysis method.

If a requested analysis requires a method that is not implemented, the system must fail explicitly rather than substitute another method.

This behavior preserves the scientific meaning of the requested analysis.

---

## 9. Comparison Definition

Every differential-expression analysis must explicitly define:

- the cancer cohort;
- the biological comparison;
- the reference group;
- the comparison group;
- the expression representation;
- the statistical method;
- the multiple-testing procedure;
- the significance threshold.

For the executed Version 1 analyses, the comparison is:

```text
Reference group:
Solid Tissue Normal

Comparison group:
Primary Tumor
```

Effect-size direction is interpreted according to this configured comparison.

Positive and negative effect sizes are both retained.

---

## 10. Effect Size

For the executed normalized-log2 analyses, the reported expression effect is represented as log2 fold change.

Effect size is retained independently from statistical significance.

Both positive and negative expression changes may be biologically meaningful.

The differential-expression stage does not discard a statistically significant candidate solely because its effect is negative.

---

## 11. Statistical Testing

For each testable gene, the module performs a two-sided Welch two-sample t-test between the configured biological groups.

The analysis produces, where statistically testable:

- an effect size;
- a raw p-value;
- an adjusted p-value.

Genes that cannot produce a valid statistical result are handled explicitly rather than silently removed or assigned fabricated statistical values.

---

## 12. Multiple-Testing Correction

Thousands of genes are evaluated simultaneously.

Raw p-values are therefore corrected for multiple testing.

The implemented Version 1 correction method is:

**Benjamini-Hochberg false-discovery-rate correction**

The significance criterion is:

```text
adjusted p-value < 0.05
```

A gene satisfying this criterion is considered a significant differential-expression candidate.

---

## 13. Candidate Selection

Version 1 candidate selection is based on:

```text
adjusted p-value < 0.05
```

No additional default effect-size cutoff is applied.

This decision keeps the differential-expression stage focused on statistical candidate discovery.

Effect-size magnitude is retained in the output and becomes relevant during downstream prioritization, but it is not silently introduced as an additional Phase 3 significance criterion.

---

## 14. Gene Filtering

The executed Version 1 workflow does not apply an additional default gene-expression filtering threshold before candidate selection.

Filtering behavior must be explicit and configuration-driven.

The analysis must not silently introduce a new expression threshold, fold-change threshold, or candidate-selection rule.

---

## 15. Version 1 Cohort Results

The frozen Version 1 differential-expression analysis produced the following results:

| Cohort | Genes tested | Significant candidates |
| --- | ---: | ---: |
| TCGA-BRCA | 20,530 | 16,644 |
| TCGA-LUAD | 20,530 | 16,031 |
| TCGA-COAD | 20,530 | 15,468 |

These candidate counts are the frozen Phase 3 outputs used by the downstream evidence-integration and biomarker-prioritization stages.

The large number of significant genes is one of the motivations for the later BLC Mark evidence-integration and prioritization workflow.

Differential-expression significance is therefore treated as the beginning of candidate evaluation rather than its endpoint.

---

## 16. Executed Sample Counts

The executed Version 1 comparisons included the following samples:

| Cohort | Solid Tissue Normal | Primary Tumor | Included samples |
| --- | ---: | ---: | ---: |
| TCGA-BRCA | 114 | 1,097 | 1,211 |
| TCGA-LUAD | 59 | 515 | 574 |
| TCGA-COAD | 41 | 286 | 327 |

Samples excluded during validation or comparison construction were recorded rather than silently incorporated into the analysis.

The frozen runs recorded:

| Cohort | Excluded samples |
| --- | ---: |
| TCGA-BRCA | 7 |
| TCGA-LUAD | 2 |
| TCGA-COAD | 2 |

---

## 17. Non-Testable or Missing Statistical Results

Not every gene necessarily produces a valid adjusted p-value.

The frozen Version 1 runs recorded missing adjusted p-values for:

| Cohort | Missing adjusted p-values |
| --- | ---: |
| TCGA-BRCA | 280 |
| TCGA-LUAD | 333 |
| TCGA-COAD | 480 |

These values are preserved transparently.

A missing or unavailable statistical result must not be converted automatically into a statistically meaningful zero or another fabricated value.

This distinction is important for downstream reproducibility and interpretation.

---

## 18. Required Result Table

The differential-expression stage produces a machine-readable per-gene result table.

At minimum, each tested gene must be associated with the information required to interpret its statistical result, including:

- gene identifier;
- effect size;
- raw p-value;
- adjusted p-value.

Additional structured fields may be retained where required by the implementation.

Results must remain traceable to the analysis that generated them.

---

## 19. Analysis Metadata

Each differential-expression execution must produce metadata describing the analysis.

Metadata should record information required to reconstruct and interpret the run, including:

- project/version information;
- analysis identifier;
- cancer cohort;
- expression representation;
- comparison groups;
- statistical method;
- multiple-testing method;
- significance threshold;
- relevant input information;
- execution information;
- relevant software or method versions where recorded.

The metadata is part of the scientific output and not optional documentation.

---

## 20. Quality Control

The differential-expression stage must produce quality-control information sufficient to audit the execution.

QC information includes relevant values such as:

- number of genes tested;
- number of significant candidates;
- number of samples included;
- sample counts by comparison group;
- excluded samples;
- missing or non-testable statistical results;
- analysis completion status.

Quality-control information must correspond to the same analysis identifier as the statistical result table and metadata.

---

## 21. Reproducibility

Reproducibility is a core requirement of the differential-expression module.

A repeated analysis using the same validated inputs, configuration, software behavior, and frozen dependencies should produce equivalent scientific results.

The module therefore avoids hidden scientific decisions.

The following must remain explicit:

- expression representation;
- comparison groups;
- comparison direction;
- statistical method;
- significance threshold;
- multiple-testing procedure;
- filtering configuration.

The implementation must not silently change any of these values during execution.

---

## 22. Explicit Failure Behavior

Scientific failures must be visible.

The module must fail explicitly when an analysis cannot be performed under the requested configuration.

Examples include:

- missing required input files;
- invalid expression data;
- invalid metadata;
- incompatible expression representation;
- invalid comparison configuration;
- insufficient group replication;
- unsupported statistical method;
- statistical execution failure;
- invalid result generation;
- unwritable output locations.

A failed or incomplete analysis must not be presented as a successful differential-expression run.

---

## 23. No Silent Method Substitution

BLC Mark does not automatically replace an unavailable statistical method with another method.

For example, if an analysis requires a raw-count method that is not implemented, the system must not silently substitute Welch's t-test.

Likewise, the system must not silently:

- change significance thresholds;
- change comparison groups;
- reverse group definitions;
- add an effect-size cutoff;
- infer missing sample labels;
- remove problematic observations without recording them;
- change the declared expression representation.

Explicit failure is preferred over producing scientifically ambiguous results.

---

## 24. TCGA-OV Evaluation

TCGA-OV was considered during Version 1 development.

The intended Version 1 comparison was kept consistent across cancer cohorts:

```text
Primary Tumor
vs
Solid Tissue Normal
```

The selected TCGA-OV dataset did not provide the required Solid Tissue Normal samples needed to perform that comparison under the frozen Version 1 design.

The project therefore did not:

- treat recurrent tumor samples as normal controls;
- substitute another tissue category;
- introduce an external normal cohort without a defined harmonization strategy;
- change the biological comparison only for TCGA-OV.

TCGA-OV was consequently excluded from the executed Version 1 analysis.

The final executed Version 1 cohort set is therefore:

```text
TCGA-BRCA
TCGA-LUAD
TCGA-COAD
```

TCGA-OV may be reconsidered in a future version if an appropriate scientifically justified comparison design is defined.

---

## 25. Relationship to Evidence Integration

Differential expression is not the final BLC Mark output.

Significant Phase 3 candidates are passed to the evidence-integration stage.

Phase 4 evaluates additional information including:

- cancer-association evidence;
- clinical/prognostic evidence;
- cross-cancer evidence;
- functional annotation;
- pathway context.

This separation ensures that statistical evidence from the expression analysis is not confused with independent biological or clinical evidence.

---

## 26. Relationship to Biomarker Prioritization

The effect sizes and statistical results generated during differential expression contribute to downstream biomarker prioritization.

However, Phase 3 itself does not assign the final BLC Mark biomarker-prioritization score.

The final prioritization framework is implemented separately so that:

- statistical candidate discovery;
- external evidence integration;
- evidence scoring;
- final ranking

remain distinct scientific responsibilities.

---

## 27. Interpretation

A significant differential-expression result means that the analyzed expression data provide statistical evidence of a difference between the configured biological groups under the implemented statistical model.

It does not establish:

- diagnostic performance;
- prognostic validity;
- therapeutic relevance;
- causal biological function;
- treatment response;
- clinical utility;
- independent experimental validation.

Genes identified during this phase are therefore described as **differential-expression candidates**, not clinically validated biomarkers.

---

## 28. Version 1 Limitations

The implemented differential-expression workflow has several important limitations.

### RNA-seq gene expression only

Version 1 does not incorporate mutation, copy-number, methylation, proteomic, metabolomic, or other molecular layers.

### Normalized-log2 statistical design

The executed analyses use Welch's t-test on normalized log2 expression data.

A raw-count RNA-seq workflow using a count-specific model represents a different analysis design and may produce different results.

### Two-group comparisons

Version 1 is restricted to explicitly configured two-group comparisons.

More complex experimental designs are outside scope.

### No default effect-size threshold

Candidate discovery is based on adjusted statistical significance without an additional default fold-change cutoff.

This produces large candidate sets and makes downstream evidence integration and prioritization important.

### Cohort availability

Only cohorts supporting the defined Primary Tumor versus Solid Tissue Normal comparison were retained in the executed Version 1 analysis.

---

## 29. Frozen Version 1 Configuration

The executed and frozen Version 1 differential-expression configuration is summarized below.

| Parameter | Version 1 configuration |
| --- | --- |
| Data type | RNA-seq gene expression |
| Executed cohorts | TCGA-BRCA, TCGA-LUAD, TCGA-COAD |
| Comparison | Primary Tumor vs Solid Tissue Normal |
| Expression representation | Normalized log2 |
| Statistical method | Welch's unequal-variance two-sample t-test |
| Test direction | Two-sided |
| Multiple-testing correction | Benjamini-Hochberg |
| Significance threshold | Adjusted p-value < 0.05 |
| Default effect-size cutoff | None |
| Additional default gene filtering | None |
| Unsupported raw-count method | Explicit failure rather than substitution |

These settings define the completed BLC Mark Version 1 differential-expression analysis.

They should not be changed retroactively unless a genuine scientific or implementation defect is identified and the resulting analysis is explicitly versioned accordingly.

---

## 30. Version 1 Status

The differential-expression implementation and executed Version 1 results are complete.

The frozen Phase 3 outputs serve as upstream scientific inputs for:

```text
Phase 4 — Evidence Integration
        ↓
Phase 5 — Biomarker Prioritization
        ↓
Phase 6 — Scientific Outputs & Reproducibility
```

Version 1 therefore treats the differential-expression stage as a completed and reproducible component of the BLC Mark scientific workflow.
