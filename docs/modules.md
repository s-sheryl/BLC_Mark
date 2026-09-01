# BLC Mark V1 Modules

BLC Mark Version 1 is organized as a modular scientific workflow for
transcriptomic cancer biomarker discovery, evidence integration, and
transparent candidate prioritization.

The executed Version 1 analysis includes three TCGA cohorts:

- TCGA-BRCA — Breast Invasive Carcinoma
- TCGA-LUAD — Lung Adenocarcinoma
- TCGA-COAD — Colon Adenocarcinoma

Version 1 is limited to RNA-seq gene-expression analysis.

---

## 1. Data Acquisition and Dataset Management

**Purpose:**
Acquire, register, organize, and preserve provenance for public
transcriptomic datasets used by BLC Mark.

**Responsibilities:**

- dataset acquisition;
- dataset registration;
- download management;
- file-integrity handling;
- metadata management;
- source and dataset provenance.

**Output:**

- registered source datasets;
- dataset metadata;
- reproducible input records.

---

## 2. Data Validation and Preparation

**Purpose:**
Validate downloaded datasets and prepare expression and sample metadata
for downstream analysis.

**Responsibilities:**

- structural validation;
- expression-data validation;
- metadata validation;
- sample matching;
- preprocessing;
- representation tracking;
- quality-control reporting.

**Output:**

- validated analysis-ready expression data;
- corresponding sample metadata;
- preprocessing and validation records.

---

## 3. Differential Expression and Candidate Discovery

**Purpose:**
Identify genes showing statistically significant expression differences
between predefined biological groups within each cancer cohort.

**Version 1 comparison:**

Primary Tumor versus Solid Tissue Normal.

**Executed Version 1 method:**

- normalized log2 expression;
- Welch's unequal-variance two-sample t-test;
- Benjamini-Hochberg false-discovery-rate correction;
- adjusted p-value threshold of 0.05;
- no default effect-size threshold.

**Responsibilities:**

- explicit comparison configuration;
- statistical testing;
- multiple-testing correction;
- effect-size calculation;
- candidate identification;
- analysis metadata;
- quality-control reporting;
- explicit failure handling.

**Output:**

- complete gene-level differential-expression results;
- significant candidate genes;
- analysis metadata;
- QC reports.

A significant gene is treated as a candidate and not as a clinically
validated biomarker.

---

## 4. Evidence Integration

**Purpose:**
Collect additional evidence for significant Phase 3 candidate genes.

**Version 1 evidence categories include:**

- cancer-association evidence;
- clinical/prognostic evidence;
- cross-cancer evidence;
- functional annotation;
- pathway context.

Functional and pathway evidence are retained for interpretation and
provenance but do not directly contribute to the Version 1 final
prioritization score.

**Responsibilities:**

- external evidence retrieval and integration;
- gene-identifier handling;
- evidence provenance;
- evidence availability tracking;
- quality-control reporting.

**Output:**

- cohort-level evidence tables;
- evidence metadata;
- QC reports.

Unavailable evidence is distinguished from evaluated evidence that
provides no positive support.

---

## 5. Biomarker Prioritization

**Purpose:**
Rank significant candidate genes using an explicit and reproducible
evidence-scoring framework.

**Version 1 scoring components:**

1. differential-expression strength;
2. direct cancer-association evidence;
3. clinical/prognostic evidence;
4. cross-cancer recurrence.

Each component has a fixed Version 1 weight of 0.25.

**Responsibilities:**

- component normalization;
- evidence scoring;
- missing-evidence handling;
- final-score calculation;
- deterministic ranking;
- scoring provenance;
- quality-control reporting.

**Output:**

- prioritized candidate tables;
- component scores;
- final prioritization scores;
- deterministic ranks;
- metadata;
- QC reports.

A high BLC Mark score represents stronger support under the Version 1
prioritization framework. It does not represent clinical validation,
diagnostic accuracy, therapeutic efficacy, or probability of clinical
usefulness.

---

## 6. Scientific Outputs and Reproducibility

**Purpose:**
Transform the completed scientific analyses into validated,
traceable, and reproducible Version 1 outputs.

**Responsibilities:**

- final scientific-output validation;
- generation of top-candidate tables;
- generation of scientific figures;
- generation of the scientific results report;
- generation of the reproducibility manifest;
- SHA-256 integrity recording;
- reproducible regeneration of Phase 6 outputs.

**Output:**

- final biomarker tables;
- scientific figures;
- scientific results report;
- reproducibility manifest.

Phase 6 outputs can be regenerated from the existing frozen upstream
results with:

    python scripts/run_phase6.py

This command regenerates Phase 6 outputs. It does not reacquire external
evidence or rerun Phases 1–5.

---

## 7. Researcher-Facing Application

BLC Mark Version 1 includes a Streamlit application for exploring the
completed scientific results.

The application is intentionally separated from the scientific analysis
engine. It reads the validated Version 1 outputs and does not recalculate
differential expression, retrieve external evidence, or modify biomarker
prioritization.

**Current capabilities include:**

- project overview;
- biomarker exploration;
- search and filtering of candidate genes;
- cohort-level result exploration;
- cross-cohort comparison;
- methodology information;
- downloadable generated results.

The application is implemented in:

    streamlit_app.py

It can be started locally with:

    streamlit run streamlit_app.py

Keeping the application layer separate from the scientific pipeline allows
the validated analysis outputs to remain unchanged while providing a more
accessible interface for researchers exploring the results.

Future application development may extend deployment, interactivity, or
data-access capabilities without changing the frozen Version 1 scientific
analysis.
