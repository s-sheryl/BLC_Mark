# BLC Mark V1 Completion Checklist

## Scientific Scope

- [x] Define the Version 1 scientific question
- [x] Limit Version 1 to RNA-seq gene-expression analysis
- [x] Execute three TCGA cancer cohorts:
  - TCGA-BRCA
  - TCGA-LUAD
  - TCGA-COAD
- [x] Document the exclusion of TCGA-OV from the executed Version 1 cohort set

## Phase 1 — Scientific & Engineering Foundation

- [x] Define scientific scope and project boundaries
- [x] Establish reproducibility requirements
- [x] Establish explicit validation and failure behaviour
- [x] Define modular project architecture

## Phase 2 — Data Acquisition & Preparation

- [x] Acquire public cancer transcriptomic data
- [x] Register dataset metadata
- [x] Validate input datasets
- [x] Prepare analysis-ready expression and sample metadata
- [x] Preserve dataset provenance

## Phase 3 — Biomarker Discovery

- [x] Implement differential-expression analysis
- [x] Define explicit tumor-versus-normal comparisons
- [x] Use normalized log2 expression for the executed V1 analyses
- [x] Use Welch's unequal-variance two-sample t-test
- [x] Apply Benjamini-Hochberg multiple-testing correction
- [x] Retain effect sizes and adjusted p-values
- [x] Generate analysis metadata and QC reports
- [x] Preserve explicit failure behaviour
- [x] Generate candidate genes for all three executed cohorts

## Phase 4 — Evidence Integration

- [x] Integrate cancer-association evidence
- [x] Integrate clinical/prognostic evidence
- [x] Integrate cross-cancer evidence
- [x] Retain functional annotation
- [x] Retain pathway context
- [x] Preserve evidence provenance
- [x] Distinguish unavailable evidence from biological negative evidence
- [x] Generate evidence and QC outputs for all three cohorts

## Phase 5 — Biomarker Prioritization

- [x] Implement transparent evidence scoring
- [x] Use four explicitly defined scoring components
- [x] Use fixed Version 1 component weights
- [x] Preserve raw and normalized component evidence
- [x] Implement deterministic ranking
- [x] Retain candidates with unavailable scores
- [x] Generate prioritization results, metadata and QC reports
- [x] Generate ranked candidates for all three cohorts

## Phase 6 — Scientific Outputs & Reproducibility

- [x] Validate final scientific outputs
- [x] Generate final biomarker tables
- [x] Generate scientific figures
- [x] Generate the Version 1 scientific results report
- [x] Generate a SHA-256 reproducibility manifest
- [x] Provide one-command regeneration of Phase 6 outputs
- [ ] Complete final repository and documentation cleanup
- [ ] Run final release verification and full test suite
- [ ] Create final Version 1 Git commit and release tag

## Version 1 Software Quality

- [x] Modular scientific codebase
- [x] Automated tests
- [x] Explicit validation
- [x] Explicit failure behaviour
- [x] Machine-readable scientific outputs
- [x] Quality-control records
- [x] Reproducibility metadata
- [x] Deterministic biomarker ranking
- [x] Traceable evidence integration
- [ ] Final documentation review
- [ ] Final repository release verification

## Post-V1 Application Development

The following features are intentionally outside the completed Version 1
scientific engine and may be developed as a separate application layer:

- [ ] Web application
- [ ] Interactive biomarker dashboard
- [ ] Search by gene
- [ ] Search by cancer cohort
- [ ] Interactive result exploration
- [ ] User-facing result export
- [ ] Deployment for external researchers

## Future Scientific Extensions

The following are not part of Version 1:

- [ ] Independent-cohort validation
- [ ] Survival-analysis extensions
- [ ] Additional cancer cohorts
- [ ] Additional evidence sources
- [ ] Multi-omics integration
- [ ] Alternative or configurable prioritization models
- [ ] Experimental or clinical biomarker validation