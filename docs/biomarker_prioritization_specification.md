 # Phase 5 — Biomarker Prioritization Specification

Version 1.0

## 1. Ranking Unit

One significant Phase 3 candidate gene within one cancer cohort.

## 2. Eligibility

Only Phase 3 significant candidates are prioritized.

No new p-value or fold-change threshold is introduced.

## 3. Scoring Components

1. Differential-expression strength
2. Direct cancer-association evidence
3. Clinical/prognostic evidence
4. Cross-cancer recurrence

## 4. Context-Only Evidence


- NCBI functional annotation
- Reactome pathways

Context-only evidence is retained for interpretation and provenance but does not directly contribute to the Version 1 prioritization score.

## 5. Component Normalization

### 5.1 Differential-Expression Strength

The differential-expression component is the percentile rank of the absolute log2 fold change among significant Phase 3 candidates within the same cancer cohort.

The original signed effect size is retained separately for biological interpretation.

### 5.2 Cancer-Association Evidence

The Open Targets cancer-association score is used as the cancer-association component.

Expected range: 0 to 1.

### 5.3 Clinical/Prognostic Evidence

Human Protein Atlas prognostic categories are mapped as follows:

- unprognostic = 0.0
- potential prognostic = 0.5
- validated prognostic = 1.0

Favorable versus unfavorable direction does not change evidence strength. The direction is retained separately for interpretation.

### 5.4 Cross-Cancer Recurrence

Cross-cancer recurrence is mapped as follows:

- significant in one cohort = 0.0
- significant in two cohorts = 0.5
- significant in three cohorts = 1.0

This component represents additional cross-cancer support rather than basic Phase 3 eligibility.

## 6. Version 1 Weights

The four scoring components receive equal fixed weights:

- Differential-expression strength = 0.25
- Cancer-association evidence = 0.25
- Clinical/prognostic evidence = 0.25
- Cross-cancer recurrence = 0.25

The weights are fixed in Version 1.

Configurable evidence weights are deferred to Version 2.

## 7. Final Prioritization Score

For candidates whose four scoring components are evaluable, the final prioritization score is the weighted sum of the normalized component scores.

Expected range: 0 to 1.

## 8. Missing and Unavailable Evidence

Failure to find supporting evidence must not automatically be interpreted as evidence that a gene has no biological importance.

For a resolvable candidate where a source was successfully evaluated but no supporting evidence record was captured, that component provides no positive support.

Legacy unresolved identifiers of the form `?|<id>` are retained as candidates. External evidence that could not be evaluated for such identifiers must be recorded explicitly as unavailable rather than silently interpreted as biological negative evidence.

No candidate is silently discarded because external evidence is absent or unavailable.

## 9. Interpretation

A high prioritization score indicates stronger support for further investigation under the BLC Mark Version 1 prioritization framework.

The score does not represent:

- probability of clinical usefulness
- clinical validation
- diagnostic accuracy
- therapeutic efficacy
- causal cancer-driver status

## 10. Ranking

Candidates are ordered by final prioritization score from highest to lowest.

Raw component inputs and normalized component scores must remain available so that every resulting rank is explainable.

Candidates are ranked by final prioritization score from highest to lowest. Exact score ties are resolved deterministically by gene_id in ascending lexical order.

## 11. Reproducibility

Phase 5 outputs must record:

- scoring version
- fixed component weights
- raw component inputs
- normalized component scores
- final prioritization score
- rank
- evidence availability status
- cancer cohort
- Phase 3 provenance
- Phase 4 provenance
