# BLC Mark Technology Stack

## Implemented Version 1 Scientific Stack

### Programming Language
- Python 3.14

### Core Data Analysis
- pandas
- NumPy
- SciPy

### Statistical Analysis
- Welch's unequal-variance two-sample t-test
- Benjamini-Hochberg false-discovery-rate correction

### Visualization
- Matplotlib

### Data Formats
- CSV
- JSON
- Markdown
- PNG

### External Data and Evidence Access
- Public TCGA/Xena transcriptomic resources
- Open Targets Platform
- Human Protein Atlas
- NCBI functional annotation resources
- Reactome pathway resources

### Software Engineering
- pathlib
- dataclasses
- structured configuration and metadata handling
- SHA-256 file integrity tracking
- explicit validation and failure handling
- deterministic ranking

### Testing
- pytest

### Development Environment
- Visual Studio Code
- Git
- GitHub
- Python virtual environment

---

## Reproducibility

BLC Mark Version 1 records:

- analysis configuration;
- dataset and input provenance;
- software versions;
- quality-control information;
- evidence provenance;
- scoring configuration;
- output hashes.

Phase 6 scientific outputs can be regenerated with:

    python scripts/run_phase6.py

This command regenerates Phase 6 outputs from the existing frozen upstream
results and does not rerun data acquisition, differential expression,
evidence retrieval, or biomarker prioritization.

---

## Post-V1 Application Stack

The researcher-facing web application is outside the completed Version 1
scientific engine.

Technologies being considered for the application layer include:

- FastAPI for backend/API services;
- Streamlit or another suitable frontend framework;
- SQLite for lightweight local persistence;
- PostgreSQL for larger deployments;
- Plotly for interactive visualizations;
- Docker for deployment and reproducibility.

These are planned application technologies and should not be interpreted
as implemented Version 1 scientific dependencies.

---

## Future Scientific Technologies

Possible future scientific extensions may require additional tools for:

- survival analysis;
- independent-cohort validation;
- multi-omics integration;
- machine learning;
- advanced genomic analysis.

These technologies are intentionally outside the Version 1 scientific
scope unless formally added to a future specification.