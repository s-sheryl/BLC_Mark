from functools import lru_cache
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from src.evidence_integration.analysis import run_phase4_analysis
from src.evidence_integration.cross_cancer import build_cross_cancer_index
from src.evidence_integration.hpa_prognostic import load_hpa_prognostic_data
from src.evidence_integration.ncbi_gene import (
    collect_functional_evidence,
    find_human_gene_id,
)
from src.evidence_integration.open_targets_identifiers import (
    resolve_gene_symbol_to_ensembl,
)
from src.evidence_integration.reactome import load_reactome_mapping

@lru_cache(maxsize=None)
def cached_ncbi_gene_id(gene_symbol: str) -> str | None:
    return find_human_gene_id(gene_symbol)


@lru_cache(maxsize=None)
def cached_ensembl_id(gene_symbol: str) -> str | None:
    return resolve_gene_symbol_to_ensembl(gene_symbol)


@lru_cache(maxsize=None)
def cached_functional_evidence(
    gene_symbol: str,
    retrieved_at: str,
):
    return collect_functional_evidence(
        gene_symbol,
        retrieved_at=retrieved_at,
    )
def main():
    temp_path = Path("results/phase4_smoke_input.csv")

    rows = [
        {
            "gene_id": "TP53",
            "tested": True,
            "effect_size": 1.0,
            "effect_size_label": "difference_in_group_means",
            "raw_p_value": 0.001,
            "adjusted_p_value": 0.01,
            "significant": True,
            "missing_reason": None,
        },
        {
            "gene_id": "BRCA1",
            "tested": True,
            "effect_size": 1.0,
            "effect_size_label": "difference_in_group_means",
            "raw_p_value": 0.001,
            "adjusted_p_value": 0.01,
            "significant": True,
            "missing_reason": None,
        },
        {
            "gene_id": "HIF3A",
            "tested": True,
            "effect_size": 1.0,
            "effect_size_label": "difference_in_group_means",
            "raw_p_value": 0.001,
            "adjusted_p_value": 0.01,
            "significant": True,
            "missing_reason": None,
        },
    ]

    pd.DataFrame(rows).to_csv(temp_path, index=False)

    reactome_mapping = load_reactome_mapping(
        "data/external/reactome/NCBI2Reactome_V97.txt"
    )

    hpa_index = load_hpa_prognostic_data(
        "data/external/hpa/cancer_prognostic_data_HPA25.1.tsv.zip"
    )

    cross_cancer_index = build_cross_cancer_index(
        {
            "TCGA-BRCA": {"TP53", "BRCA1", "HIF3A"},
            "TCGA-LUAD": {"TP53", "HIF3A"},
            "TCGA-COAD": {"TP53"},
        }
    )

    retrieved_at = datetime.now(timezone.utc).isoformat()

    result = run_phase4_analysis(
        temp_path,
        "TCGA-BRCA",
        cross_cancer_index=cross_cancer_index,
        reactome_mapping=reactome_mapping,
        hpa_index=hpa_index,
        retrieved_at=retrieved_at,
        reactome_version="97",
        hpa_version="25.1",
        open_targets_version="26.06",
        resolve_ncbi_gene_id=cached_ncbi_gene_id,
        resolve_ensembl_id=cached_ensembl_id,
        collect_functional=cached_functional_evidence,
        open_targets_page_size=10,
    )

    print("Candidate count:", result.candidate_count)
    print(
        "Unresolved identifiers:",
        result.unresolved_identifier_count,
    )

    for profile in result.profiles:
        print()
        print("GENE:", profile.gene_id)
        print("EVIDENCE:", len(profile.evidence_records))

        for record in profile.evidence_records[:10]:
            print(
                record.evidence_type.value,
                "|",
                record.source,
                "|",
                record.evidence_id,
            )


if __name__ == "__main__":
    main()