"""
BLC Mark - Dataset Registry

Purpose:
    Single source of truth for every dataset BLC Mark pulls in V1.
    Nothing downloads data or touches the network from this file -- it's
    pure metadata. data_collection.py will import DATASET_REGISTRY and
    iterate over it instead of hardcoding one cancer type at a time.

Why this exists as its own file (and not just a dict inside
data_collection.py):
    The Dataset Explorer feature (see project_vision.md) needs exactly
    this information -- source, cancer type, sample count, platform,
    publication -- to populate its per-cancer dataset table. Keeping it
    here means the dashboard and the collection script both read from
    the same place, so they can't drift out of sync.

Scope reminder (V1):
    Only these three cancers. Do not add more without updating the
    scope decision in master_blueprint.md first.
"""

from dataclasses import dataclass, field
from enum import Enum


class CancerType(str, Enum):
    """The three TCGA cohorts in scope for BLC Mark V1."""

    BRCA = "TCGA-BRCA"  # Breast Invasive Carcinoma
    LUAD = "TCGA-LUAD"  # Lung Adenocarcinoma
    COAD = "TCGA-COAD"  # Colon Adenocarcinoma


class DataType(str, Enum):
    """Assay categories we care about for biomarker discovery."""

    GENE_EXPRESSION_RNASEQ = "gene_expression_rnaseq"
    SOMATIC_MUTATION = "somatic_mutation"
    CLINICAL = "clinical"
    COPY_NUMBER = "copy_number"


@dataclass(frozen=True)
class DatasetEntry:
    """
    Metadata for one downloadable dataset.

    Deliberately does NOT include a hardcoded .gz download URL yet.
    Xena's direct file URLs are queried through their API (XenaQuery)
    rather than guessed by hand -- typing one out from memory risks
    silently pointing at a stale or wrong file, which is worse than no
    URL at all for a biology project. The `datapage_url` lets a human
    verify the dataset by eye in the meantime.

    TODO(data_collection.py): resolve `resolved_download_url` at
    runtime via the Xena hub API before this entry is used to fetch
    a real file. Track that resolution logic in the collection
    script, not here -- this file stays pure metadata.
    """

    cancer_type: CancerType
    data_type: DataType
    cohort_name: str          # Human-readable Xena cohort name
    xena_dataset_id: str      # e.g. "TCGA.BRCA.sampleMap/HiSeqV2"
    xena_host: str            # e.g. "https://tcga.xenahubs.net"
    datapage_url: str         # Browsable page to verify the dataset by eye
    platform: str             # Sequencing/array platform
    expected_samples: int | None = None  # None until confirmed post-download
    notes: str = ""


# ---------------------------------------------------------------------------
# The registry itself.
#
# Sample counts and download URLs are left as None / TODO where we haven't
# independently verified them against the live Xena hub yet. Filling these
# in with made-up numbers would violate the "never fake biological results"
# rule -- better to leave an honest gap than a confident-looking lie.
# ---------------------------------------------------------------------------

DATASET_REGISTRY: list[DatasetEntry] = [
    DatasetEntry(
        cancer_type=CancerType.BRCA,
        data_type=DataType.GENE_EXPRESSION_RNASEQ,
        cohort_name="TCGA Breast Cancer (BRCA)",
        xena_dataset_id="TCGA.BRCA.sampleMap/HiSeqV2",
        xena_host="https://tcga.xenahubs.net",
        datapage_url=(
            "https://xenabrowser.net/datapages/"
            "?cohort=TCGA%20Breast%20Cancer%20(BRCA)"
        ),
        platform="Illumina HiSeq RNASeqV2 (RSEM, log2 normalized)",
        expected_samples=None,  # TODO: confirm against live hub metadata
        notes="RSEM-normalized gene-level expression, PANCAN pipeline.",
    ),
    DatasetEntry(
        cancer_type=CancerType.LUAD,
        data_type=DataType.GENE_EXPRESSION_RNASEQ,
        cohort_name="TCGA Lung Adenocarcinoma (LUAD)",
        xena_dataset_id="TCGA.LUAD.sampleMap/HiSeqV2",
        xena_host="https://tcga.xenahubs.net",
        datapage_url=(
            "https://xenabrowser.net/datapages/"
            "?cohort=TCGA%20Lung%20Adenocarcinoma%20(LUAD)"
        ),
        platform="Illumina HiSeq RNASeqV2 (RSEM, log2 normalized)",
        expected_samples=None,
        notes="RSEM-normalized gene-level expression, PANCAN pipeline.",
    ),
    DatasetEntry(
        cancer_type=CancerType.COAD,
        data_type=DataType.GENE_EXPRESSION_RNASEQ,
        cohort_name="TCGA Colon Adenocarcinoma (COAD)",
        xena_dataset_id="TCGA.COAD.sampleMap/HiSeqV2",
        xena_host="https://tcga.xenahubs.net",
        datapage_url=(
            "https://xenabrowser.net/datapages/"
            "?cohort=TCGA%20Colon%20Cancer%20(COAD)"
        ),
        platform="Illumina HiSeq RNASeqV2 (RSEM, log2 normalized)",
        expected_samples=None,
        notes="RSEM-normalized gene-level expression, PANCAN pipeline.",
    ),
    # TODO: add SOMATIC_MUTATION and CLINICAL entries for each cohort once
    # expression download is working end-to-end. Doing mutation + clinical
    # for all 3 cancers before expression works for even 1 would spread
    # Phase 2 too thin -- get one data type flowing first.
]


def get_datasets_for_cancer(cancer_type: CancerType) -> list[DatasetEntry]:
    """
    Return every registered dataset for a given cancer type.

    Args:
        cancer_type: One of the three in-scope CancerType values.

    Returns:
        A list of DatasetEntry objects (empty list if none registered yet).
    """
    return [entry for entry in DATASET_REGISTRY if entry.cancer_type == cancer_type]


def get_datasets_by_type(data_type: DataType) -> list[DatasetEntry]:
    """
    Return every registered dataset of a given assay type, across all
    cancers. Useful later for "run differential expression on everything
    we have" style batch jobs.
    """
    return [entry for entry in DATASET_REGISTRY if entry.data_type == data_type]
def validate_registry() -> bool:
    """
    Validate dataset registry consistency.

    Checks:
    - Registry is not empty
    - Dataset IDs are unique
    - All entries have required metadata
    """

    if not DATASET_REGISTRY:
        print("Registry validation failed: no datasets found.")
        return False

    dataset_ids = [
        entry.xena_dataset_id
        for entry in DATASET_REGISTRY
    ]

    if len(dataset_ids) != len(set(dataset_ids)):
        print("Registry validation failed: duplicate dataset IDs found.")
        return False

    for entry in DATASET_REGISTRY:
        if not entry.cohort_name:
            print("Registry validation failed: missing cohort name.")
            return False

        if not entry.xena_host:
            print("Registry validation failed: missing Xena host.")
            return False

    print("Dataset registry validation passed.")
    return True



def __init__(self):
        print("Dataset Manager initialized.")

if __name__ == "__main__":

    validate_registry()


    # Quick manual sanity check -- not a substitute for real tests,
    # just lets you eyeball the registry from the command line.
    print(f"Total registered datasets: {len(DATASET_REGISTRY)}\n")
    for cancer in CancerType:
        entries = get_datasets_for_cancer(cancer)
        print(f"{cancer.value}: {len(entries)} dataset(s)")
        for entry in entries:
            print(f"    - {entry.xena_dataset_id} ({entry.data_type.value})")
