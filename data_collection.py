"""
BLC Mark - Dataset Manager

Purpose:
    Coordinate dataset management for the BLC Mark platform.

Responsibilities:
    - Read the dataset registry
    - Check dataset availability
    - Coordinate dataset downloads
    - Coordinate dataset validation
    - Report dataset status

Dataset definitions live in dataset_registry.py.
"""

from pathlib import Path

from dataset_registry import DATASET_REGISTRY


class DatasetManager:
    """
    Manage datasets used by BLC Mark.
    """

    def __init__(self):
        """Initialize the Dataset Manager."""
        self.project_root = Path(__file__).resolve().parent.parent
        self.raw_data_dir = self.project_root / "data" / "raw"
        self.dataset_registry = DATASET_REGISTRY

    def list_registered_datasets(self):
        """
        Display all datasets currently registered.
        """
        print("\nRegistered Datasets")
        print("-" * 35)

        for index, dataset in enumerate(self.dataset_registry, start=1):
            print(f"{index}. {dataset.cancer_type.value}")
            print(f"   Data Type : {dataset.data_type.value}")
            print(f"   Platform  : {dataset.platform}")
            print()

    def check_dataset_exists(self, filename):
        """
        Check whether a dataset file exists in the raw data directory.

        Parameters
        ----------
        filename : str
            Dataset filename.

        Returns
        -------
        bool
            True if the file exists, otherwise False.
        """
        dataset_path = self.raw_data_dir / filename
        return dataset_path.exists()


if __name__ == "__main__":
    manager = DatasetManager()

    print("BLC Mark Dataset Manager")
    print("-" * 35)
    print(f"Project Root       : {manager.project_root}")
    print(f"Raw Data Directory : {manager.raw_data_dir}")
    print(f"Registered Datasets: {len(manager.dataset_registry)}")

    manager.list_registered_datasets()

    print("Checking for BRCA dataset...")
    exists = manager.check_dataset_exists("TCGA-BRCA_expression.tsv")
    print(f"Dataset exists: {exists}")

