from pathlib import Path

import pytest

from src.configuration_manager import ProjectConfig
from src.exceptions import ConfigurationError


def test_default_configuration_constructs():
    config = ProjectConfig()

    assert config.project_root.is_absolute()
    assert config.data_dir == config.project_root / "data"
    assert config.downloads_dir == config.data_dir / "downloads"
    assert config.processed_dir == config.data_dir / "processed"
    assert config.results_dir == config.project_root / "results"
    assert config.logs_dir == config.project_root / "logs"
    assert config.temp_dir == config.project_root / "tmp"


def test_default_network_configuration():
    config = ProjectConfig()

    assert config.xena_host.startswith(("http://", "https://"))
    assert config.timeout_seconds == 30
    assert config.max_retries == 3
    assert config.chunk_size_bytes == 1024 * 1024
    assert config.backoff_base_seconds == 1.0
    assert config.verify_checksums_by_default is False


def test_configuration_is_frozen():
    config = ProjectConfig()

    with pytest.raises(AttributeError):
        config.timeout_seconds = 60


def test_relative_project_root_is_rejected():
    with pytest.raises(ConfigurationError, match="absolute path"):
        ProjectConfig(project_root=Path("relative/path"))


def test_invalid_data_dir_type_is_rejected():
    with pytest.raises(TypeError, match="data_dir"):
        ProjectConfig(data_dir="data")


def test_zero_timeout_is_rejected():
    with pytest.raises(ConfigurationError, match="timeout_seconds"):
        ProjectConfig(timeout_seconds=0)


def test_negative_max_retries_is_rejected():
    with pytest.raises(ConfigurationError, match="max_retries"):
        ProjectConfig(max_retries=-1)


def test_zero_chunk_size_is_rejected():
    with pytest.raises(ConfigurationError, match="chunk_size_bytes"):
        ProjectConfig(chunk_size_bytes=0)


def test_zero_backoff_is_rejected():
    with pytest.raises(ConfigurationError, match="backoff_base_seconds"):
        ProjectConfig(backoff_base_seconds=0)


def test_invalid_xena_host_is_rejected():
    with pytest.raises(ConfigurationError, match="http:// or https://"):
        ProjectConfig(xena_host="not-a-url")


def test_empty_xena_host_is_rejected():
    with pytest.raises(ConfigurationError, match="must not be empty"):
        ProjectConfig(xena_host="   ")


def test_invalid_xena_host_type_is_rejected():
    with pytest.raises(TypeError, match="xena_host"):
        ProjectConfig(xena_host=123)


def test_invalid_checksum_flag_type_is_rejected():
    with pytest.raises(TypeError, match="verify_checksums_by_default"):
        ProjectConfig(verify_checksums_by_default="yes")


def test_custom_configuration_is_accepted(tmp_path):
    config = ProjectConfig(
        project_root=tmp_path,
        data_dir=tmp_path / "custom_data",
        downloads_dir=tmp_path / "custom_downloads",
        processed_dir=tmp_path / "custom_processed",
        results_dir=tmp_path / "custom_results",
        logs_dir=tmp_path / "custom_logs",
        temp_dir=tmp_path / "custom_tmp",
        xena_host="https://example.org",
        timeout_seconds=60,
        max_retries=5,
        chunk_size_bytes=4096,
        backoff_base_seconds=2.5,
        verify_checksums_by_default=True,
    )

    assert config.project_root == tmp_path
    assert config.data_dir == tmp_path / "custom_data"
    assert config.downloads_dir == tmp_path / "custom_downloads"
    assert config.processed_dir == tmp_path / "custom_processed"
    assert config.results_dir == tmp_path / "custom_results"
    assert config.logs_dir == tmp_path / "custom_logs"
    assert config.temp_dir == tmp_path / "custom_tmp"
    assert config.xena_host == "https://example.org"
    assert config.timeout_seconds == 60
    assert config.max_retries == 5
    assert config.chunk_size_bytes == 4096
    assert config.backoff_base_seconds == 2.5
    assert config.verify_checksums_by_default is True


def test_construction_does_not_create_directories(tmp_path):
    data_dir = tmp_path / "data"
    downloads_dir = data_dir / "downloads"
    processed_dir = data_dir / "processed"
    results_dir = tmp_path / "results"
    logs_dir = tmp_path / "logs"
    temp_dir = tmp_path / "tmp"

    ProjectConfig(
        project_root=tmp_path,
        data_dir=data_dir,
        downloads_dir=downloads_dir,
        processed_dir=processed_dir,
        results_dir=results_dir,
        logs_dir=logs_dir,
        temp_dir=temp_dir,
    )

    assert not data_dir.exists()
    assert not downloads_dir.exists()
    assert not processed_dir.exists()
    assert not results_dir.exists()
    assert not logs_dir.exists()
    assert not temp_dir.exists()