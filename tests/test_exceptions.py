"""Tests for src.differential_expression.exceptions."""

from src.exceptions import BLCMarkError
from src.differential_expression.exceptions import (
    DEError,
    DEValidationError,
    InvalidExpressionMatrixError,
    InvalidMetadataError,
    SampleMismatchError,
    InsufficientReplicationError,
    InvalidConfigurationError,
    UnsupportedMethodError,
    StatisticalMethodError,
    MultipleTestingError,
    QualityControlError,
    ReproducibilityError,
    ResultWritingError,
)


def test_de_error_subclasses_root():
    assert issubclass(DEError, BLCMarkError)


def test_validation_subclasses():
    assert issubclass(DEValidationError, DEError)
    assert issubclass(InvalidExpressionMatrixError, DEValidationError)
    assert issubclass(InvalidMetadataError, DEValidationError)
    assert issubclass(SampleMismatchError, DEValidationError)
    assert issubclass(InsufficientReplicationError, DEValidationError)


def test_other_categories_subclass_de_error_directly():
    for cls in (
        InvalidConfigurationError,
        UnsupportedMethodError,
        StatisticalMethodError,
        MultipleTestingError,
        QualityControlError,
        ReproducibilityError,
        ResultWritingError,
    ):
        assert issubclass(cls, DEError)


def test_catching_root_blc_mark_error_catches_de_failures():
    try:
        raise InvalidExpressionMatrixError("boom")
    except BLCMarkError as error:
        assert isinstance(error, InvalidExpressionMatrixError)
    else:
        raise AssertionError("Expected BLCMarkError to catch the DE exception.")


def test_all_exports_are_distinct_classes():
    classes = [
        DEError,
        DEValidationError,
        InvalidExpressionMatrixError,
        InvalidMetadataError,
        SampleMismatchError,
        InsufficientReplicationError,
        InvalidConfigurationError,
        UnsupportedMethodError,
        StatisticalMethodError,
        MultipleTestingError,
        QualityControlError,
        ReproducibilityError,
        ResultWritingError,
    ]
    assert len(classes) == len(set(classes))