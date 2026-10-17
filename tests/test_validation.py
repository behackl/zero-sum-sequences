"""Shared exact-integer rules; component tests check their API boundaries."""

from fractions import Fraction

import pytest

from zero_sum_sequences._validation import _non_negative_integer, _positive_integer

from groups import IndexableInteger


@pytest.mark.parametrize(
    "validate", [_non_negative_integer, _positive_integer], ids=["nonnegative", "positive"]
)
@pytest.mark.parametrize(
    "value", [True, 1.0, Fraction(3, 2), "1", None],
    ids=["boolean", "integral-float", "rational", "string", "none"],
)
def test_integer_validation_requires_the_exact_index_protocol(validate, value):
    with pytest.raises(ValueError, match="integer"):
        validate(value, name="count")


@pytest.mark.parametrize(
    "validate", [_non_negative_integer, _positive_integer], ids=["nonnegative", "positive"]
)
@pytest.mark.parametrize("value", [2, IndexableInteger(2)], ids=["integer", "indexable"])
def test_integer_validation_accepts_exact_indices(validate, value):
    result = validate(value, name="count")
    assert result == 2
    assert type(result) is int


@pytest.mark.parametrize(
    "validate", [_non_negative_integer, _positive_integer], ids=["nonnegative", "positive"]
)
def test_integer_validation_rejects_negative_counts(validate):
    with pytest.raises(ValueError, match="integer"):
        validate(-1, name="count")


def test_zero_is_accepted_only_by_nonnegative_validation():
    assert _non_negative_integer(0, name="count") == 0
    with pytest.raises(ValueError, match="positive integer"):
        _positive_integer(0, name="count")
