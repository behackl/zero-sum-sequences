"""Automatic formula recognition has bounded work; exact queries retain support."""

import pytest

from zero_sum_sequences import AdditiveSequenceSpace, davenport_constant
from zero_sum_sequences import group_arithmetic as arithmetic


class MetadataOnly:
    def __init__(self, p=101):
        self.p = p

    def __call__(self, value):
        return int(value)

    def zero(self):
        return 0

    def is_finite(self):
        return True

    def cardinality(self):
        return self.p ** 3

    def exponent(self):
        return self.p

    def invariants(self):
        return (self.p,) * 3

    def __iter__(self):
        raise AssertionError("metadata recognition must not enumerate group elements")


@pytest.mark.parametrize("n, expected", [(1, False), (2, True), (3, True), (4, True), (6, False),
                                         (49, True), (121, True), (143, False), (1024, True)])
def test_exact_prime_power_recognition(n, expected):
    assert arithmetic._is_prime_power(n) is expected


@pytest.mark.parametrize("value", [1_000_000_007, 2**100])
def test_trial_budget_limits_divisibility_tests(value):
    class Counted(int):
        divisions = 0

        def __mod__(self, divisor):
            Counted.divisions += 1
            return int(self) % divisor

        def __floordiv__(self, divisor):
            return Counted(int(self) // divisor)

    assert arithmetic._is_prime_power(Counted(value), max_steps=5) is None
    assert Counted.divisions <= 5


def test_automatic_unknown_is_not_a_false_formula_or_an_exact_search(monkeypatch):
    monkeypatch.setattr(arithmetic, "_AUTOMATIC_TRIAL_LIMIT", 4)

    def forbidden(*args, **kwargs):
        raise AssertionError("formula recognition must not start a constant search")

    monkeypatch.setattr("zero_sum_sequences.invariants.zero_sum_constant", forbidden)
    space = AdditiveSequenceSpace(MetadataOnly())
    expected_bound = 101 * (1 + (101 ** 2 - 1).bit_length())
    assert space.davenport_bound == expected_bound
    assert space._group_arithmetic().exact_davenport_constant is None
    assert "addition_table" not in space._group_arithmetic().__dict__
    # An explicit exact query retries recognition without the automatic budget.
    assert davenport_constant(space) == 301
    assert space.davenport_bound == 301


def test_an_explicit_bound_is_not_replaced_after_full_formula_recognition(monkeypatch):
    monkeypatch.setattr(arithmetic, "_AUTOMATIC_TRIAL_LIMIT", 4)
    space = AdditiveSequenceSpace(MetadataOnly(), davenport_bound=1)
    assert space.davenport_bound == 1
    assert space._arithmetic_context is None
    assert davenport_constant(space) == 301
    assert space.davenport_bound == 1
