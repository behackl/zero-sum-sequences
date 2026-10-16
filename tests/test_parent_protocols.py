"""Cross-layer arithmetic contracts for ordinary and foreign parent labels."""

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    FiniteAdditiveGroup,
    davenport_constant,
    zero_sum_constant,
)


class Label(int):
    def is_zero(self):
        return int(self) == 0


class ForeignC3:
    def __call__(self, term):
        return int(term) % 3

    def zero(self):
        return 0

    def is_finite(self):
        return True

    def __iter__(self):
        return iter(range(3))

    def add(self, left, right):
        assert 0 <= left < 3 and 0 <= right < 3
        return left + right  # The parent canonicalizes the result.


@pytest.mark.parametrize("parent", [
    lambda: FiniteAdditiveGroup(
        map(Label, [1, 2, 3]), zero=Label(1),
        coerce=lambda x: Label((int(x) - 1) % 3 + 1),
        add=lambda a, b: Label((int(a) + int(b) - 2) % 3 + 1),
    ),
    ForeignC3,
])
def test_parent_arithmetic_agrees_across_predicates_catalogues_and_constants(parent):
    parent = parent()
    space = AdditiveSequenceSpace(parent)
    nonzero = tuple(t for t in parent if t != parent.zero())
    a, b = nonzero
    pair, pure_a, pure_b = space([a, b]), space([a] * 3), space([b] * 3)
    for sequence in (pair, pure_a, pure_b):
        assert sequence.total() == parent.zero()
        assert sequence.is_zero_sum()
        assert sequence.is_atom()
    assert not space([a]).is_zero_sum()
    assert not (pair + pair).is_atom()
    assert space().is_zero_sum()
    assert space._arithmetic_context is None
    catalogue = space.enumerate_atom_catalogue()
    assert set(catalogue) == {pair, pure_a, pure_b}
    sequence = pure_a + pure_b
    assert sequence.length_set() == sequence.length_set(atom_catalogue=catalogue) == {2, 3}
    assert davenport_constant(space, symmetry=False) == 3
    assert zero_sum_constant(space, {1, 2, 3}, symmetry=False).value == 3


def test_native_group_totals_do_not_coerce_an_already_canonical_result():
    coercions = []

    def coerce(term):
        coercions.append(term)
        return int(term) % 3

    parent = FiniteAdditiveGroup(range(3), zero=0, add=lambda a, b: (a + b) % 3, coerce=coerce)
    space = AdditiveSequenceSpace(parent, davenport_bound=3)
    sequence = space([1, 2])
    coercions.clear()
    assert sequence.total() == 0
    assert len(coercions) == 3 * len(sequence)  # FiniteAdditiveGroup.add's own checks only.
    assert parent._arithmetic_context is space._arithmetic_context is None


def test_totals_do_not_require_a_finite_parent():
    class Integers:
        def __call__(self, term):
            return int(term)

        def zero(self):
            return 0

    space = AdditiveSequenceSpace(Integers())
    assert space([-2, 1, 1]).is_atom()
    assert space._arithmetic_context is None
