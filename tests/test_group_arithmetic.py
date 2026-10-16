"""Lazy group metadata, cached lookup data, and automatic Davenport bounds."""

from math import gcd, lcm

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    FiniteAdditiveGroup,
    davenport_constant,
    exponent,
    zero_sum_constant,
)

from groups import cyclic_group


class ForeignParent:
    """A non-package parent supplying native finite-group metadata."""

    def __init__(self, group, *, invariants=None):
        self.group = group
        self.invariant_factors = invariants
        self.additions = self.iterations = self.exponent_queries = self.order_queries = 0

    def __call__(self, value):
        return self.group(value)

    def zero(self):
        return self.group.zero()

    def is_finite(self):
        return True

    def cardinality(self):
        return len(self.group)

    def exponent(self):
        self.exponent_queries += 1
        return self.group.exponent()

    def element_order(self, term):
        self.order_queries += 1
        return self.group.element_order(term)

    def negate(self, term):
        return self.group.negate(term)

    def invariants(self):
        if self.invariant_factors is None:
            raise NotImplementedError
        return self.invariant_factors

    def add(self, left, right):
        self.additions += 1
        return self.group.add(left, right)

    def __iter__(self):
        self.iterations += 1
        return iter(self.group)


@pytest.mark.parametrize("moduli", [(1,), (6,), (2, 3), (2, 4), (2, 2, 4)])
def test_coordinate_arithmetic_does_not_build_tables(moduli):
    group = FiniteAdditiveGroup.cyclic_product(*moduli)
    assert group._arithmetic_context is None
    assert group.order() == len(group)
    assert group.exponent() == lcm(*moduli)
    for g in group:
        assert group.element_order(g) == lcm(*(n // gcd(x, n) for x, n in zip(g, moduli)))
        assert group.add(g, group.negate(g)) == group.zero()
    assert "addition_table" not in group._group_arithmetic().__dict__


def test_generic_element_orders_and_inverses_are_cached():
    additions = []

    def add(left, right):
        additions.append((left, right))
        return (left + right) % 5

    group = FiniteAdditiveGroup(range(5), zero=0, add=add)
    assert group.element_order(2) == 5
    assert group.negate(2) == 3
    before = len(additions)
    assert group.element_order(2) == 5
    assert group.negate(2) == 3
    assert len(additions) == before
    assert group.exponent() == 5
    before = len(additions)
    assert group.exponent() == 5
    assert len(additions) == before
    assert "addition_table" not in group._group_arithmetic().__dict__


def test_spaces_share_the_group_context_and_lookup_tables():
    group = cyclic_group(3)
    left = AdditiveSequenceSpace(group)
    right = AdditiveSequenceSpace(group)
    assert left._arithmetic_context is right._arithmetic_context is None
    zero_sum_constant(left, {3}, symmetry=False)
    arithmetic = left._group_arithmetic()
    addition = arithmetic.addition_table
    difference = arithmetic.difference_table
    zero_sum_constant(right, {2, 3}, symmetry=False)
    assert right._group_arithmetic() is arithmetic
    assert arithmetic is group._group_arithmetic()
    assert arithmetic.addition_table is addition
    assert arithmetic.difference_table is difference


@pytest.mark.parametrize(
    "moduli, expected",
    [
        ((1,), 1),
        ((6,), 6),
        ((2, 3), 6),  # not already an invariant-factor decomposition
        ((2, 2, 3), 7),  # isomorphic to C2 x C6
        ((4, 6), 13),  # isomorphic to C2 x C12
        ((2, 4, 8), 12),
        ((3, 3, 9), 13),
        ((2, 6, 6), 30),  # mixed-prime rank 3: use a bound, not D*
    ],
)
def test_structural_or_certified_bounds(moduli, expected):
    group = FiniteAdditiveGroup.cyclic_product(*moduli)
    space = AdditiveSequenceSpace(group)
    assert space._arithmetic_context is None
    assert space.davenport_bound == expected
    arithmetic = space._group_arithmetic()
    assert "addition_table" not in arithmetic.__dict__
    assert space.davenport_bound == expected


def test_standard_lower_bound_is_not_mistaken_for_an_upper_bound():
    # This family occurs in Geroldinger-Liebmann-Philipp's counterexamples.
    group = FiniteAdditiveGroup.cyclic_product(2, 6, 6, 6, 6)
    space = AdditiveSequenceSpace(group)
    assert space.davenport_bound > 1 + 1 + 4 * 5
    assert group._group_arithmetic().exact_davenport_constant is None
    assert group._automorphism_group is None


def test_generic_bound_does_not_compute_orders_or_build_tables():
    def forbidden(*args):
        raise AssertionError("bound derivation must not perform group additions")

    group = FiniteAdditiveGroup(range(9), zero=0, add=forbidden)
    space = AdditiveSequenceSpace(group)
    assert space.davenport_bound == 9
    assert space._group_arithmetic()._exponent is None
    assert "addition_table" not in space._group_arithmetic().__dict__


def test_foreign_native_metadata_is_used_without_enumeration_or_tables():
    parent = ForeignParent(FiniteAdditiveGroup.cyclic_product(*([2] * 7)))
    space = AdditiveSequenceSpace(parent)
    assert parent.exponent_queries == 0
    assert space.davenport_bound == 14
    assert exponent(space) == 2
    assert exponent(space) == 2
    assert parent.exponent_queries == 1
    assert parent.iterations == parent.additions == 0
    arithmetic = space._group_arithmetic()
    term = (1,) + (0,) * 6
    assert arithmetic.element_order(term) == arithmetic.element_order(term) == 2
    assert parent.order_queries == 1
    assert arithmetic.negate(term) == term
    assert "addition_table" not in arithmetic.__dict__


def test_new_exponent_metadata_refreshes_a_cached_bound_without_enumeration():
    class DelayedExponent(ForeignParent):
        exponent_available = False

        def exponent(self):
            self.exponent_queries += 1
            if not self.exponent_available:
                raise NotImplementedError
            return 2

    group = FiniteAdditiveGroup(range(8), zero=0, add=lambda a, b: a ^ b)
    parent = DelayedExponent(group)
    automatic = AdditiveSequenceSpace(parent)
    explicit = AdditiveSequenceSpace(parent, davenport_bound=1)
    assert automatic.davenport_bound == 8

    parent.exponent_available = True
    assert exponent(automatic) == 2
    assert automatic.davenport_bound == 6
    assert explicit.davenport_bound == 1
    assert parent.iterations == parent.additions == parent.order_queries == 0

    queries = parent.exponent_queries
    assert exponent(automatic) == 2
    assert automatic.davenport_bound == 6
    assert parent.exponent_queries == queries


def test_foreign_decomposition_uses_a_valid_exact_formula():
    parent = ForeignParent(FiniteAdditiveGroup.cyclic_product(4, 6), invariants=(4, 6))
    space = AdditiveSequenceSpace(parent)
    assert space.davenport_bound == davenport_constant(space) == 13
    assert parent.iterations == parent.additions == 0


def test_foreign_context_is_retained_by_the_space():
    parent = ForeignParent(cyclic_group(3))
    space = AdditiveSequenceSpace(parent)
    zero_sum_constant(space, {3}, symmetry=False)
    before = parent.additions
    zero_sum_constant(space, {2, 3}, symmetry=False)
    assert parent.additions == before
    assert parent.iterations == 1
    assert parent.exponent_queries == 1


def test_exact_computation_tightens_an_automatic_bound_and_is_shared(monkeypatch):
    group = FiniteAdditiveGroup(
        [(a, b) for a in range(2) for b in range(2)],
        zero=(0, 0),
        add=lambda g, h: tuple((a + b) % 2 for a, b in zip(g, h)),
    )
    space = AdditiveSequenceSpace(group)
    other = AdditiveSequenceSpace(group)
    assert space.davenport_bound == other.davenport_bound == 4
    assert davenport_constant(space, symmetry=False) == 3
    assert space.davenport_bound == other.davenport_bound == 3

    def forbidden(*args, **kwargs):
        raise AssertionError("cached constants must not be recomputed")

    monkeypatch.setattr("zero_sum_sequences.invariants.zero_sum_constant", forbidden)
    assert davenport_constant(other, symmetry=False) == 3


def test_explicit_bound_remains_authoritative_and_does_not_initialize_arithmetic():
    group = cyclic_group(3)
    space = AdditiveSequenceSpace(group, davenport_bound=7)
    assert space.davenport_bound == 7
    assert group._arithmetic_context is None
    assert davenport_constant(space, symmetry=False) == 3
    assert space.davenport_bound == 7


def test_automatic_bound_enumerates_complete_catalogues_and_factorizations():
    group = FiniteAdditiveGroup.cyclic_product(2, 4)
    automatic = AdditiveSequenceSpace(group)
    explicit = AdditiveSequenceSpace(group, davenport_bound=5)
    assert [tuple(a) for a in automatic.enumerate_atom_catalogue()] == [
        tuple(a) for a in explicit.enumerate_atom_catalogue()
    ]
    assert automatic([(0, 1)] * 8).length_set() == {2}


def test_unknown_parent_requires_a_bound_only_when_needed():
    class Integers:
        def __call__(self, value):
            return int(value)

        def zero(self):
            return 0

        def is_finite(self):
            return False

    space = AdditiveSequenceSpace(Integers())
    assert space([1, -1]).is_atom()
    assert "davenport_bound=None" in repr(space)
    with pytest.raises(ValueError, match="supply davenport_bound explicitly"):
        _ = space.davenport_bound
    explicit = AdditiveSequenceSpace(Integers(), davenport_bound=2)
    assert explicit([1, -1]).length_set() == {1}


def test_repeated_foreign_elements_are_rejected():
    class RepeatedParent(ForeignParent):
        def __iter__(self):
            return iter([0, 1, 1, 2])

    space = AdditiveSequenceSpace(RepeatedParent(cyclic_group(3)))
    with pytest.raises(ValueError, match="repeated elements"):
        zero_sum_constant(space, {3}, symmetry=False)


def test_missing_inverses_fail_without_an_infinite_order_loop():
    parent = FiniteAdditiveGroup([0, 1], zero=0, add=max)
    with pytest.raises(ValueError, match="finite additive orders"):
        parent.element_order(1)


def test_foreign_parent_without_native_orders_uses_and_caches_addition():
    parent = ForeignParent(cyclic_group(5))
    parent.exponent = parent.element_order = parent.negate = parent.invariants = None
    space = AdditiveSequenceSpace(parent)
    assert space.davenport_bound == 5
    assert parent.additions == parent.iterations == 0
    arithmetic = space._group_arithmetic()
    assert arithmetic.element_order(2) == 5
    assert arithmetic.negate(2) == 3
    assert exponent(space) == 5
    before = parent.additions
    assert exponent(space) == 5
    assert arithmetic.element_order(2) == 5
    assert parent.additions == before
    assert "addition_table" not in arithmetic.__dict__


def test_native_not_implemented_falls_back_to_additive_arithmetic():
    class PartialParent(ForeignParent):
        def element_order(self, term):
            raise NotImplementedError

        def negate(self, term):
            raise NotImplementedError

        def exponent(self):
            raise NotImplementedError

    parent = PartialParent(FiniteAdditiveGroup.cyclic_product(2, 4))
    space = AdditiveSequenceSpace(parent)
    assert space.davenport_bound == 8
    assert parent.additions == parent.iterations == 0
    arithmetic = space._group_arithmetic()
    assert arithmetic.negate((1, 1)) == (1, 3)
    assert arithmetic.element_order((1, 1)) == 4
    assert exponent(space) == 4
    assert zero_sum_constant(space, {4}, symmetry=False).value == 9


def test_metadata_queries_do_not_require_an_iterable_parent():
    class NonIterableParent(ForeignParent):
        __iter__ = None

    space = AdditiveSequenceSpace(NonIterableParent(cyclic_group(5)))
    assert space.davenport_bound == davenport_constant(space) == exponent(space) == 5
    with pytest.raises(TypeError, match="iterable parent"):
        zero_sum_constant(space, {5}, symmetry=False)


def test_native_decomposition_must_match_the_cardinality():
    parent = ForeignParent(cyclic_group(5), invariants=(3,))
    with pytest.raises(ValueError, match="do not match"):
        davenport_constant(AdditiveSequenceSpace(parent), symmetry=False)


def test_custom_addition_takes_precedence_over_native_integer_negation():
    group = FiniteAdditiveGroup(
        range(4), zero=0, add=lambda g, h: g ^ h, coerce=lambda g: int(g) % 4
    )
    assert group.negate(1) == 1  # not (-1) % 4 == 3
    assert group.element_order(1) == group.exponent() == 2
    space = AdditiveSequenceSpace(group)
    assert davenport_constant(space, symmetry=False) == 3
    assert len(space.enumerate_atom_catalogue()) == 4


def test_foreign_custom_addition_does_not_use_the_labels_additive_orders():
    class Label(int):
        def additive_order(self):
            raise AssertionError("the label's native arithmetic is not the group operation")

    class Parent:
        def __call__(self, value):
            return Label(int(value) % 4)

        def zero(self):
            return Label(0)

        def add(self, left, right):
            return Label(left ^ right)

        def is_finite(self):
            return True

        def __iter__(self):
            return iter(map(Label, range(4)))

    space = AdditiveSequenceSpace(Parent())
    arithmetic = space._group_arithmetic()
    assert arithmetic.negate(Label(1)) == Label(1)
    assert arithmetic.element_order(Label(1)) == 2
    assert exponent(space) == 2
    assert davenport_constant(space, symmetry=False) == 3


def test_foreign_native_term_arithmetic_is_used_when_no_parent_operation_exists():
    class Modular(int):
        queries = 0

        def __add__(self, other):
            return Modular((int(self) + int(other)) % 3)

        def __neg__(self):
            return Modular(-int(self) % 3)

        def additive_order(self):
            Modular.queries += 1
            return 3 if self else 1

    class Parent:
        def __call__(self, term):
            return Modular(int(term) % 3)

        def zero(self):
            return Modular(0)

        def is_finite(self):
            return True

        def __iter__(self):
            return iter(map(Modular, range(3)))

    space = AdditiveSequenceSpace(Parent())
    arithmetic = space._group_arithmetic()
    assert arithmetic.negate(Modular(1)) == Modular(2)
    assert arithmetic.element_order(Modular(1)) == 3
    assert arithmetic.element_order(Modular(1)) == 3
    assert Modular.queries == 1
    assert exponent(space) == 3
    assert Modular.queries == 2
    assert "addition_table" not in arithmetic.__dict__
    assert zero_sum_constant(space, {3}, symmetry=False).value == 5
