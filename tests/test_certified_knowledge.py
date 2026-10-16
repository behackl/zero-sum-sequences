"""Only mathematically certified completed work may tighten automatic bounds."""

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace, AtomCatalogue, FiniteAdditiveGroup,
    davenport_constant, exponent, zero_sum_constant,
)


def generic_group(order=4):
    return FiniteAdditiveGroup(range(order), zero=0, add=lambda a, b: a ^ b)


def test_full_automatic_catalogue_records_exact_D_and_shares_it(monkeypatch):
    group = generic_group()
    space, other = AdditiveSequenceSpace(group), AdditiveSequenceSpace(group)
    explicit = AdditiveSequenceSpace(group, davenport_bound=1)
    assert space.davenport_bound == other.davenport_bound == 4
    catalogue = space.enumerate_atom_catalogue()
    assert max(map(len, catalogue)) == 3
    assert group._group_arithmetic().exact_davenport_constant == 3
    assert space.davenport_bound == other.davenport_bound == 3
    assert explicit.davenport_bound == 1

    def forbidden(*args, **kwargs):
        raise AssertionError("the completed catalogue already determined D")

    monkeypatch.setattr("zero_sum_sequences.invariants.zero_sum_constant", forbidden)
    assert davenport_constant(other) == 3


@pytest.mark.parametrize("support", [[], [0], [1], [1, 2], [1, 2, 3]])
def test_restricted_catalogues_do_not_claim_ambient_completeness(support):
    group = generic_group()
    space = AdditiveSequenceSpace(group)
    space.enumerate_atom_catalogue(support=support)
    assert group._group_arithmetic().exact_davenport_constant is None
    assert space.davenport_bound == 4


@pytest.mark.parametrize("bound", [1, 2, 4])
def test_explicitly_bounded_enumeration_does_not_guess_D(bound):
    group = generic_group()
    space = AdditiveSequenceSpace(group, davenport_bound=bound)
    space.enumerate_atom_catalogue()
    assert group._group_arithmetic().exact_davenport_constant is None
    assert space.davenport_bound == bound


@pytest.mark.parametrize("lengths", [{1, 2, 3}, {1, 2, 3, 10**12}])
def test_general_constant_can_certify_D_without_a_second_search(lengths, monkeypatch):
    space = AdditiveSequenceSpace(generic_group())
    result = zero_sum_constant(space, lengths, symmetry=False)
    assert result.value == 3
    assert space._group_arithmetic().exact_davenport_constant == 3

    def forbidden(*args, **kwargs):
        raise AssertionError("the previous constant search already determined D")

    monkeypatch.setattr("zero_sum_sequences.invariants.zero_sum_constant", forbidden)
    assert davenport_constant(space) == space.davenport_bound == 3


@pytest.mark.parametrize("lengths", [{1, 2}, {2}, {3}, {1, 2, 9}])
def test_other_constants_must_not_be_mistaken_for_D(lengths):
    space = AdditiveSequenceSpace(generic_group())
    zero_sum_constant(space, lengths, symmetry=False)
    assert space._group_arithmetic().exact_davenport_constant is None


def test_new_exponent_information_refreshes_automatic_bounds_without_tables():
    group = generic_group(8)
    space, other = AdditiveSequenceSpace(group), AdditiveSequenceSpace(group)
    explicit = AdditiveSequenceSpace(group, davenport_bound=8)
    assert space.davenport_bound == other.davenport_bound == 8
    assert exponent(space) == 2
    assert space.davenport_bound == other.davenport_bound == 6
    assert explicit.davenport_bound == 8
    assert group._group_arithmetic().exact_davenport_constant is None
    assert "addition_table" not in group._group_arithmetic().__dict__


def test_loaded_catalogue_metadata_does_not_certify_completeness(tmp_path):
    source = AdditiveSequenceSpace(generic_group())
    path = tmp_path / "catalogue.jsonl"
    original = source.enumerate_atom_catalogue()
    original.to_jsonl(path)
    target = AdditiveSequenceSpace(generic_group())
    loaded = AtomCatalogue.from_jsonl(path, target)
    assert loaded.digest() == original.digest()
    assert target._arithmetic_context is None
    assert target.davenport_bound == 4
    assert target._group_arithmetic().exact_davenport_constant is None


def test_foreign_parent_exact_knowledge_remains_space_local():
    class Parent:
        def __call__(self, term):
            return int(term)

        def zero(self):
            return 0

        def add(self, a, b):
            return a ^ b

        def is_finite(self):
            return True

        def __iter__(self):
            return iter(range(4))

    parent = Parent()
    first, second = AdditiveSequenceSpace(parent), AdditiveSequenceSpace(parent)
    assert second.davenport_bound == 4
    assert zero_sum_constant(first, {1, 2, 3}, symmetry=False).value == 3
    assert first.davenport_bound == 3
    assert second.davenport_bound == 4
    assert second._group_arithmetic().exact_davenport_constant is None


def test_trivial_full_catalogue_certifies_D_one_without_nonzero_atoms():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup([0], zero=0, add=lambda a, b: 0))
    assert len(space.enumerate_atom_catalogue()) == 0
    assert space._group_arithmetic().exact_davenport_constant == 1
    assert space.davenport_bound == 1
