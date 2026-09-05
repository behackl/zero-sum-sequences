"""Catalogue generation on a subset uses arithmetic in the ambient group."""

from collections import Counter
from itertools import combinations, combinations_with_replacement

import pytest

from zero_sum_sequences import AdditiveSequence, AdditiveSequenceSpace, FiniteAdditiveGroup

from groups import cyclic_group


@pytest.mark.parametrize("modulus", [3, 4, 5])
def test_every_support_matches_filtering_the_complete_catalogue(modulus):
    space = AdditiveSequenceSpace(cyclic_group(modulus), davenport_bound=modulus)
    full = space.enumerate_atom_catalogue()
    terms = tuple(space.base_parent)
    for size in range(len(terms) + 1):
        for support in combinations(terms, size):
            catalogue = space.enumerate_atom_catalogue(support=iter(support))
            expected = tuple(atom for atom in full if set(atom.support) <= set(support))
            assert catalogue.space is space
            assert tuple(catalogue) == expected


def test_non_subgroup_support_uses_ambient_inverses():
    space = AdditiveSequenceSpace(cyclic_group(5), davenport_bound=5)
    support = (1, 2)  # Not closed under sums or negation.
    catalogue = space.enumerate_atom_catalogue(support=support)
    expected = {
        space(terms)
        for length in range(2, 6)
        for terms in combinations_with_replacement(support, length)
        if space(terms).is_atom()
    }

    assert set(catalogue) == expected
    assert space([1, 2, 2]) in catalogue
    assert space([1] * 5) in catalogue
    assert space([2] * 5) in catalogue


def test_support_is_coerced_deduplicated_and_zero_is_ignored():
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    support = (term for term in [0, 4, 1, 4, 3])

    assert tuple(space.enumerate_atom_catalogue(support=support)) == (space([1] * 3),)
    assert tuple(space.enumerate_atom_catalogue(support=space([1, 1]))) == (
        space([1] * 3),
    )


@pytest.mark.parametrize("support", [[], [0], [0, 0]])
def test_empty_or_zero_only_support_has_no_reduced_atoms(support):
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)

    assert tuple(space.enumerate_atom_catalogue(support=support)) == ()


def test_support_generation_does_not_enumerate_outside_atoms(monkeypatch):
    space = AdditiveSequenceSpace(cyclic_group(5), davenport_bound=5)
    original = AdditiveSequence.is_atom

    def check_supported_candidate(sequence):
        assert set(sequence.support) <= {1}
        return original(sequence)

    monkeypatch.setattr(AdditiveSequence, "is_atom", check_supported_candidate)

    assert tuple(space.enumerate_atom_catalogue(support=[1])) == (space([1] * 5),)


def test_four_point_support_catalogue_and_factorizations():
    group = FiniteAdditiveGroup.cyclic_product(3, 3, 3)
    space = AdditiveSequenceSpace(group, davenport_bound=7)
    e1, e2, e3, e4 = (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)
    support = [e1, e2, e3, e4]
    # Lists are coerced to the ambient group's tuple representation.
    catalogue = space.enumerate_atom_catalogue(support=[list(g) for g in support])
    long_atom = space([e1, e1, e2, e2, e3, e3, e4])
    other_atom = space([e1, e2, e3, e4, e4])

    assert Counter(map(len, catalogue)) == {3: 4, 5: 1, 7: 1}
    assert set(catalogue) == {space([g] * 3) for g in support} | {
        long_atom, other_atom
    }
    assert (3 * long_atom).length_set(atom_catalogue=catalogue) == {3, 5, 7}
    assert (long_atom + other_atom).length_set(atom_catalogue=catalogue) == {2, 4}


def test_support_generation_respects_the_bound():
    space = AdditiveSequenceSpace(cyclic_group(5), davenport_bound=4)

    assert tuple(space.enumerate_atom_catalogue(support=[1])) == ()


def test_invalid_support_is_rejected():
    group = FiniteAdditiveGroup(range(3), zero=0, add=lambda a, b: (a + b) % 3)
    space = AdditiveSequenceSpace(group, davenport_bound=3)

    with pytest.raises(ValueError, match="not an element"):
        space.enumerate_atom_catalogue(support=[3])
    with pytest.raises(TypeError):
        space.enumerate_atom_catalogue(support=42)
