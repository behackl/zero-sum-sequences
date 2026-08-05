from collections import Counter

import pytest
from sage.all import AdditiveAbelianGroup, ZZ, Zmod

from zero_sum_sequences import AdditiveSequenceSpace, AtomCatalogue


def sequence_space():
    return AdditiveSequenceSpace(Zmod(3), davenport_bound=3)


def test_catalogue_support_index_finds_exact_divisors():
    space = sequence_space()
    one, two = (space.base_parent(value) for value in range(1, 3))
    catalogue = AtomCatalogue(
        space,
        [
            space([one, two]),
            space([one, one, one]),
            space([two, two, two]),
        ],
    )
    sequence = space([one, one, one, two, two, two])

    assert list(catalogue.divisors(sequence)) == [
        space([one, two]),
        space([one, one, one]),
        space([two, two, two]),
    ]


def test_catalogue_rejects_a_sequence_from_another_space():
    space = sequence_space()
    other_space = sequence_space()
    one, two = (space.base_parent(value) for value in range(1, 3))
    catalogue = AtomCatalogue(space, [space([one, two])])

    with pytest.raises(TypeError, match="different spaces"):
        list(catalogue.divisors(other_space([other_space.base_parent(0)])))


def test_catalogue_rejects_zero_and_non_atoms():
    space = sequence_space()
    zero = space.base_parent(0)
    one = space.base_parent(1)

    with pytest.raises(ValueError, match="must not contain zero"):
        AtomCatalogue(space, [space([zero])])
    with pytest.raises(ValueError, match="must be atoms"):
        AtomCatalogue(space, [space([one, one])])


@pytest.mark.parametrize(
    ("invariants", "davenport_bound", "expected_counts"),
    [
        ([2, 2], 3, {2: 3, 3: 1}),
        ([2, 4], 5, {2: 5, 3: 9, 4: 16, 5: 8}),
    ],
)
def test_space_enumerates_complete_catalogue_for_finite_groups(
    invariants, davenport_bound, expected_counts
):
    group = AdditiveAbelianGroup(invariants)
    space = AdditiveSequenceSpace(
        group,
        davenport_bound=davenport_bound,
    )

    catalogue = space.enumerate_atom_catalogue()

    assert catalogue.space is space
    assert Counter(map(len, catalogue)) == Counter(expected_counts)
    assert all(atom.is_atom() for atom in catalogue)
    assert all(group.zero() not in atom for atom in catalogue)


def test_space_enumerates_empty_catalogue_for_trivial_group():
    group = AdditiveAbelianGroup([])
    space = AdditiveSequenceSpace(group, davenport_bound=1)

    catalogue = space.enumerate_atom_catalogue()

    assert catalogue.space is space
    assert tuple(catalogue) == ()


def test_enumeration_respects_the_configured_bound():
    space = AdditiveSequenceSpace(Zmod(3), davenport_bound=2)

    assert tuple(space.enumerate_atom_catalogue()) == (space([1, 2]),)


def test_enumeration_rejects_an_infinite_parent():
    space = AdditiveSequenceSpace(ZZ, davenport_bound=3)

    with pytest.raises(ValueError, match="requires a finite parent"):
        space.enumerate_atom_catalogue()


def test_enumeration_rejects_a_noniterable_parent():
    class FiniteNonIterableParent:
        def __call__(self, value):
            return value

        def zero(self):
            return 0

        def is_finite(self):
            return True

    space = AdditiveSequenceSpace(
        FiniteNonIterableParent(),
        davenport_bound=1,
    )

    with pytest.raises(TypeError, match="requires an iterable parent"):
        space.enumerate_atom_catalogue()
