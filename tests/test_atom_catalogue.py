import pytest
from sage.all import Zmod

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
