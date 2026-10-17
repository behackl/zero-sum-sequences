from collections import Counter

import pytest

from zero_sum_sequences import (
    AdditiveSequence,
    AdditiveSequenceSpace,
    AtomCatalogue,
    FiniteAdditiveGroup,
)

from groups import cyclic_group


def sequence_space():
    return AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)


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


def test_catalogue_ignores_target_terms_absent_from_its_index():
    space = sequence_space()
    atom = space([1, 1, 1])
    catalogue = AtomCatalogue(space, [atom])

    assert list(catalogue.divisors(space([1, 1, 1, 2]))) == [atom]


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
    ("modulus", "davenport_bound", "expected_counts"),
    [
        (2, 2, {2: 1}),
        (3, 3, {2: 1, 3: 2}),
        (4, 4, {2: 2, 3: 2, 4: 2}),
    ],
)
def test_space_enumerates_complete_catalogue_for_finite_groups(
    modulus, davenport_bound, expected_counts
):
    group = cyclic_group(modulus)
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
    group = cyclic_group(1)
    space = AdditiveSequenceSpace(group, davenport_bound=1)

    catalogue = space.enumerate_atom_catalogue()

    assert catalogue.space is space
    assert tuple(catalogue) == ()


def test_enumeration_uses_additive_inverses_for_tuple_elements():
    group = FiniteAdditiveGroup.cyclic_product(2, 2)
    space = AdditiveSequenceSpace(group, davenport_bound=3)

    catalogue = space.enumerate_atom_catalogue()

    assert Counter(map(len, catalogue)) == Counter({2: 3, 3: 1})


def test_enumeration_respects_the_configured_bound():
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=2)

    assert tuple(space.enumerate_atom_catalogue()) == (space([1, 2]),)


def test_enumeration_rejects_an_infinite_parent():
    class InfiniteIntegerParent:
        def __call__(self, value):
            return int(value)

        def zero(self):
            return 0

        def is_finite(self):
            return False

    space = AdditiveSequenceSpace(InfiniteIntegerParent(), davenport_bound=3)

    with pytest.raises(ValueError, match="requires a finite parent"):
        space.enumerate_atom_catalogue()


def test_enumeration_does_not_require_an_iterator_length_hint():
    class Iterator:
        def __init__(self):
            self.terms = iter(range(3))

        def __iter__(self):
            return self

        def __next__(self):
            return next(self.terms)

        def __len__(self):
            raise NotImplementedError("iterator length unavailable")

    class Parent:
        def __call__(self, term):
            return int(term) % 3

        def zero(self):
            return 0

        def add(self, left, right):
            return (left + right) % 3

        def is_finite(self):
            return True

        def cardinality(self):
            return 3

        def exponent(self):
            return 3

        def __iter__(self):
            return Iterator()

    space = AdditiveSequenceSpace(Parent())
    assert tuple(space.enumerate_atom_catalogue()) == (
        space([1, 2]), space([1, 1, 1]), space([2, 2, 2]),
    )


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


class CountedMasks(dict):
    def __init__(self, masks):
        super().__init__(masks)
        self.lookups = self.visits = 0

    def get(self, key, default=None):
        self.lookups += 1
        return super().get(key, default)

    def items(self):
        for entry in super().items():
            self.visits += 1
            yield entry


def test_wide_support_scans_only_occupied_masks():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(17))
    catalogue = AtomCatalogue(space, [space([(i,), (17 - i,)]) for i in range(1, 9)])
    masks = catalogue._by_support_mask = CountedMasks(catalogue._by_support_mask)
    sequence = space([(i,) for i in range(1, 17)])
    assert list(catalogue.divisors(sequence)) == list(catalogue)
    assert masks.lookups + masks.visits <= 8


def test_narrow_support_limits_mask_work():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(7))
    atoms = [space([(i,)] * 7) for i in range(1, 7)]
    atoms += [space([(i,), (7 - i,)]) for i in range(1, 4)]
    catalogue = AtomCatalogue(space, atoms)
    masks = catalogue._by_support_mask = CountedMasks(catalogue._by_support_mask)
    sequence = space([(1,)] * 7 + [(6,)] * 7)
    assert list(catalogue.divisors(sequence)) == [a for a in catalogue if a.divides(sequence)]
    assert masks.lookups + masks.visits <= 3


def test_empty_or_disjoint_support_has_no_divisors():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    catalogue = AtomCatalogue(space, [space([(1,)] * 3)])
    assert list(catalogue.divisors(space())) == []
    assert list(catalogue.divisors(space([(2,)] * 3))) == []
    assert list(AtomCatalogue(space, []).divisors(space([(1,)]))) == []


def test_enumeration_does_not_repeat_accepted_atom_tests(monkeypatch):
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    original = AdditiveSequence.is_atom
    calls = []

    def counted(self):
        calls.append(self)
        return original(self)

    monkeypatch.setattr(AdditiveSequence, "is_atom", counted)
    catalogue = space.enumerate_atom_catalogue()
    assert len(catalogue) == 3
    assert all(count <= 1 for count in Counter(calls).values())


def test_subcatalogues_do_not_reverify_atoms(monkeypatch):
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    catalogue = space.enumerate_atom_catalogue()

    def forbidden(*args):
        raise AssertionError("trusted subsets must not repeat atom tests")

    monkeypatch.setattr(AdditiveSequence, "is_atom", forbidden)
    assert catalogue.filter(lambda atom: True).atoms == catalogue.atoms
    assert catalogue.restrict(max_length=2).atoms == (catalogue[0],)
