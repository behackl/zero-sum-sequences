"""Structural work limits for catalogue lookup and trusted atom production."""

import pytest

from zero_sum_sequences import AdditiveSequence, AdditiveSequenceSpace, AtomCatalogue, FiniteAdditiveGroup


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


def test_narrow_support_uses_submask_lookup():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(7))
    atoms = [space([(i,)] * 7) for i in range(1, 7)]
    atoms += [space([(i,), (7 - i,)]) for i in range(1, 4)]
    catalogue = AtomCatalogue(space, atoms)
    masks = catalogue._by_support_mask = CountedMasks(catalogue._by_support_mask)
    sequence = space([(1,)] * 7 + [(6,)] * 7)
    assert list(catalogue.divisors(sequence)) == [a for a in catalogue if a.divides(sequence)]
    assert masks.lookups == 3
    assert masks.visits == 0


def test_empty_or_disjoint_support_has_no_divisors():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    catalogue = AtomCatalogue(space, [space([(1,)] * 3)])
    assert list(catalogue.divisors(space())) == []
    assert list(catalogue.divisors(space([(2,)] * 3))) == []
    assert list(AtomCatalogue(space, []).divisors(space([(1,)]))) == []


def test_enumeration_verifies_each_accepted_atom_once(monkeypatch):
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    original = AdditiveSequence.is_atom
    calls = []

    def counted(self):
        calls.append(self)
        return original(self)

    monkeypatch.setattr(AdditiveSequence, "is_atom", counted)
    catalogue = space.enumerate_atom_catalogue()
    assert len(catalogue) == 3
    assert calls == list(catalogue)


def test_subcatalogues_do_not_reverify_atoms(monkeypatch):
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    catalogue = space.enumerate_atom_catalogue()

    def forbidden(*args):
        raise AssertionError("trusted subsets must not repeat atom tests")

    monkeypatch.setattr(AdditiveSequence, "is_atom", forbidden)
    assert catalogue.filter(lambda atom: True).atoms == catalogue.atoms
    assert catalogue.restrict(max_length=2).atoms == (catalogue[0],)


def test_public_catalogue_construction_still_verifies():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    with pytest.raises(ValueError, match="must be atoms"):
        AtomCatalogue(space, [space([(1,)])])
    with pytest.raises(ValueError, match="contain zero"):
        AtomCatalogue(space, [space([(0,)])])
