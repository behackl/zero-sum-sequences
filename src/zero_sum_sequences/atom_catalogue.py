"""Reusable indexed collections of atoms."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Iterator
from typing import Generic

from .additive_sequence import AdditiveSequence, AdditiveSequenceSpace, Element


class AtomCatalogue(Generic[Element]):
    """A collection of nonzero atoms indexed for fast divisor queries.

    Factorizations in this package are taken in the reduced block monoid:
    the identity element may be a length-one atom mathematically, but it is
    deliberately not a factorization atom.
    """

    def __init__(
        self,
        space: AdditiveSequenceSpace[Element],
        atoms: Iterable[AdditiveSequence[Element]],
    ) -> None:
        if not isinstance(space, AdditiveSequenceSpace):
            raise TypeError("expected an AdditiveSequenceSpace")
        self.space = space
        zero = space.base_parent.zero()
        unique_atoms = set()
        for atom in atoms:
            if not isinstance(atom, AdditiveSequence):
                raise TypeError("atom catalogue entries must be additive sequences")
            if atom.parent() is not space:
                raise TypeError("atom catalogue entries use a different space")
            if not atom.is_atom():
                raise ValueError("atom catalogue entries must be atoms")
            if zero in atom:
                raise ValueError("atom catalogue entries must not contain zero")
            unique_atoms.add(atom)
        self.atoms = tuple(
            sorted(unique_atoms, key=lambda atom: (len(atom), tuple(atom)))
        )
        self.terms = tuple(
            sorted({term for atom in self.atoms for term in atom.support})
        )
        self._term_index = {
            term: index for index, term in enumerate(self.terms)
        }
        by_support_mask: dict[int, list[AdditiveSequence[Element]]] = defaultdict(list)
        for atom in self.atoms:
            support_mask = 0
            for term in atom.support:
                support_mask |= 1 << self._term_index[term]
            by_support_mask[support_mask].append(atom)
        self._by_support_mask = {
            mask: tuple(mask_atoms)
            for mask, mask_atoms in by_support_mask.items()
        }

    def __iter__(self) -> Iterator[AdditiveSequence[Element]]:
        return iter(self.atoms)

    def __len__(self) -> int:
        return len(self.atoms)

    def divisors(
        self,
        sequence: AdditiveSequence[Element],
    ) -> Iterator[AdditiveSequence[Element]]:
        """Yield catalogue atoms dividing ``sequence`` in canonical order."""

        if sequence.parent() is not self.space:
            raise TypeError("atom catalogue and sequence use different spaces")

        target_counts = sequence.multiplicities
        target_mask = 0
        for term in sequence.support:
            index = self._term_index.get(term)
            if index is None:
                return
            target_mask |= 1 << index

        candidates = []
        submask = target_mask
        while submask:
            candidates.extend(self._by_support_mask.get(submask, ()))
            submask = (submask - 1) & target_mask
        for atom in sorted(candidates, key=lambda item: (len(item), tuple(item))):
            if all(
                target_counts.get(term, 0) >= count
                for term, count in atom.multiplicities.items()
            ):
                yield atom
