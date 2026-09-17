"""Reusable indexed collections of atoms."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Generic

from .additive_sequence import AdditiveSequence, AdditiveSequenceSpace, Element

CATALOGUE_SCHEMA = "atom-catalogue-v1"


@dataclass(frozen=True)
class AtomOrbit(Generic[Element]):
    """One automorphism orbit of catalogue atoms.

    ``representative`` is the orbit minimum, ``indices`` the catalogue
    indices of all members (sorted), ``stabilizer_order`` the order of the
    stabilizer of the representative.
    """

    representative: AdditiveSequence[Element]
    indices: tuple[int, ...]
    stabilizer_order: int

    def __len__(self) -> int:
        return len(self.indices)


class AtomCatalogue(Generic[Element]):
    """A collection of nonzero atoms with a stable order and fast divisor queries.

    Atoms are stored in the total order of sequences (length first, then the
    sorted term list); ``index(atom)`` and ``catalogue[i]`` refer to this
    order, which is therefore a stable identifier namespace for any
    serialized data that refers to atoms.

    Factorizations in this package are taken in the reduced block monoid:
    the identity element may be a length-one atom mathematically, but it is
    deliberately not a factorization atom.
    """

    def __init__(
        self,
        space: AdditiveSequenceSpace[Element],
        atoms: Iterable[AdditiveSequence[Element]],
        *,
        annotations: Iterable[tuple[AdditiveSequence[Element], dict]] | None = None,
        verify: bool = True,
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
            if verify:
                if not atom.is_atom():
                    raise ValueError("atom catalogue entries must be atoms")
                if zero in atom:
                    raise ValueError("atom catalogue entries must not contain zero")
            unique_atoms.add(atom)
        self.atoms = tuple(sorted(unique_atoms))
        self._index = {atom: position for position, atom in enumerate(self.atoms)}
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
        self.annotations: dict[AdditiveSequence[Element], dict] = {}
        for atom, annotation in annotations or ():
            if atom not in self._index:
                raise ValueError("annotation for a sequence outside the catalogue")
            self.annotations[atom] = dict(annotation)
        self._orbit_cache: tuple | None = None  # (group, atom -> (representative, word))

    # indexing

    def __iter__(self) -> Iterator[AdditiveSequence[Element]]:
        return iter(self.atoms)

    def __len__(self) -> int:
        return len(self.atoms)

    def __getitem__(self, position: int) -> AdditiveSequence[Element]:
        return self.atoms[position]

    def __contains__(self, atom: object) -> bool:
        return atom in self._index

    def index(self, atom: AdditiveSequence[Element]) -> int:
        """The position of ``atom`` in the catalogue order."""

        try:
            return self._index[atom]
        except KeyError:
            raise ValueError("the sequence is not in the catalogue") from None

    def __repr__(self) -> str:
        return f"AtomCatalogue({len(self)} atoms)"

    # sub-catalogues

    def filter(
        self, predicate: Callable[[AdditiveSequence[Element]], bool]
    ) -> AtomCatalogue[Element]:
        """The catalogue of the atoms satisfying ``predicate``."""

        kept = [atom for atom in self.atoms if predicate(atom)]
        return AtomCatalogue(
            self.space,
            kept,
            annotations=((atom, self.annotations[atom]) for atom in kept if atom in self.annotations),
        )

    def restrict(
        self,
        *,
        support: Iterable[Element] | None = None,
        max_length: int | None = None,
    ) -> AtomCatalogue[Element]:
        """The atoms supported on ``support`` and of length at most ``max_length``."""

        allowed = None if support is None else set(self.space(support).support)
        return self.filter(
            lambda atom: (allowed is None or set(atom.support) <= allowed)
            and (max_length is None or len(atom) <= max_length)
        )

    # divisors

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
            if index is not None:
                target_mask |= 1 << index

        candidates = []
        submask = target_mask
        while submask:
            candidates.extend(self._by_support_mask.get(submask, ()))
            submask = (submask - 1) & target_mask
        for atom in sorted(candidates):
            if all(
                target_counts.get(term, 0) >= count
                for term, count in atom.multiplicities.items()
            ):
                yield atom

    # automorphisms

    def _group(self, group):
        return self.space.automorphism_group() if group is None else group

    def permutation(self, element, *, group=None) -> tuple[int, ...]:
        """An automorphism as a permutation of catalogue indices.

        ``result[i]`` is the index of the image of ``self[i]``.  The
        catalogue must be closed under the automorphism.
        """

        group = self._group(group)
        images = []
        for atom in self.atoms:
            image = group.apply(element, atom)
            if image not in self._index:
                raise ValueError("the catalogue is not closed under the automorphism")
            images.append(self._index[image])
        return tuple(images)

    def _orbit_data(self, group) -> dict:
        """``atom -> (representative, word from the representative)``, cached."""

        if self._orbit_cache is None or self._orbit_cache[0] is not group:
            mapping = {}
            for representative, members in group.orbit_representatives(self.atoms).items():
                words = group.orbit_words(representative)
                if len(words) != len(members):
                    raise ValueError("the catalogue is not closed under the group")
                for atom in members:
                    mapping[atom] = (representative, words[atom])
            self._orbit_cache = (group, mapping)
        return self._orbit_cache[1]

    def orbits(self, *, group=None) -> tuple[AtomOrbit[Element], ...]:
        """The automorphism orbits of the catalogue, sorted by representative."""

        group = self._group(group)
        members: dict = defaultdict(list)
        for atom, (representative, _) in self._orbit_data(group).items():
            members[representative].append(self._index[atom])
        return tuple(
            AtomOrbit(representative, tuple(sorted(indices)), len(group.stabilizer(representative)))
            for representative, indices in sorted(members.items())
        )

    def representative(self, atom: AdditiveSequence[Element], *, group=None):
        """The orbit minimum of a catalogue atom (precomputed for all atoms)."""

        return self._orbit_data(self._group(group))[self._atom(atom)][0]

    def transporter(self, atom: AdditiveSequence[Element], *, group=None):
        """An automorphism mapping ``atom`` to its orbit representative.

        It is the inverse of the breadth-first generator word from the
        representative, so it is deterministic for a fixed generator list
        and cheap for all atoms at once (needs ``compose`` and ``inverse``).
        """

        group = self._group(group)
        _, word = self._orbit_data(group)[self._atom(atom)]
        return group.inverse(group.materialize(word))

    def _atom(self, atom):
        if atom not in self._index:
            raise ValueError("the sequence is not in the catalogue")
        return atom

    # serialization

    def records(self, *, annotate: Callable[[AdditiveSequence[Element]], dict] | None = None):
        """Yield one JSON-compatible record per atom, in catalogue order."""

        for position, atom in enumerate(self.atoms):
            record = {
                "index": position,
                "length": len(atom),
                "sequence": atom.encode(),
            }
            extra = dict(self.annotations.get(atom, {}))
            if annotate is not None:
                extra.update(annotate(atom))
            for key in extra:
                if key in record:
                    raise ValueError(f"annotation key {key!r} is reserved")
            record.update(extra)
            yield record

    def digest(self) -> str:
        """SHA-256 of the atoms (lengths and sequences only, no annotations)."""

        digest = hashlib.sha256()
        for atom in self.atoms:
            digest.update(
                json.dumps(
                    {"length": len(atom), "sequence": atom.encode()},
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode("utf-8")
            )
            digest.update(b"\n")
        return digest.hexdigest()

    def to_jsonl(
        self,
        path: str | Path,
        *,
        annotate: Callable[[AdditiveSequence[Element]], dict] | None = None,
    ) -> int:
        """Write the catalogue as JSON lines (metadata line, then one line per atom)."""

        path = Path(path)
        metadata = {
            "schema": CATALOGUE_SCHEMA,
            "count": len(self.atoms),
            "davenport_bound": self.space.davenport_bound,
            "digest": self.digest(),
        }
        with path.open("w", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(metadata, sort_keys=True, separators=(",", ":")) + "\n")
            for record in self.records(annotate=annotate):
                stream.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
        return len(self.atoms)

    @classmethod
    def from_jsonl(
        cls,
        path: str | Path,
        space: AdditiveSequenceSpace[Element],
        *,
        verify: bool = True,
    ) -> AtomCatalogue[Element]:
        """Read a catalogue written by :meth:`to_jsonl`.

        The order, the indices and the digest are always checked; with
        ``verify=True`` every record is additionally checked to be an atom
        of ``space`` (which is what makes the file trustworthy, at the cost
        of one atom test per record).  Extra record keys are kept as
        annotations.
        """

        path = Path(path)
        with path.open(encoding="utf-8") as stream:
            lines = [line for line in stream if line.strip()]
        if not lines:
            raise ValueError("empty catalogue file")
        metadata = json.loads(lines[0])
        if metadata.get("schema") != CATALOGUE_SCHEMA:
            raise ValueError("not an atom catalogue file")
        atoms = []
        annotations = []
        previous = None
        for position, line in enumerate(lines[1:]):
            record = json.loads(line)
            try:
                atom = space.decode(record["sequence"])
                if record["index"] != position or record["length"] != len(atom):
                    raise ValueError("index or length mismatch")
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"invalid catalogue record {position}") from error
            if previous is not None and not previous < atom:
                raise ValueError(f"catalogue record {position} is out of order")
            previous = atom
            atoms.append(atom)
            extra = {k: v for k, v in record.items() if k not in {"index", "length", "sequence"}}
            if extra:
                annotations.append((atom, extra))
        if len(atoms) != metadata.get("count"):
            raise ValueError("catalogue count does not match its metadata")
        catalogue = cls(space, atoms, annotations=annotations, verify=verify)
        if catalogue.digest() != metadata.get("digest"):
            raise ValueError("catalogue digest does not match its metadata")
        return catalogue
