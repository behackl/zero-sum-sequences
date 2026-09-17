"""A memoizing, optionally group-aware front end to factorization solvers.

A :class:`FactorizationCache` answers length-set questions about many
sequences over one atom catalogue.  It caches complete length sets without
bound (they are small), keeps a bounded number of solver objects so that
follow-up questions (minimum, maximum, witnesses) about a recently seen
sequence do not rebuild the remainder graph, and can be told the
automorphism group of the space: length sets are automorphism invariants,
so queries are then answered on canonical orbit representatives and one
solver run serves a whole orbit.  Witnesses are transported back through
the inverse automorphism.

The cached length sets can be saved to and loaded from JSON lines, bound to
the digest of the catalogue they were computed against.
"""

from __future__ import annotations

import json
from collections import OrderedDict
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Generic

from .additive_sequence import AdditiveSequence, Element
from .atom_catalogue import AtomCatalogue
from .factorization import Factorization, FactorizationSolver

CACHE_SCHEMA = "length-set-cache-v1"


class AnchoredCanonicalizer(Generic[Element]):
    """Cheap partial canonicalization: transport a distinguished atom divisor.

    A full canonical form needs an orbit traversal of the whole sequence,
    which for large automorphism groups costs far more than a factorization.
    Instead, this canonicalizer picks the largest catalogue atom dividing the
    sequence (the *anchor*), applies the precomputed automorphism that maps
    the anchor to its orbit representative, and optionally reduces the result
    under the stabilizer of that representative (a small subgroup, so its
    orbit traversal is cheap).  The result is an automorphic image of the
    input, hence a valid cache key for automorphism-invariant questions.  It
    is not an orbit invariant: sequences of one orbit may anchor differently
    and then receive different keys, which only costs cache efficiency.
    """

    def __init__(self, catalogue: AtomCatalogue[Element], group, *, reduce: bool = True) -> None:
        self.catalogue = catalogue
        self.group = group
        self.reduce = reduce
        self._stabilizers: dict = {}

    def _stabilizer(self, representative):
        stabilizer = self._stabilizers.get(representative)
        if stabilizer is None:
            stabilizer = self._stabilizers[representative] = self.group.stabilizer(representative)
        return stabilizer

    def __call__(self, sequence: AdditiveSequence[Element]):
        """Return ``(image, automorphism)`` with ``image == automorphism(sequence)``."""

        anchor = max(self.catalogue.divisors(sequence), default=None)
        if anchor is None:
            return sequence, self.group.identity
        transporter = self.catalogue.transporter(anchor, group=self.group)
        image = self.group.apply(transporter, sequence)
        if not self.reduce:
            return image, transporter
        stabilizer = self._stabilizer(self.catalogue.representative(anchor, group=self.group))
        form = stabilizer.canonical_form(image, witness="word")
        return form.image, self.group.compose(form.automorphism, transporter)


class FactorizationCache(Generic[Element]):
    """Cached factorization queries over one catalogue.

    Parameters
    ----------
    catalogue:
        The atom catalogue used by every solver.  It must be complete for the
        sequences that will be queried; the cache does not check this.
    group:
        Optional :class:`AutomorphismGroup` of the space.  When given, every
        query is canonicalized first and witnesses are transported back.
    canonicalize:
        How to canonicalize when a group is given: a callable returning
        either an automorphic image of the sequence or ``(image,
        automorphism)``.  The default is the full canonical form (one orbit
        traversal of the sequence per query, worthwhile when solving is the
        expensive part); :class:`AnchoredCanonicalizer` is the cheap
        alternative for large groups.
    maxsize:
        Number of solver objects kept alive (least recently used eviction).
        Length sets are always kept.
    """

    def __init__(
        self,
        catalogue: AtomCatalogue[Element],
        *,
        group=None,
        canonicalize=None,
        maxsize: int | None = 1024,
    ) -> None:
        if not isinstance(catalogue, AtomCatalogue):
            raise TypeError("expected an AtomCatalogue")
        if canonicalize is not None and group is None:
            raise ValueError("canonicalize requires a group")
        if canonicalize is not None and not callable(canonicalize):
            raise TypeError("canonicalize must be callable")
        self._canonicalize = canonicalize
        if maxsize is not None and (isinstance(maxsize, bool) or maxsize < 1):
            raise ValueError("maxsize must be a positive integer or None")
        self.catalogue = catalogue
        self.space = catalogue.space
        self.group = group
        self._maxsize = maxsize
        self._length_sets: dict[AdditiveSequence[Element], frozenset[int]] = {}
        self._solvers: OrderedDict[AdditiveSequence[Element], FactorizationSolver] = OrderedDict()
        self._hits = 0
        self._misses = 0
        self._solver_builds = 0

    # keys

    def _check(self, sequence: AdditiveSequence[Element]) -> None:
        if not isinstance(sequence, AdditiveSequence):
            raise TypeError("expected an additive sequence")
        if sequence.parent() is not self.space:
            raise TypeError("the sequence does not belong to the cache's space")

    def _key(self, sequence):
        """The cache key: the canonical form when a group is configured."""

        if self.group is None:
            return sequence
        if self._canonicalize is not None:
            result = self._canonicalize(sequence)
            return result[0] if isinstance(result, tuple) else result
        return self.group.canonical(sequence)

    def _form(self, sequence):
        """``(key, automorphism)`` with ``key == automorphism(sequence)``."""

        if self.group is None:
            return sequence, None
        if self._canonicalize is not None:
            result = self._canonicalize(sequence)
            if isinstance(result, tuple):
                return result
            return result, self.group.transporter(sequence, result)
        form = self.group.canonical_form(sequence, witness="word")
        return form.image, form.automorphism

    def _transport(self, automorphism, factorization):
        if automorphism is None or factorization is None:
            return factorization
        inverse = self.group.inverse(automorphism)
        return tuple(self.group.apply(inverse, atom) for atom in factorization)

    def product(self, atoms: Iterable[AdditiveSequence[Element]]) -> AdditiveSequence[Element]:
        """The product (sum of multisets) of ``atoms`` in the cache's space."""

        result = self.space(())
        for atom in atoms:
            self._check(atom)
            result = result + atom
        return result

    # solvers

    def solver(self, sequence: AdditiveSequence[Element]) -> FactorizationSolver[Element]:
        """The (cached) solver of the canonical form of ``sequence``."""

        self._check(sequence)
        return self._solver_for_key(self._key(sequence))

    def _solver_for_key(self, key) -> FactorizationSolver[Element]:
        solver = self._solvers.get(key)
        if solver is not None:
            self._solvers.move_to_end(key)
            return solver
        solver = FactorizationSolver(key, atom_catalogue=self.catalogue)
        self._solver_builds += 1
        self._solvers[key] = solver
        if self._maxsize is not None:
            while len(self._solvers) > self._maxsize:
                self._solvers.popitem(last=False)
        return solver

    # queries

    def _cached(self, sequence):
        """``(key, cached length set or None)``."""

        self._check(sequence)
        key = self._key(sequence)
        return key, self._length_sets.get(key)

    def length_set(self, sequence: AdditiveSequence[Element]) -> frozenset[int]:
        """The complete set of factorization lengths of ``sequence``."""

        key, cached = self._cached(sequence)
        if cached is not None:
            self._hits += 1
            return cached
        self._misses += 1
        lengths = frozenset(self._solver_for_key(key).length_set())
        self._length_sets[key] = lengths
        return lengths

    def minimum(self, sequence: AdditiveSequence[Element]) -> int | None:
        """The minimum factorization length, or ``None`` if there is none."""

        key, cached = self._cached(sequence)
        if cached is not None:
            return min(cached, default=None)
        return self._solver_for_key(key).minimum_factorization_length()

    def maximum(self, sequence: AdditiveSequence[Element]) -> int | None:
        """The maximum factorization length, or ``None`` if there is none."""

        key, cached = self._cached(sequence)
        if cached is not None:
            return max(cached, default=None)
        return self._solver_for_key(key).maximum_factorization_length()

    def has_length(self, sequence: AdditiveSequence[Element], length: int) -> bool:
        key, cached = self._cached(sequence)
        if cached is not None:
            return length in cached
        return self._solver_for_key(key).has_factorization_of_length(length)

    def witness(
        self, sequence: AdditiveSequence[Element], length: int
    ) -> Factorization[Element] | None:
        """A factorization of ``sequence`` of the given length, or ``None``."""

        self._check(sequence)
        key, automorphism = self._form(sequence)
        return self._transport(
            automorphism, self._solver_for_key(key).factorization_witness(length)
        )

    def witnesses(self, sequence: AdditiveSequence[Element]) -> dict[int, Factorization[Element]]:
        """One factorization per attained length."""

        self._check(sequence)
        key, automorphism = self._form(sequence)
        found = self._solver_for_key(key).factorization_witnesses()
        return {length: self._transport(automorphism, f) for length, f in found.items()}

    def factorizations(self, sequence: AdditiveSequence[Element]) -> Iterator[Factorization[Element]]:
        """All factorizations of ``sequence`` (transported from the canonical form)."""

        self._check(sequence)
        key, automorphism = self._form(sequence)
        for factorization in self._solver_for_key(key).factorizations():
            yield self._transport(automorphism, factorization)

    # bookkeeping

    def statistics(self) -> dict[str, int]:
        return {
            "cached_length_sets": len(self._length_sets),
            "cached_solvers": len(self._solvers),
            "hits": self._hits,
            "misses": self._misses,
            "solver_builds": self._solver_builds,
        }

    def clear(self) -> None:
        self._length_sets.clear()
        self._solvers.clear()

    # persistence

    def preload(self, records: Iterable[tuple[AdditiveSequence[Element], Iterable[int]]]) -> int:
        """Insert known length sets; keys are canonicalized like queries."""

        count = 0
        for sequence, lengths in records:
            self._check(sequence)
            self._length_sets[self._key(sequence)] = frozenset(int(x) for x in lengths)
            count += 1
        return count

    def records(self) -> Iterator[dict]:
        """The cached length sets as JSON-compatible records, in sequence order."""

        for sequence in sorted(self._length_sets):
            yield {
                "sequence": sequence.encode(),
                "length_set": sorted(self._length_sets[sequence]),
            }

    def to_jsonl(self, path: str | Path) -> int:
        """Write the cached length sets, bound to the catalogue digest."""

        path = Path(path)
        metadata = {
            "schema": CACHE_SCHEMA,
            "catalogue_digest": self.catalogue.digest(),
            "group_order": None if self.group is None else len(self.group),
            "count": len(self._length_sets),
        }
        with path.open("w", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(metadata, sort_keys=True, separators=(",", ":")) + "\n")
            for record in self.records():
                stream.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
        return metadata["count"]

    @classmethod
    def from_jsonl(
        cls,
        path: str | Path,
        catalogue: AtomCatalogue[Element],
        *,
        group=None,
        maxsize: int | None = 1024,
    ) -> FactorizationCache[Element]:
        """Read cached length sets written by :meth:`to_jsonl`.

        The catalogue digest must match; the stored group order, if any,
        must match the supplied group so that keys are canonical in the same
        sense.
        """

        path = Path(path)
        with path.open(encoding="utf-8") as stream:
            lines = [line for line in stream if line.strip()]
        if not lines:
            raise ValueError("empty length-set cache")
        metadata = json.loads(lines[0])
        if metadata.get("schema") != CACHE_SCHEMA:
            raise ValueError("not a length-set cache file")
        if metadata.get("catalogue_digest") != catalogue.digest():
            raise ValueError("the cache was computed against a different catalogue")
        stored_order = metadata.get("group_order")
        if stored_order != (None if group is None else len(group)):
            raise ValueError("the cache was computed with a different automorphism group")
        cache = cls(catalogue, group=group, maxsize=maxsize)
        space = catalogue.space
        count = cache.preload(
            (space.decode(json.loads(line)["sequence"]), json.loads(line)["length_set"])
            for line in lines[1:]
        )
        if count != metadata.get("count"):
            raise ValueError("cache count does not match its metadata")
        return cache
