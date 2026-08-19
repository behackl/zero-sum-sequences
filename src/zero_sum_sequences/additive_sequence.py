"""Immutable finite sequences over a configured additive parent."""

from __future__ import annotations

import itertools
from collections import Counter
from collections.abc import Callable, Iterable, Iterator, Mapping
from copy import copy
from functools import reduce
from operator import add as default_add
from operator import index
from typing import TYPE_CHECKING, Generic, Self, TypeVar

if TYPE_CHECKING:
    from .atom_catalogue import AtomCatalogue

Element = TypeVar("Element")
TargetElement = TypeVar("TargetElement")


def _non_negative_integer(value: object, *, name: str) -> int:
    """Return an exact indexable integer, rejecting booleans."""

    if isinstance(value, bool):
        raise ValueError(f"{name} must be a non-negative integer")
    try:
        integer = index(value)
    except TypeError:
        raise ValueError(f"{name} must be a non-negative integer") from None
    if integer < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return integer


def _positive_integer(value: object, *, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive integer")
    try:
        integer = index(value)
    except TypeError:
        raise ValueError(f"{name} must be a positive integer") from None
    if integer < 1:
        raise ValueError(f"{name} must be a positive integer")
    return integer


def _immutable_term(term: Element) -> Element:
    """Return a hashable term, copying supported mutable elements if necessary."""

    try:
        hash(term)
    except TypeError:
        immutable_copy = copy(term)
        set_immutable = getattr(immutable_copy, "set_immutable", None)
        if not callable(set_immutable):
            raise TypeError("additive-sequence terms must be hashable") from None
        set_immutable()
        hash(immutable_copy)
        return immutable_copy
    return term


class AdditiveSequenceSpace(Generic[Element]):
    """Configured constructor for sequences over one additive parent.

    Parameters
    ----------
    parent:
        A callable parent containing every sequence term and providing
        ``zero()``.  Elements must be hashable and mutually orderable.  They
        may implement addition themselves, or the parent may provide
        ``add(left, right)``.
    davenport_bound:
        A positive upper bound for the parent's Davenport constant.
    """

    __slots__ = ("_parent", "_davenport_bound", "_automorphism_action")

    def __init__(self, parent, *, davenport_bound: int) -> None:
        if not callable(parent):
            raise TypeError("the additive parent must be callable")
        zero = getattr(parent, "zero", None)
        if not callable(zero):
            raise TypeError("the additive parent must provide zero()")
        self._parent = parent
        self._davenport_bound = _positive_integer(
            davenport_bound,
            name="Davenport bound",
        )
        self._automorphism_action = None

    @property
    def base_parent(self):
        """Return the configured additive parent of the sequence terms."""

        return self._parent

    @property
    def davenport_bound(self) -> int:
        """Return the configured upper bound for the Davenport constant."""

        return self._davenport_bound

    def __call__(self, terms: Iterable[Element] = ()) -> AdditiveSequence[Element]:
        """Construct a sequence, coercing every term through the parent."""

        if isinstance(terms, AdditiveSequence) and terms.parent() is self:
            return terms
        return AdditiveSequence(self, terms)

    def from_multiplicities(
        self, multiplicities: Mapping[Element, int]
    ) -> AdditiveSequence[Element]:
        """Construct a sequence from non-negative term multiplicities."""

        terms: list[Element] = []
        for term, count in multiplicities.items():
            try:
                count = _non_negative_integer(count, name="multiplicities")
            except ValueError:
                raise ValueError(
                    "multiplicities must be non-negative integers"
                ) from None
            terms.extend(itertools.repeat(term, count))
        return self(terms)

    def enumerate_atom_catalogue(self) -> AtomCatalogue[Element]:
        """Exhaustively enumerate reduced atoms through the configured bound.

        The base parent must be a finite iterable additive group.  For each
        candidate length, sorted prefixes are completed by their uniquely
        determined final term.  The identity singleton is omitted, matching
        the reduced-factorization convention of :class:`AtomCatalogue`.  The
        result is complete only when ``davenport_bound`` is a valid upper
        bound for the base parent.
        """

        from .atom_catalogue import AtomCatalogue

        is_finite = getattr(self._parent, "is_finite", None)
        if not callable(is_finite) or not is_finite():
            raise ValueError(
                "atom catalogue enumeration requires a finite parent"
            )

        try:
            parent_terms = tuple(self._parent)
        except TypeError:
            raise TypeError(
                "atom catalogue enumeration requires an iterable parent"
            ) from None

        terms = tuple(
            sorted(
                _immutable_term(self._parent(term))
                for term in parent_terms
            )
        )
        zero = _immutable_term(self._parent(self._parent.zero()))
        operation = getattr(self._parent, "add", default_add)
        if not callable(operation):
            operation = default_add

        inverse = {}
        for term in terms:
            for candidate in terms:
                total = _immutable_term(
                    self._parent(operation(term, candidate))
                )
                if total == zero:
                    inverse[term] = candidate
                    break
            else:
                raise ValueError(
                    "atom catalogue enumeration requires additive inverses"
                )

        nonzero_terms = tuple(term for term in terms if term != zero)
        term_index = {
            term: position for position, term in enumerate(nonzero_terms)
        }

        atoms = []
        for length in range(2, self._davenport_bound + 1):
            for prefix_indices in itertools.combinations_with_replacement(
                range(len(nonzero_terms)), length - 1
            ):
                prefix = tuple(
                    nonzero_terms[position] for position in prefix_indices
                )
                prefix_total = _immutable_term(
                    self._parent(reduce(operation, prefix, zero))
                )
                final = inverse[prefix_total]
                final_position = term_index.get(final)
                if (
                    final_position is None
                    or final_position < prefix_indices[-1]
                ):
                    continue
                candidate = self((*prefix, final))
                if candidate.is_atom():
                    atoms.append(candidate)
        return AtomCatalogue(self, atoms)

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}({self._parent!r}, "
            f"davenport_bound={self._davenport_bound})"
        )


class AdditiveSequence(Generic[Element]):
    """An immutable finite multiset in an :class:`AdditiveSequenceSpace`.

    Construct sequences by calling their configured space. Terms are coerced
    through its parent, copied when necessary to make them immutable, and
    stored as a sorted multiplicity table.
    """

    __slots__ = ("_space", "_items", "_length", "_hash")

    def __init__(
        self,
        space: AdditiveSequenceSpace[Element],
        terms: Iterable[Element] = (),
    ) -> None:
        if not isinstance(space, AdditiveSequenceSpace):
            raise TypeError("construct sequences through AdditiveSequenceSpace")
        self._space = space
        counts = Counter(
            _immutable_term(space.base_parent(term))
            for term in terms
        )
        try:
            items = tuple(sorted(counts.items()))
        except TypeError as error:
            raise TypeError(
                "additive-sequence terms must be mutually orderable"
            ) from error
        self._items: tuple[tuple[Element, int], ...] = items
        self._length = sum(count for _, count in items)
        self._hash = hash((space, items))

    def parent(self) -> AdditiveSequenceSpace[Element]:
        """Return the configured sequence space."""

        return self._space

    @property
    def support(self) -> tuple[Element, ...]:
        """The distinct terms, in canonical order."""

        return tuple(term for term, _ in self._items)

    @property
    def multiplicities(self) -> dict[Element, int]:
        """A copy of the term-to-multiplicity mapping."""

        return dict(self._items)

    def multiplicity(self, term: Element) -> int:
        """Return how often ``term`` occurs."""

        coerced = _immutable_term(self._space.base_parent(term))
        return next(
            (count for candidate, count in self._items if candidate == coerced),
            0,
        )

    def total(self):
        """Return the additive sum of the terms in the base parent."""

        parent = self._space.base_parent
        operation = getattr(parent, "add", default_add)
        if not callable(operation):
            operation = default_add
        return reduce(operation, self, parent.zero())

    def is_zero_sum(self) -> bool:
        """Return whether the sum of the sequence is zero."""

        total = self.total()
        is_zero = getattr(total, "is_zero", None)
        if callable(is_zero):
            return bool(is_zero())
        return total == self._space.base_parent.zero()

    def map_terms(
        self,
        mapping: Callable[[Element], TargetElement],
        *,
        target_space: AdditiveSequenceSpace[TargetElement] | None = None,
    ) -> AdditiveSequence[TargetElement]:
        """Apply ``mapping`` to every term and construct the image sequence.

        By default the image is constructed in this sequence's space.  A
        different ``target_space`` may be supplied when the mapping takes
        terms to another additive parent.  The result is canonicalized by the
        target space, so equal images are combined with their multiplicities.
        """

        if not callable(mapping):
            raise TypeError("mapping must be callable")
        if target_space is None:
            target_space = self._space
        elif not isinstance(target_space, AdditiveSequenceSpace):
            raise TypeError("target_space must be an AdditiveSequenceSpace")
        return target_space(mapping(term) for term in self)

    def orbit(self, *, action=None):
        """Return this sequence's finite automorphism orbit.

        ``action`` may be an :class:`AutomorphismAction` or an iterable of
        callable term maps.  When omitted, the action is resolved lazily from
        the additive parent and cached on the sequence space.
        """

        from .orbits import orbit

        return orbit(self, action=action)

    def is_in_same_orbit(
        self,
        other: AdditiveSequence[Element],
        *,
        action=None,
    ) -> bool:
        """Return whether ``other`` is in this sequence's action orbit."""

        from .orbits import is_in_same_orbit

        return is_in_same_orbit(self, other, action=action)

    def orbit_witness(
        self,
        other: AdditiveSequence[Element],
        *,
        action=None,
    ):
        """Return a generator word mapping this sequence to ``other``.

        The empty word witnesses equality.  ``None`` is returned when the
        sequences are not in the same orbit.
        """

        from .orbits import orbit_witness

        return orbit_witness(self, other, action=action)

    def subsequences(
        self,
        *,
        nonempty: bool = True,
        proper: bool = False,
        max_length: int | None = None,
    ) -> Iterator[Self]:
        """Yield each distinct subsequence once, in deterministic order."""

        if max_length is not None and max_length < 0:
            raise ValueError("maximum length must be non-negative")

        if max_length is None:
            count_vectors = itertools.product(
                *(range(count + 1) for _, count in self._items)
            )
        else:

            def bounded_count_vectors(
                index: int, remaining: int, prefix: tuple[int, ...] = ()
            ) -> Iterator[tuple[int, ...]]:
                if index == len(self._items):
                    yield prefix
                    return
                maximum_count = min(self._items[index][1], remaining)
                for count in range(maximum_count + 1):
                    yield from bounded_count_vectors(
                        index + 1,
                        remaining - count,
                        (*prefix, count),
                    )

            count_vectors = bounded_count_vectors(0, max_length)

        for chosen_counts in count_vectors:
            chosen_length = sum(chosen_counts)
            if nonempty and chosen_length == 0:
                continue
            if proper and chosen_length == len(self):
                continue
            yield self._space.from_multiplicities(
                {
                    term: count
                    for (term, _), count in zip(self._items, chosen_counts)
                    if count
                }
            )

    def is_atom(self) -> bool:
        """Return whether this is a nonempty minimal zero-sum sequence."""

        if not self or not self.is_zero_sum():
            return False
        return not any(
            subsequence.is_zero_sum()
            for subsequence in self.subsequences(
                nonempty=True,
                proper=True,
                max_length=len(self) // 2,
            )
        )

    def factorization_solver(self, *, atom_catalogue=None):
        """Return a solver for reduced factorizations into nonzero atoms.

        The sequence must not contain the identity term.  The identity remains
        a length-one atom in the full block monoid, but is excluded from the
        reduced factorization convention used by this package.
        """

        from .factorization import FactorizationSolver

        return FactorizationSolver(self, atom_catalogue=atom_catalogue)

    def factorizations(self, *, atom_catalogue=None):
        """Yield every unordered reduced factorization exactly once."""

        yield from self.factorization_solver(
            atom_catalogue=atom_catalogue
        ).factorizations()

    def length_set(self, *, atom_catalogue=None) -> set[int]:
        """Return the complete set of attained reduced factorization lengths."""

        return self.factorization_solver(atom_catalogue=atom_catalogue).length_set()

    def factorization_witnesses(self, *, atom_catalogue=None):
        """Return one deterministic factorization per attained length."""

        return self.factorization_solver(
            atom_catalogue=atom_catalogue
        ).factorization_witnesses()

    def factorization_digraph(self, *, atom_catalogue=None):
        """Return the memoized remainder DAG as a NetworkX directed graph."""

        return self.factorization_solver(atom_catalogue=atom_catalogue).digraph()

    def divides(self, other: AdditiveSequence[Element]) -> bool:
        """Return whether this sequence is a subsequence of ``other``."""

        self._require_same_space(other)
        other_counts = other.multiplicities
        return all(other_counts.get(term, 0) >= count for term, count in self._items)

    def _require_same_space(self, other: object) -> AdditiveSequence[Element]:
        if not isinstance(other, AdditiveSequence):
            raise TypeError("expected an additive sequence")
        if self._space is not other._space:
            raise TypeError("additive sequences belong to different spaces")
        return other

    def __iter__(self) -> Iterator[Element]:
        for term, count in self._items:
            yield from itertools.repeat(term, count)

    def __len__(self) -> int:
        return self._length

    def __bool__(self) -> bool:
        return self._length != 0

    def __contains__(self, term: object) -> bool:
        try:
            coerced = _immutable_term(self._space.base_parent(term))
        except (TypeError, ValueError):
            return False
        return any(candidate == coerced for candidate, _ in self._items)

    def __add__(self, other: object) -> Self:
        if not isinstance(other, AdditiveSequence):
            return NotImplemented
        self._require_same_space(other)
        counts = self.multiplicities
        for term, count in other._items:
            counts[term] = counts.get(term, 0) + count
        return self._space.from_multiplicities(counts)

    def __sub__(self, other: object) -> Self:
        if not isinstance(other, AdditiveSequence):
            return NotImplemented
        self._require_same_space(other)
        if not other.divides(self):
            raise ValueError(f"{other} is not a subsequence of {self}")
        counts = self.multiplicities
        for term, count in other._items:
            counts[term] -= count
        return self._space.from_multiplicities(counts)

    def __mul__(self, repetitions: object) -> Self:
        try:
            repetitions = index(repetitions)
        except TypeError:
            return NotImplemented
        if repetitions < 0:
            raise ValueError("sequence repetitions must be non-negative")
        return self._space.from_multiplicities(
            {term: count * repetitions for term, count in self._items}
        )

    def __rmul__(self, repetitions: object) -> Self:
        return self * repetitions

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AdditiveSequence):
            return NotImplemented
        return self._space is other._space and self._items == other._items

    def __hash__(self) -> int:
        return self._hash

    def __repr__(self) -> str:
        return f"{self._space!r}({list(self)!r})"

    def __str__(self) -> str:
        if not self:
            return "1"
        return " · ".join(
            str(term) if count == 1 else f"({term})^{count}"
            for term, count in self._items
        )
