"""Immutable finite sequences over a configured additive parent."""

from __future__ import annotations

import functools
import itertools
from collections import Counter
from collections.abc import Callable, Iterable, Iterator, Mapping
from functools import reduce
from operator import index
from typing import TYPE_CHECKING, Generic, Self, TypeVar

from ._exceptions import _ArithmeticUnavailable, _ParentIterationUnavailable
from ._parent_ops import _addition, _immutable_term, _is_finite, _iter_parent
from ._validation import _non_negative_integer, _positive_integer

if TYPE_CHECKING:
    from .atom_catalogue import AtomCatalogue

Element = TypeVar("Element")
TargetElement = TypeVar("TargetElement")


class AdditiveSequenceSpace(Generic[Element]):
    """Configured constructor for sequences over one additive parent.

    Parameters
    ----------
    parent:
        A callable parent containing every sequence term and providing
        ``zero()``.  Elements must be hashable and mutually orderable.  They
        may implement addition themselves, or the parent may provide
        ``add(left, right)``. Addition returns canonical parent elements;
        ``zero()`` defines the identity regardless of the labels used.
    davenport_bound:
        Optional positive upper bound for the parent's Davenport constant.
        An explicit value is used as supplied; its sufficiency is not checked.
        If omitted, :attr:`davenport_bound` derives a bound when needed,
        without starting an exhaustive constant search.
    """

    __slots__ = (
        "_parent",
        "_addition_op",
        "_zero",
        "_davenport_bound",
        "_arithmetic_context",
        "_automorphism_action",
        "_automorphism_group",
    )

    def __init__(self, parent, *, davenport_bound: int | None = None) -> None:
        if not callable(parent):
            raise TypeError("the additive parent must be callable")
        zero = getattr(parent, "zero", None)
        if not callable(zero):
            raise TypeError("the additive parent must provide zero()")
        self._parent = parent
        self._addition_op = None
        self._zero = None
        self._davenport_bound = (
            _positive_integer(davenport_bound, name="Davenport bound")
            if davenport_bound is not None else None
        )
        self._arithmetic_context = None
        self._automorphism_action = None
        self._automorphism_group = None

    def factorization_cache(self, catalogue, *, group=None, canonicalize=None, maxsize: int | None = 1024):
        """A :class:`FactorizationCache` over ``catalogue`` for this space."""

        from .factorization_cache import FactorizationCache

        if catalogue.space is not self:
            raise TypeError("the catalogue belongs to a different space")
        return FactorizationCache(
            catalogue, group=group, canonicalize=canonicalize, maxsize=maxsize
        )

    def encode_term(self, term: Element) -> object:
        """JSON-compatible encoding of one term, delegated to the parent."""

        encoder = getattr(self._parent, "encode_term", None)
        if callable(encoder):
            return encoder(term)
        from .parents import default_encode_term

        return default_encode_term(term)

    def decode_term(self, data: object) -> Element:
        """Inverse of :meth:`encode_term`, delegated to the parent."""

        decoder = getattr(self._parent, "decode_term", None)
        if callable(decoder):
            return _immutable_term(decoder(data))
        return _immutable_term(self._parent(data))

    def decode(self, data: Iterable[object]) -> AdditiveSequence[Element]:
        """Construct a sequence from a list of encoded terms."""

        return self(self.decode_term(item) for item in data)

    def automorphism_group(self):
        """Return the materialized automorphism group of the additive parent.

        The parent may supply the group itself; otherwise it is the closure
        of the parent's automorphism generators.  The result is cached.
        """

        from .automorphisms import automorphism_group

        return automorphism_group(self)

    @property
    def base_parent(self):
        """Return the configured additive parent of the sequence terms."""

        return self._parent

    def _add_operation(self):
        """Resolve canonical addition once, without finite-group metadata."""

        if self._addition_op is None:
            self._addition_op = _addition(self._parent)
        return self._addition_op

    def _canonical_zero(self):
        """The identity as a canonical, immutable parent element."""

        if self._zero is None:
            self._zero = _immutable_term(self._parent(self._parent.zero()))
        return self._zero

    def _group_arithmetic(self):
        """Return cached arithmetic, shared with package groups and local otherwise."""

        if self._arithmetic_context is None:
            from .parents import FiniteAdditiveGroup
            from .group_arithmetic import _GroupArithmetic

            if isinstance(self._parent, FiniteAdditiveGroup):
                self._arithmetic_context = self._parent._group_arithmetic()
            else:
                self._arithmetic_context = _GroupArithmetic(self._parent)
        return self._arithmetic_context

    @property
    def davenport_bound(self) -> int:
        """Return an upper bound for ``D(G)``, not necessarily its exact value.

        An explicit bound is returned as supplied, without checking that it
        suffices. Otherwise use a cached exact constant or a supported exact
        formula from :func:`davenport_constant`. If neither is available,
        use the group order ``n`` and a known exponent ``e`` to return
        ``min(n, e * (1 + ceil(log2(n / e))))``. With no known exponent,
        return ``n``. Formula recognition has bounded cost; inconclusive
        cases use the conservative estimate. Explicit exact queries can
        spend more time recognizing a formula.

        The logarithmic estimate ``D(G) <= e * (1 + ln(n / e))`` is due to
        van Emde Boas and Kruyswijk (:cite:p:`van-emde-boas-kruyswijk-1969`,
        Theorem 7.1, p. 21; see also :cite:p:`gao-geroldinger-2006`, Theorem
        3.6(1)). The base-two ceiling is a conservative relaxation chosen
        here, not the original estimate.

        Automatic bounds are cached and can tighten when completed searches
        certify ``D(G)`` or new exponent information becomes available. An
        explicit bound is never replaced. Resolution does not enumerate
        element orders or launch an exhaustive constant search; finding the group order may require
        enumerating its elements. Missing capabilities raise :class:`ValueError`
        requesting an explicit bound. Invalid metadata and provider errors
        propagate rather than triggering boundless discovery.
        """

        if self._davenport_bound is not None:
            return self._davenport_bound
        try:
            return self._group_arithmetic().davenport_bound()
        except _ArithmeticUnavailable as error:
            raise ValueError(
                "cannot derive a Davenport bound for this parent; "
                "supply davenport_bound explicitly"
            ) from error

    def _optional_davenport_bound(self) -> int | None:
        """Return an available bound for pruning or metadata, or None."""

        if self._davenport_bound is not None:
            return self._davenport_bound
        try:
            return self._group_arithmetic().davenport_bound()
        except _ArithmeticUnavailable:
            return None

    def __call__(self, terms: Iterable[Element] = ()) -> AdditiveSequence[Element]:
        """Construct a sequence, coercing every term through the parent."""

        if isinstance(terms, AdditiveSequence) and terms.parent() is self:
            return terms
        return AdditiveSequence(self, terms)

    def from_multiplicities(
        self, multiplicities: Mapping[Element, int]
    ) -> AdditiveSequence[Element]:
        """Construct directly from nonnegative counts, merging coerced aliases.

        Zero-count keys are ignored without coercion. Positive-count keys
        are coerced once each; repeated terms are not expanded.
        """

        counts: dict[Element, int] = {}
        for term, count in multiplicities.items():
            count = _non_negative_integer(count, name="multiplicity")
            if count:
                term = _immutable_term(self._parent(term))
                counts[term] = counts.get(term, 0) + count
        try:
            items = tuple(sorted(counts.items()))
        except TypeError as error:
            raise TypeError(
                "additive-sequence terms must be mutually orderable"
            ) from error
        return AdditiveSequence._from_items(self, items)

    def enumerate_atom_catalogue(
        self, *, support: Iterable[Element] | None = None
    ) -> AtomCatalogue[Element]:
        """Enumerate atoms on a nonzero support up to :attr:`davenport_bound`.

        The base parent must be a finite iterable abelian group. By default,
        use every nonzero element. A supplied ``support`` is coerced through
        the parent, with duplicates and zero removed. It need not be a
        subgroup: sums and inverses are taken in the ambient group.

        The bound is resolved when enumeration starts. A derived bound gives
        a complete catalogue on the chosen support; an explicit bound must
        cover all atom lengths there. Such a catalogue need not be complete
        for sequences using terms outside the support. The singleton zero
        is excluded, as in :class:`AtomCatalogue`. Full ambient enumeration
        with an automatic bound also records the exact Davenport constant.

        The search completes each sorted prefix with its uniquely determined
        final term, then checks minimality. It can be expensive even when the
        bound is small.
        """

        from .atom_catalogue import AtomCatalogue

        if not _is_finite(self._parent):
            raise ValueError(
                "atom catalogue enumeration requires a finite parent"
            )

        try:
            parent_terms = tuple(_iter_parent(self._parent))
        except _ParentIterationUnavailable:
            raise TypeError(
                "atom catalogue enumeration requires an iterable parent"
            ) from None

        terms = tuple(
            sorted(
                _immutable_term(self._parent(term))
                for term in parent_terms
            )
        )
        zero = self._canonical_zero()
        operation = self._add_operation()

        if support is None:
            allowed_terms = terms
        else:
            allowed_terms = self(support).support
            if not set(allowed_terms).issubset(terms):
                raise ValueError("support terms must belong to the ambient group")
        nonzero_terms = tuple(term for term in allowed_terms if term != zero)

        # Compute inverses in the ambient group, not in the chosen support.
        arithmetic = self._group_arithmetic()
        completions = {arithmetic.negate(term): term for term in nonzero_terms}
        term_index = {
            term: position for position, term in enumerate(nonzero_terms)
        }

        atoms = []
        for length in range(2, self.davenport_bound + 1):
            for prefix_indices in itertools.combinations_with_replacement(
                range(len(nonzero_terms)), length - 1
            ):
                prefix = tuple(
                    nonzero_terms[position] for position in prefix_indices
                )
                prefix_total = _immutable_term(
                    self._parent(reduce(operation, prefix, zero))
                )
                if prefix_total not in completions:
                    continue
                final = completions[prefix_total]
                if term_index[final] < prefix_indices[-1]:
                    continue
                candidate = self((*prefix, final))
                if candidate.is_atom():
                    atoms.append(candidate)
        catalogue = AtomCatalogue(self, atoms, verify=False)
        if support is None and self._davenport_bound is None:
            # A complete ambient catalogue determines D(G); D = 1 when trivial.
            arithmetic._record_exact_davenport(max(map(len, catalogue), default=1))
        return catalogue

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}({self._parent!r}, "
            f"davenport_bound={self._davenport_bound})"
        )


@functools.total_ordering
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

    @classmethod
    def _from_items(
        cls, space: AdditiveSequenceSpace[Element], items: tuple[tuple[Element, int], ...]
    ) -> AdditiveSequence[Element]:
        """Fast constructor from sorted ``(term, count)`` items of canonical terms."""

        sequence = cls.__new__(cls)
        sequence._space = space
        sequence._items = items
        sequence._length = sum(count for _, count in items)
        sequence._hash = hash((space, items))
        return sequence

    def _map_canonical(self, mapping: Callable[[Element], Element]) -> AdditiveSequence[Element]:
        """Image under a term map whose values are canonical parent elements."""

        counts: dict[Element, int] = {}
        for term, count in self._items:
            image = mapping(term)
            counts[image] = counts.get(image, 0) + count
        return AdditiveSequence._from_items(self._space, tuple(sorted(counts.items())))

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

        space = self._space
        return reduce(space._add_operation(), self, space._canonical_zero())

    def is_zero_sum(self) -> bool:
        """Return whether the sum of the sequence is zero."""

        return bool(self.total() == self._space._canonical_zero())

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
        """Return a factorization solver for a sequence without identity terms.

        A supplied catalogue must contain every atom divisor for complete
        results; no Davenport bound is needed in that case. Otherwise an
        available bound caps atom discovery, with exhaustive discovery on
        this finite sequence as the fallback. See :class:`FactorizationSolver`
        for the search assumptions and cost.
        """

        from .factorization import FactorizationSolver

        return FactorizationSolver(self, atom_catalogue=atom_catalogue)

    def factorizations(
        self,
        *,
        factor_count: int | None = None,
        atom_catalogue=None,
    ):
        """Yield each unordered reduced factorization, optionally of one length."""

        yield from self.factorization_solver(
            atom_catalogue=atom_catalogue
        ).factorizations(factor_count=factor_count)

    def length_set(self, *, atom_catalogue=None) -> set[int]:
        """Return the complete set of attained reduced factorization lengths."""

        return self.factorization_solver(atom_catalogue=atom_catalogue).length_set()

    def minimum_factorization_length(self, *, atom_catalogue=None) -> int | None:
        """Return the minimum reduced factorization length, if one exists."""

        return self.factorization_solver(
            atom_catalogue=atom_catalogue
        ).minimum_factorization_length()

    def maximum_factorization_length(self, *, atom_catalogue=None) -> int | None:
        """Return the maximum reduced factorization length, if one exists."""

        return self.factorization_solver(
            atom_catalogue=atom_catalogue
        ).maximum_factorization_length()

    def has_factorization_of_length(
        self,
        factor_count: int,
        *,
        atom_catalogue=None,
    ) -> bool:
        """Return whether a reduced factorization has the requested length."""

        return self.factorization_solver(
            atom_catalogue=atom_catalogue
        ).has_factorization_of_length(factor_count)

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

    def __radd__(self, other: object) -> Self:
        # sum() starts with the integer zero when no explicit start is given.
        if type(other) is int and other == 0:
            return self
        return NotImplemented

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

    def encode(self) -> list[object]:
        """The sorted term list with multiplicities, terms encoded by the space."""

        return [self._space.encode_term(term) for term in self]

    def __lt__(self, other: object) -> bool:
        """Total order: shorter sequences first, then the sorted term lists."""

        if not isinstance(other, AdditiveSequence):
            return NotImplemented
        self._require_same_space(other)
        return (self._length, tuple(self)) < (other._length, tuple(other))

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

        factors = []
        for term, count in self._items:
            factor = str(term)
            if count > 1:
                if not (factor.startswith("(") and factor.endswith(")")):
                    factor = f"({factor})"
                factor = f"{factor}^{count}"
            factors.append(factor)
        return " · ".join(factors)
