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
        ``add(left, right)``. Sums and ``zero()`` are canonicalized through
        the parent, so ``parent(parent.zero())`` is the zero element ``0``
        whatever labels the parent uses.
    davenport_bound:
        Optional positive upper bound for the parent's Davenport constant.
        An explicit value is used as supplied; its sufficiency is not checked
        up front, but factorization queries raise when it provably fails
        (see :class:`FactorizationSolver`).
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
        "_pickle_token",
        "__weakref__",
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
        self._pickle_token = None

    def __reduce__(self):
        """Pickle the parent and explicit bound; caches are rebuilt lazily.

        Unpickling in a process that holds this space (or an earlier copy)
        returns that object, so sequences returned by worker processes
        compare equal to the caller's sequences.
        """

        from ._pickling import _restore, _token

        return _restore, (
            _token(self),
            type(self),
            (self._parent,),
            {"davenport_bound": self._davenport_bound},
        )

    def factorization_cache(
        self, catalogue, *, symmetry=False, canonicalize=None, maxsize: int | None = 1024
    ):
        """A :class:`FactorizationCache` over ``catalogue`` for this space."""

        from .factorization_cache import FactorizationCache

        if catalogue.space is not self:
            raise TypeError("the catalogue belongs to a different space")
        return FactorizationCache(
            catalogue, symmetry=symmetry, canonicalize=canonicalize, maxsize=maxsize
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

        return self.from_terms(self.decode_term(item) for item in data)

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
        """The zero element ``0`` as a canonical, immutable parent element."""

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
        suffices. Otherwise use a known exact value (cached, or a formula from
        :func:`davenport_constant`); else, for order ``n`` and known exponent
        ``e``, return ``min(n, e * (1 + ceil(log2(n / e))))``; else ``n``.

        The logarithmic estimate ``D(G) <= e * (1 + ln(n / e))`` is due to
        van Emde Boas and Kruyswijk (:cite:p:`van-emde-boas-kruyswijk-1969`,
        Theorem 7.1, p. 21; see also :cite:p:`gao-geroldinger-2006`, Theorem
        3.6(1)). The base-two ceiling is a conservative relaxation chosen
        here, not the original estimate.

        Deriving a bound never starts an exhaustive search, although finding
        ``n`` may enumerate the group. Derived bounds are cached and tighten
        when ``D(G)`` or the exponent becomes known. Raise :class:`ValueError`
        if no bound can be derived.
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

    def __call__(self, *terms: Element) -> AdditiveSequence[Element]:
        """Construct the sequence ``g_1 · … · g_l`` of the given terms.

        Each argument is one term, coerced through the parent; ``S()`` is the
        empty sequence ``1``. Use :meth:`from_terms` for an iterable of terms.
        A sequence of this space is returned unchanged.
        """

        if len(terms) == 1 and isinstance(terms[0], AdditiveSequence):
            if terms[0].parent() is self:
                return terms[0]
            raise TypeError(
                "the sequence belongs to a different space; "
                "use from_terms() to coerce its terms"
            )
        try:
            return AdditiveSequence(self, terms)
        except (TypeError, ValueError) as error:
            if len(terms) == 1 and _looks_like_term_iterable(self._parent, terms[0]):
                raise TypeError(
                    "pass terms as separate arguments, S(g, h, ...), "
                    "or use S.from_terms(iterable)"
                ) from error
            raise

    def from_terms(self, terms: Iterable[Element]) -> AdditiveSequence[Element]:
        """Construct the sequence of the terms of an iterable."""

        if isinstance(terms, AdditiveSequence) and terms.parent() is self:
            return terms
        return AdditiveSequence(self, terms)

    def from_multiplicities(
        self, multiplicities: Mapping[Element, int]
    ) -> AdditiveSequence[Element]:
        """Construct a sequence from term multiplicities.

        Keys are coerced through the parent, merging keys with equal images;
        zero counts are ignored.
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
        the parent, with duplicates and the term ``0`` removed. It need not be a
        subgroup: sums and inverses are taken in the ambient group.

        A derived bound gives a complete catalogue on the chosen support; an
        explicit bound must cover all atom lengths there. The atom ``0`` of
        length one is excluded, as in :class:`AtomCatalogue`. Full enumeration with a
        derived bound also records the exact Davenport constant.

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
            # Consume without probing the iterator's optional length hint.
            parent_terms = tuple(term for term in _iter_parent(self._parent))
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
            allowed_terms = self.from_terms(support).support
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
                candidate = self(*prefix, final)
                if candidate.is_atom():
                    atoms.append(candidate)
        catalogue = AtomCatalogue(self, atoms, verify=False)
        if support is None and self._davenport_bound is None:
            # A complete ambient catalogue determines D(G); D = 1 when trivial.
            arithmetic._record_exact_davenport(max(map(len, catalogue), default=1))
        return catalogue

    def __repr__(self) -> str:
        bound = (
            "" if self._davenport_bound is None
            else f", davenport_bound={self._davenport_bound}"
        )
        return f"{type(self).__name__}({self._parent!r}{bound})"


def _looks_like_term_iterable(parent, argument) -> bool:
    """Whether a rejected single argument is an iterable of terms (old call style)."""

    if isinstance(argument, (str, bytes)):
        return False
    if isinstance(argument, Iterator):
        return True
    try:
        items = list(argument)
    except TypeError:
        return False
    try:
        for item in items:
            parent(item)
    except (TypeError, ValueError, ArithmeticError):
        return False
    return True


@functools.total_ordering
class AdditiveSequence(Generic[Element]):
    """A sequence over the parent of an :class:`AdditiveSequenceSpace`.

    By a sequence, we mean a finite unordered sequence of terms, with
    repetition of terms allowed: an element ``g_1 · … · g_l`` of the free
    abelian monoid over the parent. Construct sequences by calling their
    space, ``S(g, h, ...)``. Sequences are immutable and written
    multiplicatively: ``S * T``, ``S ** k``, ``S / T`` for a divisor ``T``,
    and ``-S``. Terms are coerced through the parent, copied when necessary
    to make them immutable, and stored with their multiplicities.
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

    def __reduce__(self):
        # The hash depends on the space object, so it is recomputed on load.
        return _restore_sequence, (self._space, self._items)

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
        """Return the multiplicity ``v_g(S)`` of the term ``g``."""

        coerced = _immutable_term(self._space.base_parent(term))
        return next(
            (count for candidate, count in self._items if candidate == coerced),
            0,
        )

    def sum(self):
        """Return the sum ``σ(S)`` of the terms, in the base parent."""

        space = self._space
        return reduce(space._add_operation(), self, space._canonical_zero())

    def is_zero_sum(self) -> bool:
        """Return whether this is a zero-sum sequence, ``σ(S) = 0``."""

        return bool(self.sum() == self._space._canonical_zero())

    def map_terms(
        self,
        mapping: Callable[[Element], TargetElement],
        *,
        target_space: AdditiveSequenceSpace[TargetElement] | None = None,
    ) -> AdditiveSequence[TargetElement]:
        """Return ``φ(S) = φ(g_1) · … · φ(g_l)`` for the map ``φ = mapping``.

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
        return target_space.from_terms(mapping(term) for term in self)

    def orbit(self, *, group=None):
        """Return this sequence's finite automorphism orbit, sorted.

        ``group`` may be an :class:`AutomorphismGroup`, an
        :class:`AutomorphismAction`, or an iterable of callable generators.
        When omitted, the parent's automorphisms are resolved lazily and
        cached on the sequence space.
        """

        from .orbits import orbit

        return orbit(self, group=group)

    def is_in_same_orbit(
        self,
        other: AdditiveSequence[Element],
        *,
        group=None,
    ) -> bool:
        """Return whether ``other`` is in this sequence's orbit (see :meth:`orbit`)."""

        from .orbits import is_in_same_orbit

        return is_in_same_orbit(self, other, group=group)

    def orbit_witness(
        self,
        other: AdditiveSequence[Element],
        *,
        group=None,
    ):
        """Return a generator word mapping this sequence to ``other``.

        The empty word witnesses equality.  ``None`` is returned when the
        sequences are not in the same orbit.
        """

        from .orbits import orbit_witness

        return orbit_witness(self, other, group=group)

    def subsequences(
        self,
        *,
        nonempty: bool = True,
        proper: bool = False,
        max_length: int | None = None,
    ) -> Iterator[Self]:
        """Yield each distinct subsequence once, in deterministic order."""

        if max_length is not None:
            max_length = _non_negative_integer(max_length, name="maximum length")

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
        """Return whether this is a minimal zero-sum sequence (an atom)."""

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
        """Return a factorization solver for a sequence without the term ``0``.

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
        """Yield each unordered atom factorization, optionally of one length."""

        yield from self.factorization_solver(
            atom_catalogue=atom_catalogue
        ).factorizations(factor_count=factor_count)

    def length_set(self, *, atom_catalogue=None) -> set[int]:
        """Return the length set: the lengths of all factorizations."""

        return self.factorization_solver(atom_catalogue=atom_catalogue).length_set()

    def minimum_factorization_length(self, *, atom_catalogue=None) -> int | None:
        """Return the minimum factorization length, or ``None`` if there is none."""

        return self.factorization_solver(
            atom_catalogue=atom_catalogue
        ).minimum_factorization_length()

    def maximum_factorization_length(self, *, atom_catalogue=None) -> int | None:
        """Return the maximum factorization length, or ``None`` if there is none."""

        return self.factorization_solver(
            atom_catalogue=atom_catalogue
        ).maximum_factorization_length()

    def has_factorization_of_length(
        self,
        factor_count: int,
        *,
        atom_catalogue=None,
    ) -> bool:
        """Return whether a factorization has exactly ``factor_count`` factors."""

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

    def __mul__(self, other: object) -> Self:
        """Return the product ``S · T`` of two sequences."""

        if not isinstance(other, AdditiveSequence):
            return NotImplemented
        self._require_same_space(other)
        counts = dict(self._items)
        for term, count in other._items:
            counts[term] = counts.get(term, 0) + count
        return AdditiveSequence._from_items(self._space, tuple(sorted(counts.items())))

    def __rmul__(self, other: object) -> Self:
        # math.prod() starts with the integer one when no explicit start is given.
        if type(other) is int and other == 1:
            return self
        return NotImplemented

    def __truediv__(self, other: object) -> Self:
        """Return ``S · T^{-1}`` for a divisor ``T`` of ``S``."""

        if not isinstance(other, AdditiveSequence):
            return NotImplemented
        self._require_same_space(other)
        if not other.divides(self):
            raise ValueError(f"{other} does not divide {self}")
        counts = dict(self._items)
        for term, count in other._items:
            counts[term] -= count
        return AdditiveSequence._from_items(
            self._space,
            tuple((term, counts[term]) for term, _ in self._items if counts[term]),
        )

    def __pow__(self, exponent: object, modulo: object = None) -> Self:
        """Return ``S^k`` for a nonnegative integer ``k``."""

        if modulo is not None or isinstance(exponent, bool):
            return NotImplemented
        try:
            index(exponent)
        except TypeError:
            return NotImplemented
        exponent = _non_negative_integer(exponent, name="sequence exponent")
        return AdditiveSequence._from_items(
            self._space,
            tuple((term, count * exponent) for term, count in self._items if exponent),
        )

    def __neg__(self) -> Self:
        """Return ``-S``, the sequence of the negatives of the terms."""

        return self._map_canonical(self._space._group_arithmetic().negate)

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
        return str(self)

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


def _restore_sequence(space, items):
    return AdditiveSequence._from_items(space, items)
