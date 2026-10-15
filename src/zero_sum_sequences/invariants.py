"""Zero-sum constants of finite abelian groups.

For a set ``L`` of positive integers, ``s_L(G)`` is the smallest ``l`` such
that every sequence of length ``l`` over ``G`` has a nonempty zero-sum
subsequence whose length lies in ``L``, or infinity if there is no such
``l``. Classical invariants are special cases: the Davenport constant
``D(G) = s_N(G)``, the Erdős–Ginzburg–Ziv constant ``s(G) = s_{exp(G)}(G)``,
and ``eta(G) = s_[1, exp(G)](G)``.

Call a sequence *L-free* if it has no nonempty zero-sum subsequence with
length in ``L``. Subsequences of L-free sequences are L-free, so the
computation enumerates L-free sequences length by length, keeping one
representative per orbit of the chosen automorphism group.

For each sequence, a table records for every group element ``g`` the set of
lengths of subsequences with sum ``g``, as an integer bitset. Appending a
term ``h`` shifts the entry of ``g - h`` into that of ``g``; the lengths of
the new zero-sum subsequences are the shifted entry of ``-h``. Shared group
arithmetic is retained separately from the search.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from math import inf
from typing import Generic

from .additive_sequence import (
    AdditiveSequence,
    AdditiveSequenceSpace,
    Element,
    _immutable_term,
    _positive_integer,
)
from .automorphisms import AutomorphismGroup, AutomorphismGroupUnavailable


@dataclass(frozen=True)
class ZeroSumConstant(Generic[Element]):
    """A zero-sum constant together with its extremal representatives.

    ``value`` is ``s_L(G)``, or ``math.inf``. ``extremal`` holds the L-free
    sequences of maximal length ``value - 1``, one per orbit of the chosen
    symmetry group, each the minimum of its orbit. With no reduction it
    contains every extremal sequence; it is empty when ``value`` is infinite.
    """

    lengths: frozenset[int]
    value: int | float
    extremal: tuple[AdditiveSequence[Element], ...]


def exponent(space: AdditiveSequenceSpace) -> int:
    """The exponent of the finite base parent, without building lookup tables."""

    return space._group_arithmetic().exponent()


class _Symmetry:
    """Orbit minima and stabilizers for indexed multiplicity vectors.

    In term order, the lexicographically largest multiplicity vector in an
    orbit belongs to the smallest sequence, i.e. to the orbit minimum.
    """

    def __init__(self, table, group) -> None:
        count = len(table.elements)
        if group is None:
            self.forward = self.inverse = (tuple(range(count)),)
            return
        forward = []
        inverse = []
        for automorphism in group.elements:
            try:
                images = tuple(
                    table.index[_immutable_term(group.apply_term(automorphism, term))]
                    for term in table.elements
                )
            except KeyError:
                raise ValueError("symmetry must act on the space's ambient parent") from None
            if len(set(images)) != count or images[table.zero_index] != table.zero_index:
                raise ValueError("symmetry must act by zero-fixing bijections")
            preimages = [0] * count
            for position, image in enumerate(images):
                preimages[image] = position
            forward.append(images)
            inverse.append(tuple(preimages))
        self.forward = tuple(forward)
        self.inverse = tuple(inverse)

    def canonical(self, multiplicities: tuple[int, ...]):
        """The orbit minimum and the automorphisms mapping onto it."""

        candidates = range(len(self.inverse))
        key = []
        for position in range(len(multiplicities)):
            best, kept = -1, []
            for candidate in candidates:
                value = multiplicities[self.inverse[candidate][position]]
                if value > best:
                    best, kept = value, [candidate]
                elif value == best:
                    kept.append(candidate)
            key.append(best)
            candidates = kept
        return tuple(key), candidates

    def extension_terms(self, realising: list[int]) -> list[int]:
        """The least term of every orbit of a sequence's stabilizer.

        ``realising`` lists automorphisms ``a`` mapping the sequence onto
        its orbit minimum; the stabilizer consists of the ``a^-1 a0``.
        """

        first = self.forward[realising[0]]
        stabilizer = [
            [self.inverse[candidate][image] for image in first]
            for candidate in realising
        ]
        return [
            term
            for term in range(len(first))
            if all(table[term] >= term for table in stabilizer)
        ]


def _normalize_lengths(lengths: Iterable[int]) -> frozenset[int]:
    try:
        normalized = frozenset(
            _positive_integer(length, name="lengths") for length in lengths
        )
    except TypeError:
        raise TypeError("lengths must be an iterable of positive integers") from None
    if not normalized:
        raise ValueError("lengths must not be empty")
    return normalized


def _validate_symmetry(symmetry) -> None:
    if (
        symmetry is not True
        and symmetry is not False
        and not isinstance(symmetry, AutomorphismGroup)
    ):
        raise TypeError("symmetry must be True, False, or an AutomorphismGroup")


def _resolve_symmetry(space: AdditiveSequenceSpace, symmetry):
    if symmetry is False:
        return None
    if symmetry is not True:
        return symmetry
    try:
        return space.automorphism_group()
    except AutomorphismGroupUnavailable:
        return None


def zero_sum_constant(
    space: AdditiveSequenceSpace[Element],
    lengths: Iterable[int],
    *,
    symmetry: bool | AutomorphismGroup = True,
) -> ZeroSumConstant[Element]:
    """Compute ``s_L(G)`` for the finite base parent ``G`` of ``space``.

    ``lengths`` is a finite nonempty set ``L`` of positive integers. An
    infinite set is represented by an initial segment: if the result is at
    most ``max(L)``, it is ``s_L'(G)`` for every ``L'`` agreeing with ``L`` on
    ``[1, max(L)]``. The value is infinite exactly when ``L`` contains no
    multiple of the exponent of ``G``.

    ``symmetry=True`` discovers the space's automorphism group, falling
    back to no reduction when unavailable. ``False`` disables discovery
    and reduction. An explicit :class:`AutomorphismGroup` uses exactly
    that group. Symmetry affects the running time and representatives,
    never the value; discovery errors other than unavailability propagate.

    The configured Davenport bound of ``space`` is not used.
    """

    lengths = _normalize_lengths(lengths)
    _validate_symmetry(symmetry)
    table = space._group_arithmetic()
    e = table.exponent()
    if all(length % e for length in lengths):
        # An element of order exp(G) can be repeated indefinitely.
        return ZeroSumConstant(lengths, inf, ())
    action = _Symmetry(table, _resolve_symmetry(space, symmetry))

    count = len(table.elements)
    maximum_length = max(lengths)
    forbidden = 0
    negative = table.negative_indices
    difference = table.difference_table

    # orbit minimum -> (a sequence in the orbit, its subsequence-sum lengths)
    empty_sums = [0] * count
    empty_sums[table.zero_index] = 1
    level = {(0,) * count: ((0,) * count, tuple(empty_sums))}
    length = 0
    while True:
        # Grow bitsets with the search, not with potentially huge sparse L.
        if length + 1 in lengths:
            forbidden |= 1 << (length + 1)
        window = (1 << (min(maximum_length, length + 1) + 1)) - 1
        extensions = {}
        for multiplicities, sums in level.values():
            _, realising = action.canonical(multiplicities)
            for term in action.extension_terms(realising):
                if (sums[negative[term]] << 1) & forbidden:
                    continue
                extended = list(multiplicities)
                extended[term] += 1
                extended = tuple(extended)
                key, _ = action.canonical(extended)
                if key not in extensions:
                    extensions[key] = (
                        extended,
                        tuple(
                            sums[total]
                            | ((sums[difference[total][term]] << 1) & window)
                            for total in range(count)
                        ),
                    )
        if not extensions:
            break
        level = extensions
        length += 1

    extremal = tuple(
        space.from_multiplicities(
            {
                table.elements[term]: multiplicity
                for term, multiplicity in enumerate(key)
                if multiplicity
            }
        )
        for key in sorted(level, reverse=True)
    )
    return ZeroSumConstant(lengths, length + 1, extremal)


def davenport_constant(
    space: AdditiveSequenceSpace, *, symmetry: bool | AutomorphismGroup = True
) -> int:
    """Return and cache the exact Davenport constant of the base parent.

    Use a cached value or a supported exact formula when possible; otherwise
    search exhaustively. The caller's configured bound is not used.
    """

    _validate_symmetry(symmetry)
    arithmetic = space._group_arithmetic()
    known = arithmetic.known_davenport_constant()
    if known is not None:
        return known
    bound = arithmetic.davenport_bound()
    value = zero_sum_constant(space, range(1, bound + 1), symmetry=symmetry).value
    arithmetic.exact_davenport_constant = value
    return value


def erdos_ginzburg_ziv_constant(
    space: AdditiveSequenceSpace, *, symmetry: bool | AutomorphismGroup = True
) -> int:
    """The Erdős–Ginzburg–Ziv constant ``s(G) = s_{exp(G)}(G)``."""

    return zero_sum_constant(space, {exponent(space)}, symmetry=symmetry).value


def eta_constant(
    space: AdditiveSequenceSpace, *, symmetry: bool | AutomorphismGroup = True
) -> int:
    """The constant ``eta(G) = s_[1, exp(G)](G)``."""

    return zero_sum_constant(space, range(1, exponent(space) + 1), symmetry=symmetry).value
