"""Exact zero-sum constants and extremal sequences for finite abelian groups.

The classical constants and notation follow :cite:p:`gao-geroldinger-2006`,
Definition 2.1. :func:`zero_sum_constant` allows any finite set of positive
subsequence lengths. Its search includes the identity element, unlike the
package's factorization queries.

For related computational approaches, see :cite:p:`zerdoum-2021`, Chapters
IV–V, and :cite:p:`geroldinger-liebmann-philipp-2012`, Section 4.
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
    """An exact zero-sum constant and the longest sequences avoiding it.

    Attributes
    ----------
    lengths : frozenset[int]
        The allowed lengths of a nonempty zero-sum subsequence.
    value : int or float
        The least length forcing such a subsequence, or ``math.inf`` if no
        finite length suffices. This is the exact constant, not an upper bound.
    extremal : tuple[AdditiveSequence, ...]
        Sequences of length ``value - 1`` with no such subsequence, sorted in
        sequence order. With symmetry reduction, this contains the minimum
        of each orbit under the chosen group; otherwise it contains every
        extremal sequence. It is empty when ``value`` is infinite.
    """

    lengths: frozenset[int]
    value: int | float
    extremal: tuple[AdditiveSequence[Element], ...]


def exponent(space: AdditiveSequenceSpace) -> int:
    """Return the exact exponent of the finite abelian base group.

    This is the least common multiple of its element orders. Available group
    metadata is used first; otherwise the elements and their orders are
    enumerated. The result is cached. No Davenport bound is used.
    """

    return space._group_arithmetic().exponent()


class _Symmetry:
    """Orbit minima and stabilizers for indexed multiplicity vectors.

    In term order, the lexicographically largest multiplicity vector in an
    orbit belongs to the smallest sequence, i.e. to the orbit minimum.
    """

    def __init__(self, table, group) -> None:
        """Index the given group action, using the trivial action for None."""

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
        """Return a canonical vector and the indices of automorphisms realizing it."""

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
        """Choose one term per orbit of the stabilizer determined by ``realising``."""

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
    """Deduplicate allowed lengths, rejecting empty or nonpositive inputs."""

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
    """Reject options other than booleans and explicit automorphism groups."""

    if (
        symmetry is not True
        and symmetry is not False
        and not isinstance(symmetry, AutomorphismGroup)
    ):
        raise TypeError("symmetry must be True, False, or an AutomorphismGroup")


def _resolve_symmetry(space: AdditiveSequenceSpace, symmetry):
    """Return the requested group, or None when reduction is disabled or unavailable."""

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
    """Return the exact constant ``s_L(G)`` and its extremal sequences.

    This is the least length forcing a nonempty zero-sum subsequence whose
    length belongs to ``L``, or ``math.inf`` if no finite length suffices.
    ``lengths`` must be a finite, nonempty iterable of positive integers;
    duplicates are ignored. Sequences may use every element of the finite
    abelian base group, including zero.

    No Davenport bound is used or derived. The search extends avoiding
    sequences one term at a time, tracking their possible subsequence sums
    and lengths. Since every subsequence of an avoiding sequence also avoids
    the target, the first empty search level gives the exact constant.
    Finite searches require an iterable parent and can be expensive.

    ``symmetry=True`` discovers an automorphism group, using no reduction
    when one is unavailable. ``False`` disables discovery and reduction.
    An explicit :class:`AutomorphismGroup` uses precisely that group of
    automorphisms of the base group. Symmetry changes the representatives
    and search cost, not the value. Discovery errors other than
    :class:`AutomorphismGroupUnavailable` propagate.

    Notes
    -----
    The value is infinite exactly when ``L`` contains no multiple of the
    group exponent. To study an infinite allowed set, pass a finite
    truncation: a result at most ``max(L)`` also applies to every allowed
    set agreeing with ``L`` up to that cutoff.

    Examples
    --------
    >>> from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup
    >>> from zero_sum_sequences import zero_sum_constant
    >>> space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    >>> result = zero_sum_constant(space, {3}, symmetry=False)
    >>> result.value, len(result.extremal)
    (5, 3)
    >>> zero_sum_constant(space, {1, 2}).value
    inf
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
    """Return the exact Davenport constant ``D(G)`` of the finite abelian base group.

    This is the least length forcing a nonempty zero-sum subsequence, with
    no restriction on its length. A cached value or a supported exact formula
    is used when available; otherwise :func:`zero_sum_constant` performs an
    exhaustive search. The result is cached; spaces over the same
    :class:`FiniteAdditiveGroup` share the cached value.

    The search derives its own certified upper bound to cover all relevant
    subsequence lengths. The space's configured Davenport bound is ignored,
    even if supplied explicitly. The returned value is always exact, not
    that upper bound. An exhaustive computation can be expensive;
    ``symmetry`` has the meaning described in :func:`zero_sum_constant`.

    Notes
    -----
    For known invariant factors ``1 < n_1 | ... | n_r``, the formula
    ``D(G) = 1 + sum(n_i - 1)`` is used for p-groups and groups of rank at
    most two. This is the result of Kruyswijk and Olson recorded in
    :cite:p:`gao-geroldinger-2006`, Theorem 3.1; see also
    :cite:p:`olson-1969-i,olson-1969-ii`. A group known to be cyclic has
    ``D(G) = |G|``, and the trivial group has ``D(G) = 1``.

    An exact result can tighten an automatically derived
    :attr:`AdditiveSequenceSpace.davenport_bound`, but does not replace an
    explicit bound.

    Examples
    --------
    >>> from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup
    >>> from zero_sum_sequences import davenport_constant
    >>> space = AdditiveSequenceSpace(
    ...     FiniteAdditiveGroup.cyclic_product(2, 4), davenport_bound=8)
    >>> davenport_constant(space)
    5
    >>> space.davenport_bound
    8
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
    """Return the exact Erdős–Ginzburg–Ziv constant ``s(G)``.

    This is the least length forcing a zero-sum subsequence of length
    ``exp(G)``; see :cite:p:`gao-geroldinger-2006`, Definition 2.1.
    The group exponent is computed, then :func:`zero_sum_constant` searches
    with that single allowed length. No Davenport bound is used or derived.
    The search can be expensive; ``symmetry`` is passed through unchanged.
    """

    return zero_sum_constant(space, {exponent(space)}, symmetry=symmetry).value


def eta_constant(
    space: AdditiveSequenceSpace, *, symmetry: bool | AutomorphismGroup = True
) -> int:
    """Return the exact short zero-sum constant ``eta(G)``.

    This is the least length forcing a nonempty zero-sum subsequence of
    length at most ``exp(G)``; see :cite:p:`gao-geroldinger-2006`, Definition
    2.1. The group exponent is computed, then :func:`zero_sum_constant`
    searches with allowed lengths ``1, ..., exp(G)``. No Davenport bound is
    used or derived. The search can be expensive; ``symmetry`` is passed
    through unchanged.
    """

    return zero_sum_constant(space, range(1, exponent(space) + 1), symmetry=symmetry).value
