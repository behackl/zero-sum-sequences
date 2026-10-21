"""Cached arithmetic and Davenport estimates for finite abelian parents.

Group metadata and element queries are kept separate from the addition and
subtraction tables needed by exhaustive searches.
"""

from __future__ import annotations

from functools import cached_property
from math import gcd, lcm

from ._exceptions import _ArithmeticUnavailable, _ParentIterationUnavailable
from ._parent_ops import (
    _addition,
    _canonical_addition,
    _immutable_term,
    _is_finite,
    _iter_parent,
    _native_term_arithmetic,
)
from ._validation import _positive_integer

_AUTOMATIC_TRIAL_LIMIT = 1024


def _invariant_factors(moduli):
    """Return nontrivial invariant factors from arbitrary cyclic orders."""

    factors = sorted(_positive_integer(n, name="cyclic order") for n in moduli)
    for i in range(len(factors)):
        for j in range(i + 1, len(factors)):
            factors[i], factors[j] = gcd(factors[i], factors[j]), lcm(
                factors[i], factors[j]
            )
    return tuple(n for n in factors if n > 1)


def _is_prime_power(n: int, *, max_steps: int | None = None) -> bool | None:
    """Return whether n is a prime power, or None if trial division gives up.

    ``max_steps`` bounds the number of trial divisors; once the smallest
    prime factor is found, the remaining divisions are at most ``log2(n)``.
    """

    if n < 2:
        return False
    p, steps = 2, 0
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            return n == 1
        steps += 1
        if max_steps is not None and steps >= max_steps:
            return None
        p += 1 if p == 2 else 2
    return True  # n is prime


class _GroupArithmetic:
    """Lazy arithmetic data for one finite abelian parent.

    Optional cyclic ``moduli`` give element orders and inverses directly.
    With ``native=True``, parent metadata and compatible element arithmetic
    are also used. Element queries are cached; indexed tables are built
    only when requested.
    """

    def __init__(self, parent, *, native: bool = True, moduli=None) -> None:
        """Attach a parent and its optional cyclic decomposition."""

        self.parent = parent
        self._add = _addition(parent)
        self._native = native
        self._native_terms = native and _native_term_arithmetic(parent)
        self._moduli = moduli
        self.zero = _immutable_term(parent(parent.zero()))
        self._orders = {self.zero: 1}
        self._negatives = {self.zero: self.zero}
        self._order = None
        self._exponent = None
        self._bound = None
        self.exact_davenport_constant = None

    def _require_finite(self) -> None:
        """Require the parent to declare itself finite."""

        if not _is_finite(self.parent):
            raise _ArithmeticUnavailable("group arithmetic requires a finite parent")

    def coerce(self, term):
        """Return an immutable, canonical parent element."""

        return _immutable_term(self.parent(term))

    def add(self, left, right):
        """Add two terms using the parent operation when supplied."""

        return _immutable_term(self._add(left, right))

    @cached_property
    def add_canonical(self):
        """Addition of canonical, immutable elements, coercing only the result."""

        return _canonical_addition(self.parent)

    @cached_property
    def elements(self):
        """The group's distinct elements, immutable and sorted."""

        self._require_finite()
        try:
            terms = _iter_parent(self.parent)
        except _ParentIterationUnavailable as error:
            raise _ParentIterationUnavailable(
                "indexed group arithmetic requires an iterable parent"
            ) from error
        elements = tuple(sorted(self.coerce(term) for term in terms))
        if len(set(elements)) != len(elements):
            raise ValueError("the parent contains repeated elements after coercion")
        if self.zero not in elements:
            raise ValueError("the parent does not contain its zero element")
        return elements

    @cached_property
    def index(self):
        """Map each element to its position in ``elements``."""

        return {term: i for i, term in enumerate(self.elements)}

    @cached_property
    def zero_index(self):
        """The position of the additive identity."""

        return self.index[self.zero]

    def order(self) -> int:
        """Return the group size, using metadata or enumerating its elements."""

        self._require_finite()
        if self._order is None:
            for name in ("cardinality", "order"):
                provider = getattr(self.parent, name, None)
                if callable(provider):
                    try:
                        value = provider()
                    except NotImplementedError:
                        continue
                    self._order = _positive_integer(value, name="group order")
                    break
            if self._order is None and getattr(type(self.parent), "__len__", None) is not None:
                try:
                    value = len(self.parent)
                except NotImplementedError:
                    pass
                else:
                    self._order = _positive_integer(value, name="group order")
            if self._order is None:
                try:
                    self._order = _positive_integer(len(self.elements), name="group order")
                except (_ParentIterationUnavailable, NotImplementedError) as error:
                    raise _ArithmeticUnavailable(
                        "group order requires metadata or an iterable parent"
                    ) from error
        return self._order

    def element_order(self, term) -> int:
        """Return the term's additive order, using repeated addition if needed."""

        self._require_finite()
        term = self.coerce(term)
        if term not in self._orders:
            if self._moduli is not None:
                value = lcm(*(n // gcd(g, n) for g, n in zip(term, self._moduli)))
            else:
                provider = (
                    getattr(self.parent, "element_order", None) if self._native else None
                )
                if not callable(provider):
                    provider = (
                        getattr(term, "additive_order", None) if self._native_terms else None
                    )
                    arguments = ()
                else:
                    arguments = (term,)
                value = None
                if callable(provider):
                    try:
                        value = provider(*arguments)
                    except NotImplementedError:
                        pass
                if value is None:
                    multiple = self.zero
                    for value in range(1, self.order() + 1):
                        previous = multiple
                        multiple = self.add(multiple, term)
                        if multiple == self.zero:
                            self._negatives[term] = previous
                            break
                    else:
                        raise ValueError("the parent does not provide finite additive orders")
            self._orders[term] = _positive_integer(value, name="element order")
        return self._orders[term]

    def negate(self, term):
        """Return and cache the term's additive inverse."""

        term = self.coerce(term)
        if term in self._negatives:
            return self._negatives[term]
        if self._moduli is not None:
            negative = self.coerce(tuple(-g for g in term))
            self._negatives[term] = negative
            return negative
        provider = getattr(self.parent, "negate", None) if self._native else None
        if callable(provider):
            try:
                negative = self.coerce(provider(term))
            except NotImplementedError:
                pass
            else:
                self._negatives[term] = negative
                return negative
        if self._native_terms:
            try:
                negative = self.coerce(-term)
            except (TypeError, ValueError):
                pass
            else:
                self._negatives[term] = negative
                return negative
        order = self.element_order(term)
        if term not in self._negatives:
            negative = self.zero
            for _ in range(order - 1):
                negative = self.add(negative, term)
            self._negatives[term] = negative
        return self._negatives[term]

    def _known_exponent(self):
        """Return a known exponent, or None without enumerating element orders."""

        if self._exponent is None:
            if self._moduli is not None:
                self._exponent = lcm(*self._moduli)
            elif self._native:
                provider = getattr(self.parent, "exponent", None)
                if callable(provider):
                    try:
                        value = provider()
                    except NotImplementedError:
                        pass
                    else:
                        self._exponent = _positive_integer(value, name="group exponent")
        return self._exponent

    def exponent(self) -> int:
        """Return the exact exponent, computing element orders if necessary."""

        self._require_finite()
        if self._exponent is None:
            if self._known_exponent() is None:
                self._exponent = lcm(*(self.element_order(g) for g in self.elements))
            self._bound = None  # Newly known metadata can also tighten a bound.
        return self._exponent

    def _record_exact_davenport(self, value: int) -> None:
        """Retain a certified exact value, rejecting contradictory results."""

        value = _positive_integer(value, name="Davenport constant")
        if self.exact_davenport_constant not in (None, value):
            raise ValueError("inconsistent exact Davenport constants")
        self.exact_davenport_constant = value

    def known_davenport_constant(self, *, allow_expensive: bool = False):
        """Return a certified cached or formula value, or None without searching.

        Supported formulas are those described in :func:`davenport_constant`.
        Prime-power recognition has a fixed trial budget unless an explicit
        exact query requests full recognition with ``allow_expensive=True``.
        """

        self._require_finite()
        if self.exact_davenport_constant is not None:
            return self.exact_davenport_constant
        n = self.order()
        factors = None
        if self._moduli is not None:
            factors = _invariant_factors(self._moduli)
        elif self._native:
            provider = getattr(self.parent, "invariants", None)
            if callable(provider):
                try:
                    factors = _invariant_factors(provider())
                except NotImplementedError:
                    pass
                if factors is not None:
                    product = 1
                    for factor in factors:
                        product *= factor
                    if product != n:
                        raise ValueError("cyclic invariants do not match the group order")
        if n == 1:
            value = 1
        elif factors is not None and (
            len(factors) <= 2 or _is_prime_power(
                factors[-1], max_steps=None if allow_expensive else _AUTOMATIC_TRIAL_LIMIT
            ) is True
        ):
            value = 1 + sum(factor - 1 for factor in factors)
        elif self._known_exponent() == n:
            value = n  # a finite abelian group of exponent n and order n is cyclic
        else:
            return None
        self._record_exact_davenport(value)
        return value

    def davenport_bound(self) -> int:
        """Return a cached certified upper bound, without an exhaustive search.

        See :attr:`AdditiveSequenceSpace.davenport_bound` for the estimate.
        """

        if self.exact_davenport_constant is not None:
            return self.exact_davenport_constant
        if self._bound is None:
            exact = self.known_davenport_constant()
            if exact is not None:
                return exact
            n = self.order()
            e = self._known_exponent()
            if e is None:
                self._bound = n
            elif n % e:
                raise ValueError("the group exponent must divide its order")
            else:
                # van Emde Boas–Kruyswijk (1969), Theorem 7.1, p. 21:
                # D(G) <= e * (1 + ln(n/e)). Their mu(G) is our D(G).
                # ceil(log2(n/e)) >= ln(n/e) gives a conservative integer bound.
                self._bound = min(n, e * (1 + (n // e - 1).bit_length()))
        return self._bound

    @cached_property
    def addition_table(self):
        """Indices of pairwise sums, in ``elements`` order."""

        try:
            return tuple(
                tuple(self.index[self.add(g, h)] for h in self.elements)
                for g in self.elements
            )
        except KeyError:
            raise ValueError("the parent is not closed under addition") from None

    @cached_property
    def negative_indices(self):
        """The index of each element's additive inverse."""

        try:
            indices = tuple(row.index(self.zero_index) for row in self.addition_table)
        except ValueError:
            raise ValueError("indexed group arithmetic requires additive inverses") from None
        self._negatives.update(
            (g, self.elements[i]) for g, i in zip(self.elements, indices)
        )
        return indices

    @cached_property
    def difference_table(self):
        """Indices of pairwise differences, in ``elements`` order."""

        return tuple(
            tuple(row[i] for i in self.negative_indices)
            for row in self.addition_table
        )
