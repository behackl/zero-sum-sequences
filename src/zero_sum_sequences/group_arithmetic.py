"""Cached arithmetic and Davenport estimates for finite abelian parents.

Group metadata and element queries are kept separate from the addition and
subtraction tables needed by exhaustive searches.
"""

from __future__ import annotations

from functools import cached_property
from math import gcd, lcm
from operator import add as default_add

from .additive_sequence import _immutable_term, _positive_integer


def _invariant_factors(moduli):
    """Return nontrivial invariant factors from arbitrary cyclic orders."""

    factors = sorted(_positive_integer(n, name="cyclic order") for n in moduli)
    for i in range(len(factors)):
        for j in range(i + 1, len(factors)):
            factors[i], factors[j] = gcd(factors[i], factors[j]), lcm(
                factors[i], factors[j]
            )
    return tuple(n for n in factors if n > 1)


def _is_prime_power(n: int) -> bool:
    """Return whether n is a positive power of a prime."""

    if n < 2:
        return False
    p = 2
    while p * p <= n and n % p:
        p += 1
    if n % p:
        return True  # n is prime
    while n % p == 0:
        n //= p
    return n == 1


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
        self._native = native
        self._moduli = moduli
        self.zero = _immutable_term(parent(parent.zero()))
        self._orders = {self.zero: 1}
        self._negatives = {self.zero: self.zero}
        self._order = None
        self._exponent = None
        self.exact_davenport_constant = None

    def _require_finite(self) -> None:
        """Require the parent to declare itself finite."""

        finite = getattr(self.parent, "is_finite", None)
        if not callable(finite) or not finite():
            raise ValueError("group arithmetic requires a finite parent")

    def coerce(self, term):
        """Return an immutable, canonical parent element."""

        return _immutable_term(self.parent(term))

    def add(self, left, right):
        """Add two terms using the parent operation when supplied."""

        operation = getattr(self.parent, "add", default_add)
        if not callable(operation):
            operation = default_add
        return self.coerce(operation(left, right))

    @cached_property
    def elements(self):
        """The group's distinct elements, immutable and sorted."""

        self._require_finite()
        try:
            terms = tuple(self.parent)
        except TypeError:
            raise TypeError("indexed group arithmetic requires an iterable parent") from None
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
            if self._order is None:
                try:
                    value = len(self.parent)
                except TypeError:
                    value = len(self.elements)
                self._order = _positive_integer(value, name="group order")
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
                    # Native term arithmetic is only authoritative when the
                    # parent has not supplied a different addition operation.
                    provider = (
                        getattr(term, "additive_order", None)
                        if self._native and not callable(getattr(self.parent, "add", None))
                        else None
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
        if self._native and not callable(getattr(self.parent, "add", None)):
            try:
                negative = self.coerce(-term)
            except (TypeError, ValueError):
                pass
            else:
                self._negatives[term] = negative
                return negative
        # A custom parent addition may differ from the terms' native
        # arithmetic (e.g. integer labels with XOR as the group operation).
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
        if self._known_exponent() is None:
            self._exponent = lcm(*(self.element_order(g) for g in self.elements))
        return self._exponent

    def known_davenport_constant(self):
        """Return an exact cached or formula value, or None without searching.

        Supported formulas are those described in :func:`davenport_constant`:
        cyclic groups, rank at most two, and p-groups with known invariant
        factors. The trivial group has value one.
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
            len(factors) <= 2 or _is_prime_power(factors[-1])
        ):
            value = 1 + sum(factor - 1 for factor in factors)
        elif self._known_exponent() == n:
            value = n  # a finite abelian group of exponent n and order n is cyclic
        else:
            return None
        self.exact_davenport_constant = value
        return value

    def davenport_bound(self) -> int:
        """Return a certified upper bound without an exhaustive constant search.

        Use an exact value when known; otherwise use the order/exponent
        estimate in :attr:`AdditiveSequenceSpace.davenport_bound`, falling
        back to the group order when no exponent is available.
        """

        exact = self.known_davenport_constant()
        if exact is not None:
            return exact
        n = self.order()
        e = self._known_exponent()
        if e is None:
            return n
        if n % e:
            raise ValueError("the group exponent must divide its order")
        # van Emde Boas–Kruyswijk (1969), Theorem 7.1, p. 21:
        # D(G) <= e * (1 + ln(n/e)). Their mu(G) is our D(G).
        # ceil(log2(n/e)) >= ln(n/e) gives a conservative integer bound.
        return min(n, e * (1 + (n // e - 1).bit_length()))

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
