"""Equalities between factorizations of additive sequences."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Generic

from .additive_sequence import AdditiveSequence, AdditiveSequenceSpace, Element
from .factorization import Factorization


def _factor_key(
    factor: AdditiveSequence[Element],
) -> tuple[int, tuple[Element, ...]]:
    return len(factor), tuple(factor)


def _normalize_factorization(
    factors: Iterable[AdditiveSequence[Element]],
    *,
    name: str,
) -> Factorization[Element]:
    try:
        iterator = iter(factors)
    except TypeError:
        raise TypeError(f"{name} factorization must be iterable") from None

    normalized = []
    for factor in iterator:
        if not isinstance(factor, AdditiveSequence):
            raise TypeError(f"{name} factors must be additive sequences")
        normalized.append(factor)
    return tuple(sorted(normalized, key=_factor_key))


def _resolve_space(
    source: Factorization[Element],
    target: Factorization[Element],
    space: AdditiveSequenceSpace[Element] | None,
) -> AdditiveSequenceSpace[Element]:
    if space is not None and not isinstance(space, AdditiveSequenceSpace):
        raise TypeError("space must be an AdditiveSequenceSpace")

    factors = (*source, *target)
    if space is None:
        if not factors:
            raise ValueError("cannot infer the sequence space of an empty relation")
        space = factors[0].parent()
    if any(factor.parent() is not space for factor in factors):
        raise TypeError("relation factors use different sequence spaces")
    return space


def _product_multiplicities(
    factors: Factorization[Element],
) -> Counter[Element]:
    multiplicities: Counter[Element] = Counter()
    for factor in factors:
        multiplicities.update(factor.multiplicities)
    return multiplicities


def _counter_factorization(
    counts: Counter[AdditiveSequence[Element]],
) -> Factorization[Element]:
    return tuple(sorted(counts.elements(), key=_factor_key))


@dataclass(frozen=True, slots=True, init=False)
class FactorizationRelation(Generic[Element]):
    """An oriented equality between two unordered factorizations.

    Factors are canonicalized as multisets, but the source and target sides
    are not interchanged. The constructor verifies compatible sequence spaces
    and equal sequence products. It assumes, without recomputing, that the
    supplied factors are atoms when atom factorizations are required.
    """

    source: Factorization[Element]
    target: Factorization[Element]
    product: AdditiveSequence[Element] = field(repr=False)
    _common_factor: Factorization[Element] = field(
        repr=False,
        compare=False,
        hash=False,
    )
    _reduced_source: Factorization[Element] = field(
        repr=False,
        compare=False,
        hash=False,
    )
    _reduced_target: Factorization[Element] = field(
        repr=False,
        compare=False,
        hash=False,
    )

    def __init__(
        self,
        source: Iterable[AdditiveSequence[Element]],
        target: Iterable[AdditiveSequence[Element]],
        *,
        space: AdditiveSequenceSpace[Element] | None = None,
    ) -> None:
        normalized_source = _normalize_factorization(source, name="source")
        normalized_target = _normalize_factorization(target, name="target")
        relation_space = _resolve_space(
            normalized_source,
            normalized_target,
            space,
        )

        source_product = _product_multiplicities(normalized_source)
        target_product = _product_multiplicities(normalized_target)
        if source_product != target_product:
            raise ValueError(
                "source and target factorizations have different products"
            )

        source_counts = Counter(normalized_source)
        target_counts = Counter(normalized_target)
        common_counts = source_counts & target_counts
        common_factor = _counter_factorization(common_counts)
        reduced_source = _counter_factorization(source_counts - common_counts)
        reduced_target = _counter_factorization(target_counts - common_counts)

        object.__setattr__(self, "source", normalized_source)
        object.__setattr__(self, "target", normalized_target)
        object.__setattr__(
            self,
            "product",
            relation_space.from_multiplicities(source_product),
        )
        object.__setattr__(self, "_common_factor", common_factor)
        object.__setattr__(self, "_reduced_source", reduced_source)
        object.__setattr__(self, "_reduced_target", reduced_target)

    @property
    def common_factor(self) -> Factorization[Element]:
        """Return the multiset gcd of the source and target."""

        return self._common_factor

    @property
    def reduced_source(self) -> Factorization[Element]:
        """Return the source after cancelling the common factor."""

        return self._reduced_source

    @property
    def reduced_target(self) -> Factorization[Element]:
        """Return the target after cancelling the common factor."""

        return self._reduced_target

    @property
    def is_reduced(self) -> bool:
        """Return whether the two sides have no factor in common."""

        return not self._common_factor

    @property
    def distance(self) -> int:
        """Return the standard distance between the two factorizations."""

        return max(len(self._reduced_source), len(self._reduced_target))

    def reduced(self) -> FactorizationRelation[Element]:
        """Return the relation obtained by cancelling its common factor."""

        if self.is_reduced:
            return self
        return type(self)(
            self._reduced_source,
            self._reduced_target,
            space=self.product.parent(),
        )
