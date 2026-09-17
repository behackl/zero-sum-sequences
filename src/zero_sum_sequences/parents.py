"""Small adapters for defining additive parents without a CAS dependency."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from itertools import product
from math import gcd
from operator import index
from typing import Generic, TypeVar

Element = TypeVar("Element")


def default_encode_term(term: object) -> object:
    """Encode an integer as itself and an iterable of integers as a list."""

    if isinstance(term, bool):
        raise TypeError("cannot encode a boolean term")
    if isinstance(term, int):
        return term
    try:
        return [int(coordinate) for coordinate in term]
    except TypeError:
        raise TypeError(
            f"cannot encode {term!r}; supply encode_term/decode_term on the parent"
        ) from None


def _positive_modulus(value: object) -> int:
    if isinstance(value, bool):
        raise ValueError("cyclic moduli must be positive integers")
    try:
        modulus = index(value)
    except TypeError:
        raise ValueError("cyclic moduli must be positive integers") from None
    if modulus < 1:
        raise ValueError("cyclic moduli must be positive integers")
    return modulus


def _scale_coordinate(coordinate: int, unit: int, modulus: int):
    def scaling(element):
        image = list(element)
        image[coordinate] = unit * image[coordinate] % modulus
        return tuple(image)

    return scaling


def _shear_coordinate(
    target: int,
    source: int,
    coefficient: int,
    modulus: int,
):
    def shear(element):
        image = list(element)
        image[target] = (
            image[target] + coefficient * image[source]
        ) % modulus
        return tuple(image)

    return shear


class FiniteAdditiveGroup(Generic[Element]):
    """A finite additive parent defined by elements and an operation.

    The supplied data is treated as a group presentation; associativity,
    inverses, and closure are not checked exhaustively.  Results produced by
    ``add`` are coerced and checked for membership when the operation is used.
    """

    __slots__ = (
        "_elements",
        "_element_set",
        "_zero",
        "_operation",
        "_coerce",
        "_additive_generators",
        "_automorphism_generators",
        "_automorphism_group_provider",
        "_automorphism_group",
        "_encode_term",
        "_decode_term",
    )

    @classmethod
    def cyclic_product(cls, *moduli: int) -> FiniteAdditiveGroup[tuple[int, ...]]:
        """Return the coordinate group ``C_n1 x ... x C_nr``.

        The parent includes elementary generators of its full automorphism
        group: multiplication of one coordinate by a unit and all admissible
        coordinate shears.  Consequently, sequence orbit methods work without
        an explicitly supplied action.
        """

        if not moduli:
            raise ValueError("at least one cyclic modulus is required")
        normalized_moduli = tuple(_positive_modulus(value) for value in moduli)

        def coerce(value):
            try:
                coordinates = tuple(value)
            except TypeError:
                raise ValueError(
                    "cyclic-product elements must be coordinate iterables"
                ) from None
            if len(coordinates) != len(normalized_moduli):
                raise ValueError(
                    "cyclic-product elements have the wrong number of coordinates"
                )
            try:
                return tuple(
                    index(coordinate) % modulus
                    for coordinate, modulus in zip(
                        coordinates,
                        normalized_moduli,
                    )
                )
            except TypeError:
                raise ValueError(
                    "cyclic-product coordinates must be integers"
                ) from None

        def add(left, right):
            return tuple(
                (left[position] + right[position]) % modulus
                for position, modulus in enumerate(normalized_moduli)
            )

        generators = []
        for coordinate, modulus in enumerate(normalized_moduli):
            generators.extend(
                _scale_coordinate(coordinate, unit, modulus)
                for unit in range(2, modulus)
                if gcd(unit, modulus) == 1
            )
        for target, target_modulus in enumerate(normalized_moduli):
            for source, source_modulus in enumerate(normalized_moduli):
                if target == source:
                    continue
                coefficient = target_modulus // gcd(
                    target_modulus,
                    source_modulus,
                )
                if coefficient % target_modulus:
                    generators.append(
                        _shear_coordinate(
                            target,
                            source,
                            coefficient,
                            target_modulus,
                        )
                    )

        additive_generators = tuple(
            tuple(
                1 if position == coordinate else 0
                for position in range(len(normalized_moduli))
            )
            for coordinate, modulus in enumerate(normalized_moduli)
            if modulus > 1
        )

        return cls(
            product(*(range(modulus) for modulus in normalized_moduli)),
            zero=(0,) * len(normalized_moduli),
            add=add,
            coerce=coerce,
            additive_generators=additive_generators,
            automorphism_generators=generators,
        )

    def __init__(
        self,
        elements: Iterable[Element],
        *,
        zero: Element,
        add: Callable[[Element, Element], Element],
        coerce: Callable[[object], Element] | None = None,
        additive_generators: Iterable[Element] | None = None,
        automorphism_generators: (
            Iterable[Callable[[Element], Element]] | None
        ) = None,
        automorphism_group: Callable[[], object] | None = None,
        encode_term: Callable[[Element], object] | None = None,
        decode_term: Callable[[object], Element] | None = None,
    ) -> None:
        if not callable(add):
            raise TypeError("the additive operation must be callable")
        if coerce is not None and not callable(coerce):
            raise TypeError("coerce must be callable")

        self._coerce = coerce if coerce is not None else lambda value: value
        try:
            coerced_elements = tuple(self._coerce(element) for element in elements)
            element_set = frozenset(coerced_elements)
        except TypeError:
            raise TypeError("finite-group elements must be hashable") from None
        if len(element_set) != len(coerced_elements):
            raise ValueError("finite-group elements must be distinct after coercion")

        try:
            ordered_elements = tuple(sorted(coerced_elements))
        except TypeError:
            raise TypeError(
                "finite-group elements must be mutually orderable"
            ) from None

        self._elements = ordered_elements
        self._element_set = element_set
        self._operation = add
        self._zero = self(zero)
        if additive_generators is None:
            self._additive_generators = None
        else:
            try:
                generators = tuple(self(element) for element in additive_generators)
            except TypeError:
                raise TypeError("additive generators must be iterable") from None
            if self._zero in generators:
                raise ValueError("additive generators must be nonzero")
            if len(set(generators)) != len(generators):
                raise ValueError("additive generators must be distinct")
            self._additive_generators = generators
        if automorphism_generators is None:
            self._automorphism_generators = None
        else:
            generators = tuple(automorphism_generators)
            if not all(callable(generator) for generator in generators):
                raise TypeError("automorphism generators must be callable")
            self._automorphism_generators = generators
        if automorphism_group is not None and not callable(automorphism_group):
            raise TypeError("automorphism_group must be callable")
        self._automorphism_group_provider = automorphism_group
        self._automorphism_group = None
        if (encode_term is None) != (decode_term is None):
            raise TypeError("encode_term and decode_term must be supplied together")
        if encode_term is not None and not (callable(encode_term) and callable(decode_term)):
            raise TypeError("encode_term and decode_term must be callable")
        self._encode_term = encode_term
        self._decode_term = decode_term

    def __call__(self, value: object) -> Element:
        element = self._coerce(value)
        try:
            contained = element in self._element_set
        except TypeError:
            contained = False
        if not contained:
            raise ValueError(f"{element!r} is not an element of this group")
        return element

    def zero(self) -> Element:
        """Return the additive identity."""

        return self._zero

    def add(self, left: Element, right: Element) -> Element:
        """Add two elements and return their canonical representative."""

        left = self(left)
        right = self(right)
        return self(self._operation(left, right))

    def is_finite(self) -> bool:
        return True

    def additive_generators(self):
        """Return a configured generating set for the additive group."""

        if self._additive_generators is None:
            raise NotImplementedError
        return self._additive_generators

    def automorphism_generators(self):
        """Return explicitly configured automorphism generators.

        Generic finite groups do not expose a presentation from which these
        maps could be inferred.  Supply them at construction when automatic
        orbit methods are desired.
        """

        if self._automorphism_generators is None:
            raise NotImplementedError
        return self._automorphism_generators

    def encode_term(self, term: Element) -> object:
        """Return a JSON-compatible encoding of ``term``.

        A codec supplied at construction is used; otherwise integers are
        returned as such and iterable terms as lists of integers.
        """

        if self._encode_term is not None:
            return self._encode_term(term)
        return default_encode_term(term)

    def decode_term(self, data: object) -> Element:
        """Return the element encoded by ``data`` (see :meth:`encode_term`)."""

        if self._decode_term is not None:
            return self(self._decode_term(data))
        return self(data)

    def automorphism_group(self):
        """Return the materialized automorphism group of this parent.

        A group supplied at construction takes precedence; otherwise the
        configured automorphism generators are closed under composition.
        The result is cached.
        """

        if self._automorphism_group is None:
            from .automorphisms import AutomorphismGroup

            if self._automorphism_group_provider is not None:
                group = self._automorphism_group_provider()
                if not isinstance(group, AutomorphismGroup):
                    raise TypeError(
                        "the automorphism group provider must return an "
                        "AutomorphismGroup"
                    )
            elif self._automorphism_generators is not None:
                group = AutomorphismGroup.closure(self, self._automorphism_generators)
            else:
                raise NotImplementedError
            self._automorphism_group = group
        return self._automorphism_group

    def __iter__(self) -> Iterator[Element]:
        return iter(self._elements)

    def __len__(self) -> int:
        return len(self._elements)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({list(self._elements)!r}, zero={self._zero!r})"
