"""Small adapters for defining additive parents without a CAS dependency."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from typing import Generic, TypeVar

Element = TypeVar("Element")


class FiniteAdditiveGroup(Generic[Element]):
    """A finite additive parent defined by elements and an operation.

    The supplied data is treated as a group presentation; associativity,
    inverses, and closure are not checked exhaustively.  Results produced by
    ``add`` are coerced and checked for membership when the operation is used.
    """

    __slots__ = ("_elements", "_element_set", "_zero", "_operation", "_coerce")

    def __init__(
        self,
        elements: Iterable[Element],
        *,
        zero: Element,
        add: Callable[[Element, Element], Element],
        coerce: Callable[[object], Element] | None = None,
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

    def __iter__(self) -> Iterator[Element]:
        return iter(self._elements)

    def __len__(self) -> int:
        return len(self._elements)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({list(self._elements)!r}, zero={self._zero!r})"
