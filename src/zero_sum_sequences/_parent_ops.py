"""Parent protocol primitives, separate from finite-group metadata."""

from copy import copy
from operator import add as default_add
from typing import TypeVar

from ._exceptions import _ParentIterationUnavailable
from .parents import FiniteAdditiveGroup

_MISSING = object()
Term = TypeVar("Term")


def _immutable_term(term: Term) -> Term:
    """Return a hashable term, copying supported mutable elements if needed."""

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


def _iter_parent(parent):
    """Iterate a parent, distinguishing absent protocols from provider bugs."""

    parent_type = type(parent)
    iterator = getattr(parent_type, "__iter__", _MISSING)
    if iterator is None or (
        iterator is _MISSING and getattr(parent_type, "__getitem__", None) is None
    ):
        raise _ParentIterationUnavailable("the parent is not iterable")
    return iter(parent)


def _is_finite(parent) -> bool:
    """Return whether the parent declares itself finite.

    A missing ``is_finite`` or one raising ``NotImplementedError`` means
    unknown, treated as not finite. Other provider errors propagate.
    """

    finite = getattr(parent, "is_finite", None)
    if not callable(finite):
        return False
    try:
        return bool(finite())
    except NotImplementedError:
        return False


def _native_term_arithmetic(parent) -> bool:
    """Whether the terms' own ``+``, ``-`` and orders are the group operation.

    A parent-supplied ``add(left, right)`` overrides the terms' arithmetic,
    e.g. integer labels combined by XOR.
    """

    return not callable(getattr(parent, "add", None))


def _addition(parent):
    """Return addition with canonical outputs, without finite-group queries."""

    if isinstance(parent, FiniteAdditiveGroup):
        # The package adapter already coerces inputs and results.
        return parent.add
    if not _native_term_arithmetic(parent):
        operation = parent.add
    elif type(parent).__module__.startswith("sage."):
        # Native Sage addition is closed in its parent. Trust it without
        # importing Sage or querying group metadata.
        return default_add
    else:
        operation = default_add

    def add(left, right):
        return parent(operation(left, right))

    return add
