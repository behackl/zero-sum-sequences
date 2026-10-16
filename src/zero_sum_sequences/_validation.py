"""Exact integer validation shared by parent and sequence APIs."""

from operator import index


def _non_negative_integer(value: object, *, name: str) -> int:
    """Return an exact nonnegative indexable integer, rejecting booleans."""

    if isinstance(value, bool):
        raise ValueError(f"{name} must be a non-negative integer")
    try:
        integer = index(value)
    except TypeError:
        raise ValueError(f"{name} must be a non-negative integer") from None
    if integer < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return integer


def _positive_integer(value: object, *, name: str) -> int:
    """Return an exact positive indexable integer, rejecting booleans."""

    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive integer")
    try:
        integer = index(value)
    except TypeError:
        raise ValueError(f"{name} must be a positive integer") from None
    if integer < 1:
        raise ValueError(f"{name} must be a positive integer")
    return integer
