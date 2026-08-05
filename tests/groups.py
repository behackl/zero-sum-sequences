"""Small dependency-free additive parents used by the core test suite."""

from itertools import product

from zero_sum_sequences import FiniteAdditiveGroup


def cyclic_group(modulus: int) -> FiniteAdditiveGroup[int]:
    return FiniteAdditiveGroup(
        range(modulus),
        zero=0,
        add=lambda left, right: (left + right) % modulus,
        coerce=lambda value: int(value) % modulus,
    )


def product_of_cyclic_groups(*moduli: int) -> FiniteAdditiveGroup[tuple[int, ...]]:
    return FiniteAdditiveGroup(
        product(*(range(modulus) for modulus in moduli)),
        zero=(0,) * len(moduli),
        add=lambda left, right: tuple(
            (left[index] + right[index]) % modulus
            for index, modulus in enumerate(moduli)
        ),
        coerce=tuple,
    )


class IndexableInteger:
    def __init__(self, value: int) -> None:
        self.value = value

    def __index__(self) -> int:
        return self.value
