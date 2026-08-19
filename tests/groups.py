"""Small dependency-free additive parents used by the core test suite."""

from zero_sum_sequences import FiniteAdditiveGroup


def cyclic_group(modulus: int) -> FiniteAdditiveGroup[int]:
    return FiniteAdditiveGroup(
        range(modulus),
        zero=0,
        add=lambda left, right: (left + right) % modulus,
        coerce=lambda value: int(value) % modulus,
    )


def product_of_cyclic_groups(*moduli: int) -> FiniteAdditiveGroup[tuple[int, ...]]:
    return FiniteAdditiveGroup.cyclic_product(*moduli)


class IndexableInteger:
    def __init__(self, value: int) -> None:
        self.value = value

    def __index__(self) -> int:
        return self.value
