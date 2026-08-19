from collections import deque

import pytest

from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup

from groups import cyclic_group, product_of_cyclic_groups


def test_finite_group_coerces_and_adds_canonical_elements():
    group = cyclic_group(4)

    assert tuple(group) == (0, 1, 2, 3)
    assert len(group) == 4
    assert group.zero() == 0
    assert group(5) == 1
    assert group.add(3, 3) == 2
    assert group.is_finite()


def test_parent_operation_supports_plain_tuple_elements():
    group = product_of_cyclic_groups(2, 4)
    space = AdditiveSequenceSpace(group, davenport_bound=5)

    sequence = space([(1, 1), (1, 3)])

    assert sequence.total() == (0, 0)
    assert sequence.is_zero_sum()
    assert sequence.is_atom()


def test_finite_group_validates_its_presentation():
    with pytest.raises(ValueError, match="distinct"):
        FiniteAdditiveGroup(
            [0, 1, 2],
            zero=0,
            add=lambda left, right: left + right,
            coerce=lambda value: int(value) % 2,
        )
    with pytest.raises(ValueError, match="not an element"):
        FiniteAdditiveGroup(
            [1, 2],
            zero=0,
            add=lambda left, right: left + right,
        )


def test_finite_group_checks_operation_results_when_used():
    group = FiniteAdditiveGroup(
        [0, 1],
        zero=0,
        add=lambda left, right: left + right,
    )

    with pytest.raises(ValueError, match="not an element"):
        group.add(1, 1)


def generated_permutations(group):
    """Return the permutation group generated on the finite parent."""

    elements = tuple(group)
    identity = tuple(range(len(elements)))
    generator_permutations = tuple(
        tuple(elements.index(generator(element)) for element in elements)
        for generator in group.automorphism_generators()
    )
    permutations = {identity}
    pending = deque((identity,))
    while pending:
        current = pending.popleft()
        for generator in generator_permutations:
            image = tuple(generator[position] for position in current)
            if image not in permutations:
                permutations.add(image)
                pending.append(image)
    return permutations


def test_cyclic_product_constructs_coordinate_group():
    group = FiniteAdditiveGroup.cyclic_product(2, 4)

    assert len(group) == 8
    assert group.zero() == (0, 0)
    assert group((3, 5)) == (1, 1)
    assert group.add((1, 3), (1, 2)) == (0, 1)


def test_cyclic_product_supplies_full_automorphism_generators():
    group = FiniteAdditiveGroup.cyclic_product(2, 4)
    generators = group.automorphism_generators()

    assert len(generated_permutations(group)) == 8
    assert all(generator(group.zero()) == group.zero() for generator in generators)
    assert all(
        group.add(generator(left), generator(right))
        == generator(group.add(left, right))
        for generator in generators
        for left in group
        for right in group
    )


@pytest.mark.parametrize(
    ("moduli", "expected_order"),
    [
        ((4,), 2),
        ((2, 2), 6),
        ((3, 3), 48),
        ((2, 2, 2), 168),
        ((6, 10), 48),
    ],
)
def test_cyclic_product_generators_have_expected_group_order(
    moduli,
    expected_order,
):
    group = FiniteAdditiveGroup.cyclic_product(*moduli)

    assert len(generated_permutations(group)) == expected_order


@pytest.mark.parametrize("moduli", [(), (0,), (-2,), (True,), (2.5,)])
def test_cyclic_product_validates_moduli(moduli):
    with pytest.raises(ValueError, match="modul|at least one"):
        FiniteAdditiveGroup.cyclic_product(*moduli)


def test_finite_group_validates_automorphism_generators():
    with pytest.raises(TypeError, match="automorphism generators"):
        FiniteAdditiveGroup(
            [0],
            zero=0,
            add=lambda left, right: 0,
            automorphism_generators=(object(),),
        )
