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
