import math
from fractions import Fraction

import pytest

from zero_sum_sequences import (
    AdditiveSequence,
    AdditiveSequenceSpace,
    FiniteAdditiveGroup,
)

from groups import IndexableInteger, cyclic_group


@pytest.fixture
def c3_sequences():
    return AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)


def test_space_validates_parent_and_davenport_bound():
    with pytest.raises(TypeError, match="callable"):
        AdditiveSequenceSpace(object(), davenport_bound=3)
    for bad_bound in (0, -1, 1.5, True):
        with pytest.raises(ValueError, match="positive integer"):
            AdditiveSequenceSpace(cyclic_group(3), davenport_bound=bad_bound)


def test_space_accepts_indexable_integers_for_davenport_bound():
    space = AdditiveSequenceSpace(
        cyclic_group(3),
        davenport_bound=IndexableInteger(3),
    )

    assert space.davenport_bound == 3
    assert type(space.davenport_bound) is int

    for bad_bound in (IndexableInteger(0), Fraction(3, 2)):
        with pytest.raises(ValueError, match="positive integer"):
            AdditiveSequenceSpace(cyclic_group(3), davenport_bound=bad_bound)


def test_construction_is_canonical_and_retains_space(c3_sequences):
    sequence = c3_sequences(2, 1, 1)

    assert sequence.parent() is c3_sequences
    assert tuple(map(int, sequence)) == (1, 1, 2)
    assert tuple(map(int, sequence.support)) == (1, 2)
    assert sequence.multiplicity(1) == 2
    assert sequence.multiplicity(0) == 0
    assert hash(sequence) == hash(c3_sequences(1, 2, 1))


def test_sequence_arithmetic_preserves_space(c3_sequences):
    left = c3_sequences(1, 1)
    right = c3_sequences(2)

    assert left * right == c3_sequences(1, 1, 2)
    assert right ** 2 == c3_sequences(2, 2)
    assert right ** 0 == c3_sequences()
    assert (left * right).parent() is c3_sequences
    assert right.divides(left * right)
    assert left * right / right == left

    with pytest.raises(ValueError, match="does not divide"):
        right / left
    with pytest.raises(ValueError, match="non-negative"):
        right ** (-1)


def test_negation_is_termwise(c3_sequences):
    assert -c3_sequences(1, 1, 2) == c3_sequences(2, 2, 1)
    assert -c3_sequences() == c3_sequences()


@pytest.mark.parametrize(
    "expression",
    ["S + T", "S - T", "2 * S", "S * 2", "0 + S", "S ** True", "pow(S, 2, 5)"],
)
def test_additive_and_integer_operators_are_rejected(c3_sequences, expression):
    with pytest.raises(TypeError):
        eval(expression, {"S": c3_sequences(1, 2), "T": c3_sequences(1)})


def test_sequences_multiply_with_math_prod(c3_sequences):
    sequences = [c3_sequences(1), c3_sequences(2), c3_sequences(1, 1)]

    assert math.prod(sequences) == c3_sequences(1, 1, 1, 2)
    assert math.prod([], start=c3_sequences()) == c3_sequences()
    assert sequences[0].__rmul__(True) is NotImplemented


def test_exponents_accept_indexable_integers(c3_sequences):
    right = c3_sequences(2)

    assert right ** IndexableInteger(2) == c3_sequences(2, 2)
    assert right ** IndexableInteger(0) == c3_sequences()
    with pytest.raises(ValueError, match="non-negative"):
        right ** IndexableInteger(-1)
    assert right.__pow__(Fraction(3, 2)) is NotImplemented
    assert right.__pow__(1.5) is NotImplemented


def test_constructor_takes_one_argument_per_term(c3_sequences):
    assert c3_sequences(2, 1, 1) == c3_sequences.from_terms([1, 1, 2])
    assert c3_sequences.from_terms(iter([2, 1])) == c3_sequences(1, 2)
    sequence = c3_sequences(1, 2)
    assert c3_sequences(sequence) is sequence
    for old_style in ([1, 2], [], iter([1])):
        with pytest.raises(TypeError, match="separate arguments"):
            c3_sequences(old_style)
    other = AdditiveSequenceSpace(cyclic_group(3))
    with pytest.raises(TypeError, match="different space"):
        other(sequence)
    assert other.from_terms(sequence) == other(1, 2)


@pytest.mark.parametrize("max_length", [True, 1.5, -1])
def test_subsequence_length_caps_are_exact_nonnegative_integers(c3_sequences, max_length):
    with pytest.raises(ValueError, match="maximum length must be a non-negative integer"):
        next(c3_sequences(1, 2).subsequences(max_length=max_length))


def test_arithmetic_rejects_different_sequence_spaces():
    left_space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    right_space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)

    with pytest.raises(TypeError, match="different spaces"):
        left_space(1) * right_space(2)


def test_typed_empty_sequence_uses_the_base_parent_zero(c3_sequences):
    empty = c3_sequences()

    assert empty.parent() is c3_sequences
    assert empty.sum() == c3_sequences.base_parent.zero()
    assert empty.is_zero_sum()
    assert not empty.is_atom()


def test_zero_sum_and_atom_predicates(c3_sequences):
    sequence = c3_sequences

    assert sequence(0).is_atom()
    assert sequence(1, 2).is_atom()
    assert sequence(1, 1, 1).is_atom()
    assert sequence(1, 1, 2, 2).is_zero_sum()
    assert not sequence(1, 1, 2, 2).is_atom()
    assert not sequence(1, 1).is_zero_sum()


def test_subsequences_are_distinct_and_respect_filters(c3_sequences):
    sequence = c3_sequences(1, 1, 2)

    assert list(sequence.subsequences(nonempty=False)) == [
        c3_sequences(),
        c3_sequences(2),
        c3_sequences(1),
        c3_sequences(1, 2),
        c3_sequences(1, 1),
        c3_sequences(1, 1, 2),
    ]
    assert list(sequence.subsequences(proper=True, max_length=1)) == [
        c3_sequences(2),
        c3_sequences(1),
    ]


def test_display_uses_multiplicative_notation(c3_sequences):
    assert str(c3_sequences()) == "1"
    assert str(c3_sequences(1, 1, 2)) == "(1)^2 · 2"
    assert repr(c3_sequences(1, 2)) == "1 · 2"
    group = FiniteAdditiveGroup.cyclic_product(2, 4)
    assert repr(group) == "FiniteAdditiveGroup.cyclic_product(2, 4)"
    assert repr(AdditiveSequenceSpace(group)) == (
        "AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(2, 4))"
    )
    assert repr(AdditiveSequenceSpace(group, davenport_bound=5)).endswith(
        ", davenport_bound=5)"
    )
    assert repr(c3_sequences.base_parent) == "FiniteAdditiveGroup(order=3, zero=0)"


def test_display_does_not_double_parenthesize_repeated_terms():
    space = AdditiveSequenceSpace(
        FiniteAdditiveGroup.cyclic_product(2, 4), davenport_bound=5
    )

    assert str(space((0, 1), (0, 3))) == "(0, 1) \u00b7 (0, 3)"
    assert str(space((0, 2), (0, 2))) == "(0, 2)^2"


def test_multiplicity_validation(c3_sequences):
    with pytest.raises(ValueError, match="multiplicity must be a non-negative integer"):
        c3_sequences.from_multiplicities({c3_sequences.base_parent(1): True})


def test_multiplicities_accept_indexable_integers(c3_sequences):
    one = c3_sequences.base_parent(1)

    assert c3_sequences.from_multiplicities(
        {one: IndexableInteger(2)}
    ) == c3_sequences(1, 1)
    with pytest.raises(ValueError, match="multiplicity must be a non-negative integer"):
        c3_sequences.from_multiplicities({one: Fraction(3, 2)})


def test_direct_sequence_construction_requires_a_space():
    with pytest.raises(TypeError, match="AdditiveSequenceSpace"):
        AdditiveSequence([1, 2])


def test_large_multiplicity_stays_compressed(monkeypatch):
    class Parent:
        calls = 0

        def __call__(self, term):
            self.calls += 1
            return int(term)

        def zero(self):
            return 0

    parent = Parent()
    space = AdditiveSequenceSpace(parent)

    def forbidden(*args):
        raise AssertionError("multiplicity construction must not expand terms")

    monkeypatch.setattr(AdditiveSequence, "__init__", forbidden)
    sequence = space.from_multiplicities({1: 100_000})
    assert parent.calls == 1
    assert sequence.multiplicities == {1: 100_000}
    assert len(sequence) == 100_000


def test_aliases_are_merged_and_zero_count_keys_are_not_coerced():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    result = space.from_multiplicities({(1,): 2, (4,): 3, (2,): 1, object(): 0})
    expected = space((1,)) ** 5 * space((2,))
    assert result == expected
    assert hash(result) == hash(expected)
    assert result.multiplicities == {(1,): 5, (2,): 1}
    assert result.parent() is space
    assert space.from_multiplicities({object(): 0}) == space()


def test_nonorderable_terms_keep_the_construction_error():
    class Parent:
        def __call__(self, term):
            return term

        def zero(self):
            return 0

    with pytest.raises(TypeError, match="mutually orderable"):
        AdditiveSequenceSpace(Parent()).from_multiplicities({1: 1, "a": 1})
