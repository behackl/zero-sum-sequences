import pytest
from sage.all import GF, Integer, QQ, Zmod

from zero_sum_sequences import AdditiveSequence, AdditiveSequenceSpace


@pytest.fixture
def c3_sequences():
    return AdditiveSequenceSpace(Zmod(3), davenport_bound=3)


def test_space_validates_parent_and_davenport_bound():
    with pytest.raises(TypeError, match="callable"):
        AdditiveSequenceSpace(object(), davenport_bound=3)
    for bad_bound in (0, -1, 1.5, True):
        with pytest.raises(ValueError, match="positive integer"):
            AdditiveSequenceSpace(Zmod(3), davenport_bound=bad_bound)


def test_space_accepts_sage_exact_integers_for_davenport_bound():
    space = AdditiveSequenceSpace(Zmod(3), davenport_bound=Integer(3))

    assert space.davenport_bound == 3
    assert type(space.davenport_bound) is int

    for bad_bound in (Integer(0), QQ(3) / 2):
        with pytest.raises(ValueError, match="positive integer"):
            AdditiveSequenceSpace(Zmod(3), davenport_bound=bad_bound)


def test_construction_is_canonical_and_retains_space(c3_sequences):
    sequence = c3_sequences([2, 1, 1])

    assert sequence.parent() is c3_sequences
    assert tuple(map(int, sequence)) == (1, 1, 2)
    assert tuple(map(int, sequence.support)) == (1, 2)
    assert sequence.multiplicity(1) == 2
    assert sequence.multiplicity(0) == 0
    assert hash(sequence) == hash(c3_sequences([1, 2, 1]))


def test_mutable_sage_terms_are_copied_before_becoming_immutable():
    space = AdditiveSequenceSpace(GF(3) ** 1, davenport_bound=3)
    term = space.base_parent([1])
    sequence = space([term])

    assert term.is_mutable()
    assert not next(iter(sequence)).is_mutable()
    assert sequence.multiplicity(term) == 1


def test_sequence_arithmetic_preserves_space(c3_sequences):
    left = c3_sequences([1, 1])
    right = c3_sequences([2])

    assert left + right == c3_sequences([1, 1, 2])
    assert 2 * right == c3_sequences([2, 2])
    assert 0 * right == c3_sequences()
    assert (left + right).parent() is c3_sequences
    assert right.divides(left + right)
    assert left + right - right == left

    with pytest.raises(ValueError, match="not a subsequence"):
        right - left
    with pytest.raises(ValueError, match="non-negative"):
        -1 * right


def test_sequence_arithmetic_accepts_sage_integer_repetitions(c3_sequences):
    right = c3_sequences([2])

    assert Integer(2) * right == c3_sequences([2, 2])
    assert right * Integer(0) == c3_sequences()
    with pytest.raises(ValueError, match="non-negative"):
        Integer(-1) * right
    assert right.__mul__(QQ(3) / 2) is NotImplemented
    assert right.__mul__(1.5) is NotImplemented


def test_arithmetic_rejects_different_sequence_spaces():
    left_space = AdditiveSequenceSpace(Zmod(3), davenport_bound=3)
    right_space = AdditiveSequenceSpace(Zmod(3), davenport_bound=3)

    with pytest.raises(TypeError, match="different spaces"):
        left_space([1]) + right_space([2])


def test_typed_empty_sequence_uses_the_base_parent_zero(c3_sequences):
    empty = c3_sequences()

    assert empty.parent() is c3_sequences
    assert empty.total() == c3_sequences.base_parent.zero()
    assert empty.is_zero_sum()
    assert not empty.is_atom()


def test_zero_sum_and_atom_predicates(c3_sequences):
    sequence = c3_sequences

    assert sequence([0]).is_atom()
    assert sequence([1, 2]).is_atom()
    assert sequence([1, 1, 1]).is_atom()
    assert sequence([1, 1, 2, 2]).is_zero_sum()
    assert not sequence([1, 1, 2, 2]).is_atom()
    assert not sequence([1, 1]).is_zero_sum()


def test_subsequences_are_distinct_and_respect_filters(c3_sequences):
    sequence = c3_sequences([1, 1, 2])

    assert list(sequence.subsequences(nonempty=False)) == [
        c3_sequences(),
        c3_sequences([2]),
        c3_sequences([1]),
        c3_sequences([1, 2]),
        c3_sequences([1, 1]),
        c3_sequences([1, 1, 2]),
    ]
    assert list(sequence.subsequences(proper=True, max_length=1)) == [
        c3_sequences([2]),
        c3_sequences([1]),
    ]


def test_display_is_unambiguous(c3_sequences):
    assert str(c3_sequences()) == "1"
    assert str(c3_sequences([1, 1, 2])) == "(1)^2 · 2"
    assert "AdditiveSequenceSpace" in repr(c3_sequences([1, 2]))


@pytest.mark.parametrize("bad_count", [-1, 1.5, True])
def test_multiplicity_validation(c3_sequences, bad_count):
    with pytest.raises(ValueError, match="non-negative integers"):
        c3_sequences.from_multiplicities({c3_sequences.base_parent(1): bad_count})


def test_multiplicities_accept_sage_integers_but_not_nonintegral_values(c3_sequences):
    one = c3_sequences.base_parent(1)

    assert c3_sequences.from_multiplicities({one: Integer(2)}) == c3_sequences([1, 1])
    with pytest.raises(ValueError, match="non-negative integers"):
        c3_sequences.from_multiplicities({one: QQ(3) / 2})


def test_direct_sequence_construction_requires_a_space():
    with pytest.raises(TypeError, match="AdditiveSequenceSpace"):
        AdditiveSequence([1, 2])
