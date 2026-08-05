"""Integration tests over structurally different Sage additive parents."""

from collections import Counter

from sage.all import AdditiveAbelianGroup, GF

from zero_sum_sequences import AdditiveSequenceSpace


def assert_factorization_profiles(sequence, expected_profiles):
    """Check complete factorizations by their multisets of atom lengths."""

    factorizations = list(sequence.factorizations())
    actual_profiles = Counter(
        tuple(sorted(map(len, factorization)))
        for factorization in factorizations
    )

    assert actual_profiles == Counter(expected_profiles)
    assert sequence.length_set() == {
        len(profile) for profile in expected_profiles
    }
    assert all(
        sum(factorization, sequence.parent()()) == sequence
        for factorization in factorizations
    )

    witnesses = sequence.factorization_witnesses()
    assert set(witnesses) == sequence.length_set()
    assert all(
        len(factorization) == length
        and sum(factorization, sequence.parent()()) == sequence
        for length, factorization in witnesses.items()
    )


def test_klein_four_group_has_two_distinct_factorization_lengths():
    group = GF(2) ** 2
    space = AdditiveSequenceSpace(group, davenport_bound=3)
    e1, e2 = group.basis()
    maximal_atom = space([e1, e2, e1 + e2])
    sequence = 2 * maximal_atom

    assert maximal_atom.is_atom()
    assert_factorization_profiles(
        sequence,
        [
            (3, 3),
            (2, 2, 2),
        ],
    )


def test_mixed_exponent_additive_abelian_group_is_supported():
    group = AdditiveAbelianGroup([2, 4])
    space = AdditiveSequenceSpace(group, davenport_bound=5)
    a, b, c = group((1, 0)), group((0, 1)), group((1, 3))
    atom = space([a, b, c])
    inverse_atom = space(-term for term in atom)
    sequence = atom + inverse_atom

    assert atom.is_atom()
    assert inverse_atom.is_atom()
    assert_factorization_profiles(
        sequence,
        [
            (3, 3),
            (2, 2, 2),
        ],
    )


def test_rank_two_ternary_vector_space_has_richer_factorizations():
    group = GF(3) ** 2
    space = AdditiveSequenceSpace(group, davenport_bound=5)
    e1, e2 = group.basis()
    maximal_atom = space([e1, e1, e2, e2, e1 + e2])
    inverse_atom = space(-term for term in maximal_atom)
    sequence = maximal_atom + inverse_atom

    assert maximal_atom.is_atom()
    assert inverse_atom.is_atom()
    assert_factorization_profiles(
        sequence,
        [
            (5, 5),
            (2, 4, 4),
            (2, 4, 4),
            (2, 2, 3, 3),
            (2, 2, 2, 2, 2),
        ],
    )
