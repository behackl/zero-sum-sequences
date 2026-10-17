from math import gcd, inf, lcm

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AutomorphismGroup,
    FiniteAdditiveGroup,
    ZeroSumConstant,
    davenport_constant,
    erdos_ginzburg_ziv_constant,
    eta_constant,
    exponent,
    zero_sum_constant,
)

from groups import cyclic_group


def space_of(*moduli):
    parent = FiniteAdditiveGroup.cyclic_product(*moduli)
    order = 1
    for modulus in moduli:
        order *= modulus
    return AdditiveSequenceSpace(parent, davenport_bound=order)


def is_L_free(sequence, lengths):
    return not any(
        len(subsequence) in lengths and subsequence.is_zero_sum()
        for subsequence in sequence.subsequences()
    )


def maximal_L_free_sequences(space, lengths):
    """All L-free sequences of maximal length, by exhaustive search.

    Subsequences of L-free sequences are L-free, so every L-free sequence
    extends a shorter one by a single term.
    """

    elements = tuple(space.base_parent)
    free = {space()}
    while True:
        extended = {
            sequence + space([term])
            for sequence in free
            for term in elements
        }
        extended = {sequence for sequence in extended if is_L_free(sequence, lengths)}
        if not extended:
            return sorted(free)
        free = extended


def trivial_group(space):
    group = space.automorphism_group()
    return group.subgroup([group.identity])


# agreement with exhaustive search


# Smaller groups are covered by the independent coordinate oracle.
@pytest.mark.parametrize(
    "moduli, lengths",
    [
        pytest.param((5,), {5}, id="C5-exponent"),
        pytest.param((6,), range(1, 7), id="C6-at-most-exponent"),
        pytest.param((6,), {2, 6}, id="C6-sparse"),
        pytest.param((2, 4), {4}, id="C2xC4-exponent"),
        pytest.param((2, 4), range(1, 5), id="C2xC4-at-most-exponent"),
        pytest.param((2, 4), {3, 4}, id="C2xC4-near-exponent"),
        pytest.param((2, 4), {2, 4}, id="C2xC4-sparse"),
        pytest.param((2, 4), {4, 8}, id="C2xC4-longer-multiple"),
        pytest.param((3, 3), {3}, id="C3xC3-exponent"),
        pytest.param((3, 3), range(1, 4), id="C3xC3-at-most-exponent"),
        pytest.param((3, 3), {2, 3}, id="C3xC3-near-exponent"),
        pytest.param((3, 3), {3, 6}, id="C3xC3-longer-multiple"),
    ],
)
def test_value_and_extremal_orbits_match_exhaustive_search(moduli, lengths):
    space = space_of(*moduli)
    group = space.automorphism_group()
    lengths = set(lengths)
    maximal = maximal_L_free_sequences(space, lengths)

    result = zero_sum_constant(space, lengths)

    assert result.value == len(maximal[0]) + 1
    assert result.extremal == tuple(sorted({group.canonical(s) for s in maximal}))


def test_trivial_group_yields_every_extremal_sequence():
    space = space_of(2, 4)
    lengths = {2, 4}

    result = zero_sum_constant(space, lengths, symmetry=trivial_group(space))

    assert result.value == zero_sum_constant(space, lengths).value
    assert result.extremal == tuple(sorted(maximal_L_free_sequences(space, lengths)))


def test_parent_without_automorphisms_is_enumerated_without_reduction():
    space = AdditiveSequenceSpace(cyclic_group(5), davenport_bound=5)

    result = zero_sum_constant(space, {5})

    assert result.value == 9
    assert result.extremal == (
        space([0, 0, 0, 0, 1, 1, 1, 1]),
        space([0, 0, 0, 0, 2, 2, 2, 2]),
        space([0, 0, 0, 0, 3, 3, 3, 3]),
        space([0, 0, 0, 0, 4, 4, 4, 4]),
        *(
            space([a] * 4 + [b] * 4)
            for a, b in [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
        ),
    )


def test_result_records_the_normalized_lengths():
    space = space_of(3)

    # zero is excluded (length 1), and 1*1*2*2 has zero sums only of length 2
    assert zero_sum_constant(space, [3, 3, 1]) == ZeroSumConstant(
        frozenset({1, 3}), 5, (space([(1,), (1,), (2,), (2,)]),)
    )


# closed formulas for cyclic groups


@pytest.mark.parametrize("n", range(2, 10))
def test_erdos_ginzburg_ziv_theorem(n):
    assert erdos_ginzburg_ziv_constant(space_of(n)) == 2 * n - 1


@pytest.mark.parametrize("n", range(2, 10))
def test_two_lengths_near_the_exponent(n):
    # Guillot, Schmid, Zerdoum: s_{n-j, n}(C_n) is 2n - j - 1 if j | n,
    # and 2n - 1 otherwise
    space = space_of(n)
    for j in range(1, n):
        expected = 2 * n - j - 1 if n % j == 0 else 2 * n - 1
        assert zero_sum_constant(space, {n - j, n}).value == expected


@pytest.mark.parametrize("n", range(2, 9))
@pytest.mark.parametrize("j", range(2, 5))
def test_multiples_of_a_fixed_length(n, j):
    # s_{jN}(C_n) = lcm(n, j) + gcd(n, j) - 1
    expected = lcm(n, j) + gcd(n, j) - 1
    lengths = range(j, expected + 1, j)

    assert zero_sum_constant(space_of(n), lengths).value == expected


@pytest.mark.parametrize("moduli", [(5,), (2, 4), (3, 3)])
def test_value_is_infinite_without_a_multiple_of_the_exponent(moduli):
    space = space_of(*moduli)
    e = exponent(space)

    result = zero_sum_constant(space, [k for k in range(1, 3 * e) if k % e])

    assert result.value == inf
    assert result.extremal == ()


# published values for groups of higher rank


@pytest.mark.parametrize(
    "moduli, lengths, value",
    [
        ((3, 6), range(1, 7), 10),
        ((3, 6), range(1, 8), 9),
        ((2, 2, 4), range(1, 5), 8),
        ((2, 2, 4), range(1, 6), 7),
        ((2, 2, 2, 2), range(1, 3), 16),
        ((2, 2, 2, 2), range(1, 4), 9),
        ((2, 2, 2, 2), range(1, 5), 6),
        ((2, 2, 2, 2), {2}, 17),
        ((3, 3, 3), range(1, 5), 10),
        ((3, 3, 3), range(1, 6), 9),
        ((3, 3, 3), range(1, 7), 8),
    ],
)
def test_published_values(moduli, lengths, value):
    assert zero_sum_constant(space_of(*moduli), lengths).value == value


@pytest.mark.parametrize(
    "moduli, D, eta, s",
    [((2, 2), 3, 4, 5), ((2, 4), 5, 6, 9), ((3, 3), 5, 7, 9), ((2, 2, 2), 4, 8, 9)],
)
def test_classical_constants(moduli, D, eta, s):
    space = space_of(*moduli)

    assert davenport_constant(space) == D
    assert eta_constant(space) == eta
    assert erdos_ginzburg_ziv_constant(space) == s


@pytest.mark.parametrize(
    "moduli, bound", [((2, 2), 4), ((2, 4), 6), ((3, 3), 6), ((2, 2, 2), 5), ((2, 6), 8)]
)
def test_davenport_constant_is_the_maximal_atom_length(moduli, bound):
    parent = FiniteAdditiveGroup.cyclic_product(*moduli)
    space = AdditiveSequenceSpace(parent, davenport_bound=bound)
    catalogue = space.enumerate_atom_catalogue()

    assert davenport_constant(space) == max(len(atom) for atom in catalogue)


@pytest.mark.parametrize("moduli, value", [((6,), 6), ((2, 4), 4), ((2, 3), 6), ((2, 2, 2), 2)])
def test_exponent(moduli, value):
    assert exponent(space_of(*moduli)) == value


# validation


@pytest.mark.parametrize("lengths", [[], [0, 2], [-1], [True], [1.5]])
def test_rejects_invalid_lengths(lengths):
    with pytest.raises(ValueError, match="lengths"):
        zero_sum_constant(space_of(3), lengths)


def test_rejects_non_iterable_lengths():
    with pytest.raises(TypeError, match="lengths"):
        zero_sum_constant(space_of(3), 3)


@pytest.mark.parametrize("bad", [None, "auto", 0, 1, [], object()])
def test_rejects_invalid_symmetry_even_for_an_infinite_result(bad):
    with pytest.raises(TypeError, match="symmetry"):
        zero_sum_constant(space_of(3), {2}, symmetry=bad)


@pytest.mark.parametrize(
    "query", [davenport_constant, eta_constant, erdos_ginzburg_ziv_constant]
)
@pytest.mark.parametrize("warm", [False, True], ids=["cold", "cached"])
def test_scalar_wrappers_validate_symmetry_before_fast_paths(query, warm):
    space = space_of(3)
    if warm:
        query(space)
    with pytest.raises(TypeError, match="symmetry"):
        query(space, symmetry=None)


def test_disabled_symmetry_does_not_discover_automorphisms(monkeypatch):
    def forbidden(*args):
        raise AssertionError("automorphisms must not be discovered")

    monkeypatch.setattr(AdditiveSequenceSpace, "automorphism_group", forbidden)
    space = space_of(2, 4)
    result = zero_sum_constant(space, {2, 4}, symmetry=False)
    assert result.extremal == tuple(sorted(maximal_L_free_sequences(space, {2, 4})))
    assert davenport_constant(space, symmetry=False) == 5
    assert eta_constant(space, symmetry=False) == 6
    assert erdos_ginzburg_ziv_constant(space, symmetry=False) == 9


def test_discovery_errors_are_not_suppressed(monkeypatch):
    def broken(*args):
        raise RuntimeError("broken provider")

    monkeypatch.setattr(AdditiveSequenceSpace, "automorphism_group", broken)
    with pytest.raises(RuntimeError, match="broken provider"):
        zero_sum_constant(space_of(3), {3})


def test_explicit_proper_subgroup_controls_extremal_representatives():
    space = space_of(2, 4)
    full = space.automorphism_group()
    generator = full.generators[0]
    members = {full.identity}
    current = generator
    while current not in members:
        members.add(current)
        current = full.compose(generator, current)
    subgroup = full.subgroup(members)
    assert 1 < len(subgroup) < len(full)
    lengths = {2, 4}
    result = zero_sum_constant(space, lengths, symmetry=subgroup)
    maximal = maximal_L_free_sequences(space, lengths)
    assert result.value == zero_sum_constant(space, lengths).value
    assert result.extremal == tuple(sorted({subgroup.canonical(s) for s in maximal}))


def test_explicit_symmetry_does_not_discover_automorphisms(monkeypatch):
    space = space_of(3)
    group = space.automorphism_group()

    def forbidden(*args):
        raise AssertionError("explicit symmetry must be used")

    monkeypatch.setattr(AdditiveSequenceSpace, "automorphism_group", forbidden)
    assert zero_sum_constant(space, {3}, symmetry=group).value == 5


def test_incompatible_symmetry_has_a_clear_error():
    group = AutomorphismGroup([0], identity=0, apply_term=lambda auto, g: (99,))
    with pytest.raises(ValueError, match="ambient parent"):
        zero_sum_constant(space_of(3), {3}, symmetry=group)


def test_symmetry_must_be_a_zero_fixing_bijection():
    group = AutomorphismGroup([0], identity=0, apply_term=lambda auto, g: (0,))
    with pytest.raises(ValueError, match="bijections"):
        zero_sum_constant(space_of(3), {3}, symmetry=group)


def test_large_sparse_lengths_do_not_allocate_large_bitsets():
    space = space_of(2)
    result = zero_sum_constant(space, {2, 10**12}, symmetry=False)
    assert result.value == 3
    assert result.extremal == zero_sum_constant(space, {2}, symmetry=False).extremal
    assert result.lengths == frozenset({2, 10**12})


def test_constants_ignore_an_insufficient_configured_bound():
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=1)
    assert davenport_constant(space, symmetry=False) == 3
    assert zero_sum_constant(space, {3}, symmetry=False).value == 5
    assert space.davenport_bound == 1  # explicit values are never replaced


def test_trivial_group_constants_and_extremal_sequences():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(1))
    assert space.davenport_bound == 1
    assert exponent(space) == 1
    assert davenport_constant(space) == eta_constant(space) == erdos_ginzburg_ziv_constant(space) == 1
    assert zero_sum_constant(space, {1}).extremal == (space(),)
    result = zero_sum_constant(space, {3}, symmetry=False)
    assert result.value == 3
    assert result.extremal == (space([(0,), (0,)]),)
