"""Orbit-action integration tests for optional Sage parents."""

import pytest

pytest.importorskip("sage.all")

from sage.all import GF, QQ

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AutomorphismActionUnavailable,
    OrbitWitness,
)
from zero_sum_sequences.orbits import automorphism_action


def test_finite_vector_space_action_is_discovered_and_cached():
    group = GF(3) ** 2
    space = AdditiveSequenceSpace(group, davenport_bound=5)
    nonzero = space([group.basis()[0]])
    zero = space([group.zero()])

    assert len(nonzero.orbit()) == 8
    assert zero.orbit() == (zero,)
    assert automorphism_action(space) is automorphism_action(space)


def test_sage_matrix_word_materialization_has_application_order():
    group = GF(3) ** 2
    space = AdditiveSequenceSpace(group, davenport_bound=5)
    action = automorphism_action(space)
    sequence = space(group.basis())
    word = OrbitWitness((0, 1))

    matrix = action.materialize_word(word)
    assert action.apply_word(sequence, word) == sequence.map_terms(
        lambda term: group(matrix * term)
    )


def test_infinite_vector_space_requires_an_explicit_action():
    group = QQ ** 2
    space = AdditiveSequenceSpace(group, davenport_bound=3)

    with pytest.raises(AutomorphismActionUnavailable, match="pass action="):
        space([group.basis()[0]]).orbit()


def test_zero_dimensional_vector_space_has_trivial_action():
    group = GF(3) ** 0
    space = AdditiveSequenceSpace(group, davenport_bound=1)
    sequence = space([group.zero()])

    assert sequence.orbit() == (sequence,)
    assert action_matrix_shape(space) == (0, 0)


def action_matrix_shape(space):
    identity = automorphism_action(space).materialize_word(OrbitWitness())
    return identity.nrows(), identity.ncols()
