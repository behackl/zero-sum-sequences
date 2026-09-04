from dataclasses import FrozenInstanceError

import pytest

from groups import cyclic_group
from zero_sum_sequences import (
    AdditiveSequenceSpace,
    FactorizationRelation,
)


C3 = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
A = C3([1, 1, 1])
B = C3([2, 2, 2])
P = C3([1, 2])


def test_relation_canonicalizes_cancels_and_computes_distance():
    relation = FactorizationRelation(
        source=(B, P, A),
        target=(P, P, P, P),
    )

    assert relation.source == (P, A, B)
    assert relation.target == (P, P, P, P)
    assert relation.product == A + B + P
    assert relation.common_factor == (P,)
    assert relation.reduced_source == (A, B)
    assert relation.reduced_target == (P, P, P)
    assert not relation.is_reduced
    assert relation.distance == 3

    reduced = relation.reduced()
    assert reduced.source == (A, B)
    assert reduced.target == (P, P, P)
    assert reduced.product == A + B
    assert reduced.is_reduced
    assert reduced.distance == 3


def test_relation_cancels_repeated_common_factors():
    relation = FactorizationRelation(
        source=(P, A, P, B),
        target=(P, P, P, P, P),
    )

    assert relation.common_factor == (P, P)
    assert relation.reduced_source == (A, B)
    assert relation.reduced_target == (P, P, P)


def test_reduced_relation_is_unchanged_and_orientation_is_preserved():
    relation = FactorizationRelation((B, A), (P, P, P))
    reordered = FactorizationRelation((A, B), (P, P, P))
    reversed_relation = FactorizationRelation((P, P, P), (A, B))

    assert relation == reordered
    assert hash(relation) == hash(reordered)
    assert relation != reversed_relation
    assert relation.is_reduced
    assert relation.reduced() is relation


def test_identity_relation_reduces_to_the_empty_relation():
    relation = FactorizationRelation((A, P), (P, A))

    assert relation.common_factor == (P, A)
    assert relation.distance == 0

    reduced = relation.reduced()
    assert reduced.source == ()
    assert reduced.target == ()
    assert reduced.product == C3()
    assert reduced.is_reduced
    assert reduced.reduced() is reduced


def test_empty_relation_requires_an_explicit_space():
    with pytest.raises(ValueError, match="cannot infer"):
        FactorizationRelation((), ())

    relation = FactorizationRelation((), (), space=C3)

    assert relation.product == C3()
    assert relation.is_reduced
    assert relation.distance == 0


def test_relation_rejects_unequal_products_and_incompatible_spaces():
    with pytest.raises(ValueError, match="different products"):
        FactorizationRelation((A, B), (P, P))

    other_space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    other_a = other_space([1, 1, 1])
    with pytest.raises(TypeError, match="different sequence spaces"):
        FactorizationRelation((A,), (other_a,))


def test_relation_rejects_malformed_inputs():
    with pytest.raises(TypeError, match="source factorization must be iterable"):
        FactorizationRelation(None, ())
    with pytest.raises(TypeError, match="target factors"):
        FactorizationRelation((A,), (object(),))
    with pytest.raises(TypeError, match="space must be"):
        FactorizationRelation((), (), space=object())


def test_relation_treats_factorization_entries_as_trusted():
    composite = A + B
    assert not composite.is_atom()

    relation = FactorizationRelation((composite,), (composite,))

    assert relation.product == composite
    assert relation.distance == 0


def test_relation_is_immutable():
    relation = FactorizationRelation((A, B), (P, P, P))

    with pytest.raises(FrozenInstanceError):
        relation.source = ()
