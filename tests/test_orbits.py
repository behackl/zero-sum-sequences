import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AutomorphismAction,
    AutomorphismActionUnavailable,
    FiniteAdditiveGroup,
    OrbitWitness,
)
from zero_sum_sequences.orbits import is_in_same_orbit, orbit, orbit_witness

from groups import cyclic_group, product_of_cyclic_groups


@pytest.fixture
def c3_sequences():
    return AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)


def test_map_terms_preserves_multiplicities_and_canonical_order(c3_sequences):
    sequence = c3_sequences([1, 1, 2])

    assert sequence.map_terms(lambda term: -term) == c3_sequences([1, 2, 2])
    assert sequence.map_terms(lambda term: 0) == c3_sequences([0, 0, 0])
    assert sequence.map_terms(lambda term: -term).parent() is c3_sequences


def test_map_terms_can_change_sequence_space(c3_sequences):
    target = AdditiveSequenceSpace(product_of_cyclic_groups(2), davenport_bound=2)

    image = c3_sequences([1, 2]).map_terms(
        lambda term: (term % 2,),
        target_space=target,
    )

    assert image == target([(1,), (0,)])
    assert image.parent() is target


def test_map_terms_validates_mapping_and_target_space(c3_sequences):
    with pytest.raises(TypeError, match="mapping must be callable"):
        c3_sequences([1]).map_terms(1)
    with pytest.raises(TypeError, match="target_space"):
        c3_sequences([1]).map_terms(lambda term: term, target_space=object())


def test_callable_generators_produce_deterministic_orbit_and_witnesses(
    c3_sequences,
):
    # Negation generates Aut(C_3).  The sequence is unordered, so map_terms
    # also verifies that multiplicities survive the induced action.
    negation = lambda term: -term
    action = AutomorphismAction([negation])
    sequence = c3_sequences([1, 1, 2])
    image = c3_sequences([1, 2, 2])

    assert orbit(sequence, action=action) == (sequence, image)
    assert sequence.orbit(action=action) == (sequence, image)
    assert is_in_same_orbit(sequence, image, action=action)
    assert sequence.is_in_same_orbit(image, action=action)
    witness = orbit_witness(sequence, image, action=action)
    assert witness == OrbitWitness((0,))
    assert action.apply_word(sequence, witness) == image
    assert sequence.orbit_witness(image, action=action) == witness


def test_matrix_style_generators_use_apply_term_callback():
    class Shift:
        def __init__(self, amount):
            self.amount = amount

    group = cyclic_group(3)
    space = AdditiveSequenceSpace(group, davenport_bound=3)
    action = AutomorphismAction(
        [Shift(1)],
        apply_term=lambda generator, term: term + generator.amount,
    )
    sequence = space([0, 1])

    assert orbit(sequence, action=action) == (
        space([0, 1]),
        space([1, 2]),
        space([0, 2]),
    )


def test_empty_action_has_singleton_orbit(c3_sequences):
    sequence = c3_sequences([1, 2])

    assert orbit(sequence, action=AutomorphismAction(())) == (sequence,)
    assert sequence.orbit_witness(sequence, action=()) == OrbitWitness()
    assert not is_in_same_orbit(
        sequence,
        c3_sequences([1]),
        action=AutomorphismAction(()),
    )


def test_orbit_witness_is_none_for_different_orbits(c3_sequences):
    action = AutomorphismAction([lambda term: -term])

    assert orbit_witness(
        c3_sequences([1]),
        c3_sequences([1, 1]),
        action=action,
    ) is None


def test_orbit_reports_when_action_cannot_be_discovered(c3_sequences):
    with pytest.raises(AutomorphismActionUnavailable, match="pass action="):
        c3_sequences([1]).orbit()


def test_parent_can_supply_generators_lazily():
    group = FiniteAdditiveGroup(
        range(3),
        zero=0,
        add=lambda left, right: (left + right) % 3,
        coerce=lambda value: int(value) % 3,
        automorphism_generators=(lambda term: -term,),
    )
    space = AdditiveSequenceSpace(group, davenport_bound=3)

    sequence = space([1, 1, 2])
    assert sequence.orbit() == (sequence, space([1, 2, 2]))
    assert space._automorphism_action is not None


def test_action_can_materialize_a_witness(c3_sequences):
    action = AutomorphismAction(
        [lambda term: -term],
        materialize_word=lambda word: tuple(word),
    )
    witness = OrbitWitness((0, 0))

    assert action.materialize_word(witness) == (0, 0)
    with pytest.raises(TypeError, match="cannot materialize"):
        AutomorphismAction([lambda term: term]).materialize_word(witness)


def test_finite_cas_parent_can_supply_an_automorphism_group():
    class GeneratedGroup:
        def gens(self):
            return (lambda term: -term,)

    class Parent:
        def __call__(self, value):
            return int(value) % 3

        def zero(self):
            return 0

        def is_finite(self):
            return True

        def automorphism_group(self):
            return GeneratedGroup()

    space = AdditiveSequenceSpace(Parent(), davenport_bound=3)
    sequence = space([1, 1, 2])

    assert sequence.orbit() == (sequence, space([1, 2, 2]))


def test_orbits_reject_different_sequence_spaces(c3_sequences):
    other_space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    with pytest.raises(TypeError, match="different spaces"):
        is_in_same_orbit(
            c3_sequences([1]),
            other_space([1]),
            action=[lambda term: -term],
        )
