import sys

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AutomorphismAction,
    AutomorphismActionUnavailable,
    AutomorphismGroupUnavailable,
    FiniteAdditiveGroup,
    OrbitWitness,
    zero_sum_constant,
)
from zero_sum_sequences.orbits import is_in_same_orbit, orbit, orbit_witness

from groups import cyclic_group, product_of_cyclic_groups


@pytest.fixture
def c3_sequences():
    return AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)


def test_map_terms_preserves_multiplicities_and_canonical_order(c3_sequences):
    sequence = c3_sequences(1, 1, 2)

    assert sequence.map_terms(lambda term: -term) == c3_sequences(1, 2, 2)
    assert sequence.map_terms(lambda term: 0) == c3_sequences(0, 0, 0)
    assert sequence.map_terms(lambda term: -term).parent() is c3_sequences


def test_map_terms_can_change_sequence_space(c3_sequences):
    target = AdditiveSequenceSpace(product_of_cyclic_groups(2), davenport_bound=2)

    image = c3_sequences(1, 2).map_terms(
        lambda term: (term % 2,),
        target_space=target,
    )

    assert image == target((1,), (0,))
    assert image.parent() is target


def test_map_terms_validates_mapping_and_target_space(c3_sequences):
    with pytest.raises(TypeError, match="mapping must be callable"):
        c3_sequences(1).map_terms(1)
    with pytest.raises(TypeError, match="target_space"):
        c3_sequences(1).map_terms(lambda term: term, target_space=object())


def test_callable_generators_produce_deterministic_orbit_and_witnesses(
    c3_sequences,
):
    # Negation generates Aut(C_3).  The sequence is unordered, so map_terms
    # also verifies that multiplicities survive the induced action.
    negation = lambda term: -term
    action = AutomorphismAction([negation])
    sequence = c3_sequences(1, 1, 2)
    image = c3_sequences(1, 2, 2)

    assert orbit(sequence, group=action) == (sequence, image)
    assert sequence.orbit(group=action) == (sequence, image)
    assert is_in_same_orbit(sequence, image, group=action)
    assert sequence.is_in_same_orbit(image, group=action)
    witness = orbit_witness(sequence, image, group=action)
    assert witness == OrbitWitness((0,))
    assert action.apply_word(sequence, witness) == image
    assert sequence.orbit_witness(image, group=action) == witness


def test_bound_witness_shows_images_of_additive_generators(capsys):
    group = FiniteAdditiveGroup.cyclic_product(2, 4)
    space = AdditiveSequenceSpace(group, davenport_bound=5)
    source = space((0, 1)) ** 3 * space((1, 0), (1, 1))
    target = space((0, 1), (1, 0)) * space((1, 1)) ** 3

    witness = source.orbit_witness(target)

    assert repr(witness) == "OrbitWitness(generator_indices=(1,))"
    assert witness.additive_generator_images() == (
        ((1, 0), (1, 0)),
        ((0, 1), (1, 1)),
    )
    witness.show()
    assert capsys.readouterr().out == "(1, 0) ↦ (1, 0)\n(0, 1) ↦ (1, 1)\n"


def test_unbound_witness_cannot_show_additive_generator_images():
    with pytest.raises(ValueError, match="returned by orbit_witness"):
        OrbitWitness((0,)).additive_generator_images()


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
    sequence = space(0, 1)

    assert orbit(sequence, group=action) == (
        space(0, 1),
        space(0, 2),
        space(1, 2),
    )


def test_empty_action_has_singleton_orbit(c3_sequences):
    sequence = c3_sequences(1, 2)

    assert orbit(sequence, group=AutomorphismAction(())) == (sequence,)
    assert sequence.orbit_witness(sequence, group=()) == OrbitWitness()
    assert not is_in_same_orbit(
        sequence,
        c3_sequences(1),
        group=AutomorphismAction(()),
    )


def test_equal_sequences_need_no_discoverable_action(c3_sequences):
    sequence = c3_sequences(1)

    assert sequence.is_in_same_orbit(sequence)
    assert sequence.orbit_witness(sequence) == OrbitWitness()
    assert c3_sequences._automorphism_action is None


def test_orbit_witness_is_none_for_different_orbits(c3_sequences):
    action = AutomorphismAction([lambda term: -term])

    assert orbit_witness(
        c3_sequences(1),
        c3_sequences(1, 1),
        group=action,
    ) is None


def test_orbit_reports_when_action_cannot_be_discovered(c3_sequences):
    with pytest.raises(AutomorphismActionUnavailable, match="pass group="):
        c3_sequences(1).orbit()


def test_parent_can_supply_generators_lazily():
    group = FiniteAdditiveGroup(
        range(3),
        zero=0,
        add=lambda left, right: (left + right) % 3,
        coerce=lambda value: int(value) % 3,
        automorphism_generators=(lambda term: -term,),
    )
    space = AdditiveSequenceSpace(group, davenport_bound=3)

    sequence = space(1, 1, 2)
    assert sequence.orbit() == (sequence, space(1, 2, 2))
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
    sequence = space(1, 1, 2)

    assert sequence.orbit() == (sequence, space(1, 2, 2))


def test_orbits_reject_different_sequence_spaces(c3_sequences):
    other_space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    with pytest.raises(TypeError, match="different spaces"):
        is_in_same_orbit(
            c3_sequences(1),
            other_space(1),
            group=[lambda term: -term],
        )


class _RecordingModule:
    """Stand-in for ``sage.all`` that records every attribute lookup."""

    def __init__(self):
        self.lookups = []

    def __getattr__(self, name):
        self.lookups.append(name)
        raise AttributeError(name)


def test_pure_python_parents_never_import_sage(monkeypatch):
    fake = _RecordingModule()
    monkeypatch.setitem(sys.modules, "sage", fake)
    monkeypatch.setitem(sys.modules, "sage.all", fake)
    space = AdditiveSequenceSpace(cyclic_group(3))

    with pytest.raises(AutomorphismActionUnavailable):
        space(1, 2).orbit()
    assert zero_sum_constant(space, {3}).value == 5
    assert fake.lookups == []


class _CountingGroup(FiniteAdditiveGroup):
    """A group without automorphism data that counts provider queries."""

    queries = 0

    def automorphism_generators(self):
        type(self).queries += 1
        raise NotImplementedError


def test_unavailable_symmetry_is_resolved_once_per_space():
    group = _CountingGroup(range(3), zero=0, add=lambda a, b: (a + b) % 3)
    space = AdditiveSequenceSpace(group)

    for _ in range(3):
        with pytest.raises(AutomorphismGroupUnavailable):
            space.automorphism_group()
        with pytest.raises(AutomorphismActionUnavailable):
            space(1, 2).orbit()
        assert zero_sum_constant(space, {3}).value == 5
    assert _CountingGroup.queries == 1


def test_orbit_queries_accept_a_materialized_group():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(2, 4))
    group = space.automorphism_group()
    source = space((0, 1), (0, 3))
    target = space((1, 1), (1, 3))

    assert source.orbit(group=group) == source.orbit() == tuple(sorted(group.orbit(source)))
    witness = source.orbit_witness(target, group=group)
    assert group.apply(group.materialize(witness.generator_indices), source) == target
    assert source.is_in_same_orbit(target, group=group)
    with pytest.raises(TypeError, match="group must be"):
        space.enumerate_atom_catalogue().orbits(group="not a group")
