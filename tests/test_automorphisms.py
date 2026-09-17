"""Materialized automorphism groups, stabilizers and canonical forms.

The checks are deliberately spread over several small coordinate groups
(cyclic, mixed, elementary abelian).  Where the group is small enough, the
automorphism group is cross-checked against a brute-force enumeration of all
bijections that are homomorphisms; otherwise against the known order.
"""

from __future__ import annotations

import itertools
import random
from math import gcd

import pytest

from groups import cyclic_group
from zero_sum_sequences import (
    AdditiveSequence,
    AdditiveSequenceSpace,
    Automorphism,
    AutomorphismGroup,
    AutomorphismGroupUnavailable,
    FiniteAdditiveGroup,
)


def brute_force_automorphisms(group):
    """Every permutation of the elements that is a homomorphism."""

    elements = tuple(group)
    found = set()
    for images in itertools.permutations(elements):
        phi = dict(zip(elements, images))
        if phi[group.zero()] != group.zero():
            continue
        if all(phi[group.add(a, b)] == group.add(phi[a], phi[b]) for a in elements for b in elements):
            found.add(tuple(elements.index(phi[e]) for e in elements))
    return found


def is_homomorphism(group, auto) -> bool:
    elements = tuple(group)
    return all(
        auto(group.add(a, b)) == group.add(auto(a), auto(b))
        for a in elements
        for b in elements
    )


def space_for(*moduli, bound):
    return AdditiveSequenceSpace(
        FiniteAdditiveGroup.cyclic_product(*moduli), davenport_bound=bound
    )


def random_sequence(space, length, rng):
    elements = [e for e in space.base_parent if e != space.base_parent.zero()]
    return space(rng.choice(elements) for _ in range(length))


# the group itself


@pytest.mark.parametrize("moduli", [(2, 4), (2, 2), (6,), (8,), (7,), (2, 2, 2)])
def test_closure_equals_brute_force_for_small_groups(moduli):
    parent = FiniteAdditiveGroup.cyclic_product(*moduli)
    group = parent.automorphism_group()
    assert {auto.table for auto in group} == brute_force_automorphisms(parent)


@pytest.mark.parametrize(
    "moduli, order",
    [
        ((42,), 12),  # (Z/42)^x
        ((5,), 4),
        ((2, 6), 12),  # S_3 x C_2
        ((4, 4), 96),  # GL(2, Z/4)
        ((3, 9), 108),
        ((2, 2, 2, 2), 20160),  # GL(4, 2)
        ((3, 3, 3), 11232),  # GL(3, 3)
    ],
)
def test_known_orders(moduli, order):
    assert len(FiniteAdditiveGroup.cyclic_product(*moduli).automorphism_group()) == order


def test_cyclic_group_automorphisms_are_multiplications_by_units():
    parent = FiniteAdditiveGroup.cyclic_product(42)
    group = parent.automorphism_group()
    units = {u for u in range(1, 42) if gcd(u, 42) == 1}
    multipliers = set()
    for auto in group:
        (u,) = auto((1,))
        multipliers.add(u)
        assert all(auto((x,)) == ((u * x) % 42,) for x in range(42))
    assert multipliers == units


@pytest.mark.parametrize("moduli", [(2, 4), (3, 9), (2, 6), (4, 4)])
def test_group_axioms_and_homomorphisms(moduli):
    parent = FiniteAdditiveGroup.cyclic_product(*moduli)
    group = parent.automorphism_group()
    assert group.identity.is_identity()
    assert group.elements[0] is group.identity
    assert len(set(group.elements)) == len(group)
    for auto in group:
        assert is_homomorphism(parent, auto)
        assert group.compose(auto, group.inverse(auto)).is_identity()
        assert group.inverse(auto) in group
    rng = random.Random(1)
    for _ in range(50):
        a, b = rng.choice(group.elements), rng.choice(group.elements)
        c = group.compose(a, b)
        assert c in group
        assert all(c(x) == a(b(x)) for x in parent)


def test_generators_generate_and_are_elements():
    group = FiniteAdditiveGroup.cyclic_product(2, 4).automorphism_group()
    assert all(g in group for g in group.generators)
    assert group.identity not in group.generators
    # closure of the stored generators is the whole group
    seen = {group.identity}
    pending = [group.identity]
    while pending:
        current = pending.pop()
        for g in group.generators:
            image = group.compose(g, current)
            if image not in seen:
                seen.add(image)
                pending.append(image)
    assert seen == set(group.elements)


def test_generator_images_determine_elements_and_give_the_order():
    group = FiniteAdditiveGroup.cyclic_product(2, 4).automorphism_group()
    images = [auto.generator_images() for auto in group]
    assert len(set(images)) == len(group)
    assert images[1:] == sorted(images[1:])  # canonical order after the identity
    identity_images = group.identity.generator_images()
    assert identity_images == ((1, 0), (0, 1))


# actions on sequences


@pytest.mark.parametrize(
    "moduli, bound",
    [((42,), 8), ((2, 4), 5), ((2, 2, 2, 2), 5), ((4, 4), 7), ((3, 9), 11)],
)
def test_orbit_stabilizer(moduli, bound):
    space = space_for(*moduli, bound=bound)
    group = space.automorphism_group()
    rng = random.Random(str(moduli))
    for _ in range(6):
        x = random_sequence(space, rng.randint(1, bound), rng)
        orbit = group.orbit(x)
        stabilizer = group.stabilizer(x)
        assert len(orbit) * len(stabilizer) == len(group)
        assert x in orbit and orbit == tuple(sorted(orbit))
        assert all(group.apply(s, x) == x for s in stabilizer)
        assert stabilizer.identity in stabilizer


@pytest.mark.parametrize("moduli, bound", [((2, 4), 5), ((6,), 6), ((2, 6), 7)])
def test_canonical_form_is_an_orbit_invariant(moduli, bound):
    space = space_for(*moduli, bound=bound)
    group = space.automorphism_group()
    rng = random.Random(0)
    for _ in range(8):
        x = random_sequence(space, rng.randint(1, bound), rng)
        form = group.canonical_form(x)
        assert form.image == min(group.orbit(x))
        assert group.apply(form.automorphism, x) == form.image
        assert form.permutation is None
        for auto in group:
            assert group.canonical(group.apply(auto, x)) == form.image
        # the smallest witness really is the smallest transporter
        transporters = [a for a in group if group.apply(a, x) == form.image]
        assert form.automorphism == transporters[0]
        assert group.canonical_form(x, witness="word").image == form.image
        word = group.canonical_form(x, witness="word").automorphism
        assert group.apply(word, x) == form.image


def test_transporter():
    space = space_for(2, 4, bound=5)
    group = space.automorphism_group()
    x = space([(1, 1), (1, 1), (0, 2)])
    y = space([(1, 3), (1, 3), (0, 2)])
    auto = group.transporter(x, y)
    assert auto is not None and group.apply(auto, x) == y
    assert group.transporter(x, space([(0, 1), (0, 3)])) is None
    assert group.transporter(x, x) is group.identity


def test_orbit_representatives_matches_atom_classification():
    # atoms of C_2 x C_4 up to automorphism: known counts per length
    space = space_for(2, 4, bound=5)
    group = space.automorphism_group()
    atoms = tuple(space.enumerate_atom_catalogue())
    buckets = group.orbit_representatives(atoms)
    assert sum(len(members) for members in buckets.values()) == len(atoms)
    assert all(rep in members for rep, members in buckets.items())
    # independent count: number of distinct orbit minima over the whole group
    minima = {min(group.apply(a, atom) for a in group) for atom in atoms}
    assert set(buckets) == minima
    for rep, members in buckets.items():
        assert len(members) * len(group.stabilizer(rep)) == len(group)


@pytest.mark.parametrize(
    "moduli, bound, orbits",
    [((2, 2, 2), 4, 3), ((3, 3), 5, 5), ((2, 2, 2, 2), 5, 4)],
)
def test_elementary_abelian_atom_orbits(moduli, bound, orbits):
    space = space_for(*moduli, bound=bound)
    atoms = tuple(space.enumerate_atom_catalogue())
    assert len(space.automorphism_group().orbit_representatives(atoms)) == orbits


# tuples of sequences


def test_tuple_actions_multiset_and_ordered():
    space = space_for(2, 4, bound=5)
    group = space.automorphism_group()
    x = space([(0, 1), (0, 1), (0, 2)])
    y = space([(1, 0), (1, 0)])
    z = space([(1, 1), (1, 3)])
    pair = (x, y)
    assert group.orbit(pair) == group.orbit((y, x))
    assert group.stabilizer(pair).elements == group.stabilizer((y, x)).elements
    assert group.canonical(pair) == group.canonical((y, x))
    assert group.canonical((x, y), ordered=True) != group.canonical((y, x), ordered=True) or x == y
    form = group.canonical_form((y, z, x))
    image = form.image
    raw = group.apply(form.automorphism, (y, z, x), ordered=True)
    assert image == tuple(sorted(raw))
    assert all(image[form.permutation[i]] == raw[i] for i in range(3))
    assert sorted(form.permutation) == [0, 1, 2]
    # orbit-stabilizer also for tuples
    assert len(group.orbit(pair)) * len(group.stabilizer(pair)) == len(group)


def test_tuple_objects_must_share_a_space():
    a = space_for(2, 4, bound=5)
    b = space_for(2, 4, bound=6)
    with pytest.raises(TypeError):
        a.automorphism_group().orbit((a([(0, 1)]), b([(0, 1)])))
    with pytest.raises(TypeError):
        a.automorphism_group().orbit("no")


# delegation


def test_parent_supplied_group_is_used():
    modulus = 42
    units = tuple(u for u in range(1, modulus) if gcd(u, modulus) == 1)

    def build():
        return AutomorphismGroup(
            units,
            apply_term=lambda u, term: ((u * term[0]) % modulus,),
            compose=lambda a, b: (a * b) % modulus,
            inverse=lambda a: pow(a, -1, modulus),
            identity=1,
            generators=(5, 11),
        )

    parent = FiniteAdditiveGroup(
        ((x,) for x in range(modulus)),
        zero=(0,),
        add=lambda a, b: ((a[0] + b[0]) % modulus,),
        automorphism_group=build,
    )
    space = AdditiveSequenceSpace(parent, davenport_bound=modulus)
    group = space.automorphism_group()
    assert group.elements[0] == 1 and len(group) == 12
    assert group.elements[1:] == tuple(sorted(units)[1:])
    reference = space_for(42, bound=42).automorphism_group()
    rng = random.Random(3)
    for _ in range(10):
        x = space((rng.choice(units),) for _ in range(rng.randint(1, 6)))
        mirrored = space_for(42, bound=42)(tuple(x))
        assert [tuple(s) for s in group.orbit(x)] == [tuple(s) for s in reference.orbit(mirrored)]
        assert len(group.stabilizer(x)) == len(reference.stabilizer(mirrored))
    assert space.automorphism_group() is group  # cached


def test_supplied_generators_must_generate():
    units = tuple(u for u in range(1, 42) if gcd(u, 42) == 1)
    with pytest.raises(ValueError):
        AutomorphismGroup(
            units,
            apply_term=lambda u, term: ((u * term[0]) % 42,),
            compose=lambda a, b: (a * b) % 42,
            identity=1,
            generators=(5, 41),  # 41 = 5^3, so these generate a subgroup of order 6
        )


def test_provider_must_return_a_group():
    parent = FiniteAdditiveGroup(
        range(3), zero=0, add=lambda a, b: (a + b) % 3, automorphism_group=lambda: "no"
    )
    with pytest.raises(TypeError):
        parent.automorphism_group()
    with pytest.raises(TypeError):
        FiniteAdditiveGroup(range(3), zero=0, add=lambda a, b: (a + b) % 3, automorphism_group=1)


def test_unavailable_without_generators_or_provider():
    space = AdditiveSequenceSpace(cyclic_group(7), davenport_bound=7)
    with pytest.raises(AutomorphismGroupUnavailable):
        space.automorphism_group()
    with pytest.raises(NotImplementedError):
        cyclic_group(7).automorphism_group()


def test_closure_rejects_non_bijections():
    parent = FiniteAdditiveGroup.cyclic_product(4)
    with pytest.raises(ValueError):
        AutomorphismGroup.closure(parent, [lambda x: (0,)])


def test_subgroup_and_element_api():
    group = FiniteAdditiveGroup.cyclic_product(2, 2).automorphism_group()
    assert len(group) == 6
    sub = group.subgroup([group.identity])
    assert len(sub) == 1 and sub.identity is group.identity
    auto = group.elements[1]
    assert isinstance(auto, Automorphism)
    assert auto.inverse().compose(auto).is_identity()
    assert auto == Automorphism(auto.table, auto._data) and hash(auto) == hash(auto.table)
    assert "generator_images" in repr(auto)


# the order on sequences used for canonical forms


def test_sequence_order_is_length_then_terms():
    space = space_for(6, bound=6)
    a = space([(5,)])
    b = space([(1,), (1,)])
    c = space([(1,), (2,)])
    assert a < b < c and c > b > a and a <= a and c >= c
    assert sorted([c, a, b]) == [a, b, c]
    with pytest.raises(TypeError):
        _ = a < space_for(6, bound=7)([(1,)])
    assert not isinstance(AdditiveSequence.__lt__(a, 3), bool)


# subgroups, words and transporters


@pytest.mark.parametrize("moduli, bound", [((2, 4), 5), ((4, 4), 7), ((3, 3), 5), ((42,), 8)])
def test_subgroup_generators_generate_and_are_few(moduli, bound):
    space = space_for(*moduli, bound=bound)
    group = space.automorphism_group()
    rng = random.Random(str(moduli))
    x = random_sequence(space, 3, rng)
    stabilizer = group.stabilizer(x)
    assert len(stabilizer.generators) <= max(1, len(stabilizer).bit_length())
    closure = {stabilizer.identity}
    pending = [stabilizer.identity]
    while pending:
        current = pending.pop()
        for g in stabilizer.generators:
            image = stabilizer.compose(g, current)
            if image not in closure:
                closure.add(image)
                pending.append(image)
    assert closure == set(stabilizer.elements)
    with pytest.raises(ValueError):
        group.subgroup(group.elements[1:3])  # not closed


def test_orbit_words_and_materialize():
    space = space_for(2, 4, bound=5)
    group = space.automorphism_group()
    x = space([(1, 1), (1, 1), (0, 2)])
    words = group.orbit_words(x)
    assert set(words) == set(group.orbit(x)) and words[x] == ()
    for image, word in words.items():
        assert group.apply(group.materialize(word), x) == image
        assert all(0 <= i < len(group.generators) for i in word)


@pytest.mark.parametrize("moduli, bound", [((2, 4), 5), ((6,), 6), ((2, 2, 2), 4)])
def test_catalogue_representatives_and_transporters(moduli, bound):
    space = space_for(*moduli, bound=bound)
    catalogue = space.enumerate_atom_catalogue()
    group = space.automorphism_group()
    for atom in catalogue:
        representative = catalogue.representative(atom)
        assert representative == group.canonical(atom)
        transporter = catalogue.transporter(atom)
        assert group.apply(transporter, atom) == representative
    assert {catalogue.representative(a) for a in catalogue} == {
        orbit.representative for orbit in catalogue.orbits()
    }
    with pytest.raises(ValueError):
        catalogue.transporter(catalogue[0] + catalogue[1])
