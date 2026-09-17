"""The memoizing, optionally group-aware factorization oracle."""

from __future__ import annotations

import itertools
import random

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AtomCatalogue,
    FactorizationOracle,
    FactorizationSolver,
    FiniteAdditiveGroup,
)


def space_for(*moduli, bound):
    return AdditiveSequenceSpace(
        FiniteAdditiveGroup.cyclic_product(*moduli), davenport_bound=bound
    )


def product(space, atoms):
    result = space(())
    for atom in atoms:
        result = result + atom
    return result


@pytest.fixture(params=[((2, 4), 5), ((6,), 6), ((3, 3), 5), ((2, 2, 2), 4)])
def setup(request):
    moduli, bound = request.param
    space = space_for(*moduli, bound=bound)
    catalogue = space.enumerate_atom_catalogue()
    return space, catalogue


def pairs(catalogue, rng, count):
    atoms = catalogue.atoms
    return [
        (rng.choice(atoms), rng.choice(atoms)) for _ in range(count)
    ]


def test_answers_match_the_solver(setup):
    space, catalogue = setup
    oracle = space.oracle(catalogue)
    rng = random.Random(1)
    for left, right in pairs(catalogue, rng, 25):
        sequence = left + right
        solver = FactorizationSolver(sequence, atom_catalogue=catalogue)
        assert oracle.length_set(sequence) == frozenset(solver.length_set())
        assert oracle.minimum(sequence) == solver.minimum_factorization_length()
        assert oracle.maximum(sequence) == solver.maximum_factorization_length()
        for length in range(1, 6):
            assert oracle.has_length(sequence, length) == solver.has_factorization_of_length(length)
            witness = oracle.witness(sequence, length)
            assert (witness is None) == (length not in oracle.length_set(sequence))
            if witness is not None:
                assert product(space, witness) == sequence and len(witness) == length
        assert set(oracle.witnesses(sequence)) == oracle.length_set(sequence)
        assert sum(1 for _ in oracle.factorizations(sequence)) == sum(
            1 for _ in solver.factorizations()
        )


def test_caching_and_statistics():
    space = space_for(2, 4, bound=5)
    catalogue = space.enumerate_atom_catalogue()
    oracle = space.oracle(catalogue, maxsize=2)
    a, b, c = catalogue[-1], catalogue[-2], catalogue[-3]
    x, y, z = a + b, a + c, b + c
    assert oracle.length_set(x) == oracle.length_set(x)
    stats = oracle.statistics()
    assert stats["hits"] == 1 and stats["misses"] == 1 and stats["solver_builds"] == 1
    oracle.length_set(y)
    oracle.length_set(z)
    stats = oracle.statistics()
    assert stats["cached_length_sets"] == 3 and stats["cached_solvers"] == 2  # LRU bound
    # a length set is answered from the cache even after its solver was evicted
    assert oracle.minimum(x) == min(oracle.length_set(x))
    assert oracle.statistics()["solver_builds"] == 3
    oracle.clear()
    assert oracle.statistics()["cached_length_sets"] == 0
    with pytest.raises(ValueError):
        space.oracle(catalogue, maxsize=0)
    with pytest.raises(TypeError):
        oracle.length_set(space_for(2, 4, bound=6)([(0, 1), (0, 3)]))
    with pytest.raises(TypeError):
        space.oracle(space_for(2, 4, bound=6).enumerate_atom_catalogue())


def test_group_aware_oracle_answers_orbits_from_one_solve(setup):
    space, catalogue = setup
    group = space.automorphism_group()
    plain = space.oracle(catalogue)
    aware = space.oracle(catalogue, group=group)
    rng = random.Random(2)
    for left, right in pairs(catalogue, rng, 10):
        sequence = left + right
        lengths = plain.length_set(sequence)
        assert aware.length_set(sequence) == lengths
        builds = aware.statistics()["solver_builds"]
        for element in group.elements[:: max(1, len(group) // 6)]:
            image = group.apply(element, sequence)
            assert aware.length_set(image) == lengths
            assert aware.minimum(image) == min(lengths)
            assert aware.maximum(image) == max(lengths)
            for length in lengths:
                witness = aware.witness(image, length)
                assert product(space, witness) == image and len(witness) == length
                assert all(atom in catalogue for atom in witness)
            assert all(product(space, f) == image for f in aware.factorizations(image))
        assert aware.statistics()["solver_builds"] == builds  # the orbit was one solve


def test_group_without_compose_answers_length_sets_only():
    space = space_for(2, 4, bound=5)
    catalogue = space.enumerate_atom_catalogue()
    full = space.automorphism_group()
    from zero_sum_sequences import AutomorphismGroup

    # the same elements, but without compose/inverse: no word witnesses
    restricted = AutomorphismGroup(
        full.elements, apply_term=full.apply_term, identity=full.identity
    )
    oracle = space.oracle(catalogue, group=restricted)
    sequence = catalogue[-1] + catalogue[-3]
    assert oracle.length_set(sequence) == space.oracle(catalogue).length_set(sequence)
    with pytest.raises(TypeError):  # transporting a witness needs compose/inverse
        oracle.witness(sequence, min(oracle.length_set(sequence)))


def test_product_helper(setup):
    space, catalogue = setup
    oracle = space.oracle(catalogue)
    atoms = (catalogue[0], catalogue[-1], catalogue[-1])
    assert oracle.product(atoms) == product(space, atoms)
    assert oracle.product(()) == space(())
    assert oracle.length_set(oracle.product(atoms)) == oracle.length_set(product(space, atoms))


def test_persistence_round_trip(setup, tmp_path):
    space, catalogue = setup
    group = space.automorphism_group()
    oracle = space.oracle(catalogue, group=group)
    rng = random.Random(3)
    queried = [left + right for left, right in pairs(catalogue, rng, 15)]
    expected = {sequence: oracle.length_set(sequence) for sequence in queried}
    path = tmp_path / "lengths.jsonl"
    count = oracle.to_jsonl(path)
    assert count == oracle.statistics()["cached_length_sets"]

    loaded = FactorizationOracle.from_jsonl(path, catalogue, group=group)
    for sequence, lengths in expected.items():
        assert loaded.length_set(sequence) == lengths
    assert loaded.statistics()["misses"] == 0 and loaded.statistics()["solver_builds"] == 0
    records = list(loaded.records())
    assert records == list(oracle.records())
    assert [r["sequence"] for r in records] == [
        r["sequence"] for r in sorted(records, key=lambda r: (len(r["sequence"]), r["sequence"]))
    ]
    # the file is bound to the catalogue and to the group
    with pytest.raises(ValueError):
        FactorizationOracle.from_jsonl(path, catalogue)
    with pytest.raises(ValueError):
        FactorizationOracle.from_jsonl(path, catalogue.restrict(max_length=2), group=group)


def test_preload_accepts_external_tables():
    space = space_for(6, bound=6)
    catalogue = space.enumerate_atom_catalogue()
    oracle = space.oracle(catalogue)
    sequence = catalogue[-1] + catalogue[-2]
    truth = FactorizationSolver(sequence, atom_catalogue=catalogue).length_set()
    assert oracle.preload([(sequence, sorted(truth))]) == 1
    assert oracle.length_set(sequence) == frozenset(truth)
    assert oracle.statistics()["solver_builds"] == 0
    with pytest.raises(TypeError):
        oracle.preload([("not a sequence", [1])])


def test_exhaustive_pair_table_is_consistent_with_orbit_classification():
    # every pair of atoms of C_2 x C_4: the group-aware oracle needs one solve
    # per orbit of unordered pairs and agrees with the plain oracle
    space = space_for(2, 4, bound=5)
    catalogue = space.enumerate_atom_catalogue()
    group = space.automorphism_group()
    plain = space.oracle(catalogue, maxsize=None)
    aware = space.oracle(catalogue, group=group, maxsize=None)
    pair_orbits = group.orbit_representatives(
        [(a, b) for a, b in itertools.combinations_with_replacement(catalogue.atoms, 2)]
    )
    for a, b in itertools.combinations_with_replacement(catalogue.atoms, 2):
        assert aware.length_set(a + b) == plain.length_set(a + b)
    # products of pairs in one orbit are in one orbit, so at most one solve each
    assert aware.statistics()["solver_builds"] <= len(pair_orbits)


def test_anchored_canonicalizer_gives_automorphic_images(setup):
    from zero_sum_sequences import AnchoredCanonicalizer

    space, catalogue = setup
    group = space.automorphism_group()
    rng = random.Random(4)
    for reduce in (True, False):
        anchored = AnchoredCanonicalizer(catalogue, group, reduce=reduce)
        for left, right in pairs(catalogue, rng, 10):
            sequence = left + right
            image, automorphism = anchored(sequence)
            assert group.apply(automorphism, sequence) == image
            assert automorphism in group
            # orbit-mates receive keys in the same orbit, often the same key
            mate = group.apply(group.elements[-1], sequence)
            mate_image, _ = anchored(mate)
            assert group.canonical(mate_image) == group.canonical(image)
        # a sequence with no catalogue divisor is its own key
        atom_free = space([left.support[0]])
        assert anchored(atom_free) == (atom_free, group.identity)


def test_oracle_with_anchored_canonicalizer(setup):
    from zero_sum_sequences import AnchoredCanonicalizer

    space, catalogue = setup
    group = space.automorphism_group()
    plain = space.oracle(catalogue)
    oracle = space.oracle(
        catalogue, group=group, canonicalize=AnchoredCanonicalizer(catalogue, group)
    )
    rng = random.Random(5)
    for left, right in pairs(catalogue, rng, 10):
        sequence = left + right
        lengths = plain.length_set(sequence)
        for element in group.elements[:: max(1, len(group) // 5)]:
            image = group.apply(element, sequence)
            assert oracle.length_set(image) == lengths
            witness = oracle.witness(image, min(lengths))
            assert product(space, witness) == image
    # a canonicalizer returning only the image works too (transporter is searched)
    plain_key = space.oracle(catalogue, group=group, canonicalize=lambda s: group.canonical(s))
    sequence = catalogue[-1] + catalogue[-2]
    witness = plain_key.witness(sequence, min(plain.length_set(sequence)))
    assert product(space, witness) == sequence
    with pytest.raises(ValueError):
        space.oracle(catalogue, canonicalize=lambda s: s)
    with pytest.raises(TypeError):
        space.oracle(catalogue, group=group, canonicalize="no")
