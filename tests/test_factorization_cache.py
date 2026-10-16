"""The memoizing, optionally group-aware factorization cache."""

from __future__ import annotations

import itertools
import json
import random

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AnchoredCanonicalizer,
    AtomCatalogue,
    AutomorphismGroup,
    FactorizationCache,
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
    cache = space.factorization_cache(catalogue)
    rng = random.Random(1)
    for left, right in pairs(catalogue, rng, 25):
        sequence = left + right
        solver = FactorizationSolver(sequence, atom_catalogue=catalogue)
        assert cache.length_set(sequence) == frozenset(solver.length_set())
        assert cache.minimum(sequence) == solver.minimum_factorization_length()
        assert cache.maximum(sequence) == solver.maximum_factorization_length()
        for length in range(1, 6):
            assert cache.has_length(sequence, length) == solver.has_factorization_of_length(length)
            witness = cache.witness(sequence, length)
            assert (witness is None) == (length not in cache.length_set(sequence))
            if witness is not None:
                assert product(space, witness) == sequence and len(witness) == length
        assert set(cache.witnesses(sequence)) == cache.length_set(sequence)
        assert sum(1 for _ in cache.factorizations(sequence)) == sum(
            1 for _ in solver.factorizations()
        )


def test_caching_and_statistics():
    space = space_for(2, 4, bound=5)
    catalogue = space.enumerate_atom_catalogue()
    cache = space.factorization_cache(catalogue, maxsize=2)
    a, b, c = catalogue[-1], catalogue[-2], catalogue[-3]
    x, y, z = a + b, a + c, b + c
    assert cache.length_set(x) == cache.length_set(x)
    stats = cache.statistics()
    assert stats["hits"] == 1 and stats["misses"] == 1 and stats["solver_builds"] == 1
    cache.length_set(y)
    cache.length_set(z)
    stats = cache.statistics()
    assert stats["cached_length_sets"] == 3 and stats["cached_solvers"] == 2  # LRU bound
    # a length set is answered from the cache even after its solver was evicted
    assert cache.minimum(x) == min(cache.length_set(x))
    assert cache.statistics()["solver_builds"] == 3
    cache.clear()
    assert cache.statistics()["cached_length_sets"] == 0
    with pytest.raises(ValueError):
        space.factorization_cache(catalogue, maxsize=0)
    with pytest.raises(TypeError):
        cache.length_set(space_for(2, 4, bound=6)([(0, 1), (0, 3)]))
    with pytest.raises(TypeError):
        space.factorization_cache(space_for(2, 4, bound=6).enumerate_atom_catalogue())


def test_group_aware_cache_answers_orbits_from_one_solve(setup):
    space, catalogue = setup
    group = space.automorphism_group()
    plain = space.factorization_cache(catalogue)
    aware = space.factorization_cache(catalogue, group=group)
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
    cache = space.factorization_cache(catalogue, group=restricted)
    sequence = catalogue[-1] + catalogue[-3]
    assert cache.length_set(sequence) == space.factorization_cache(catalogue).length_set(sequence)
    with pytest.raises(TypeError):  # transporting a witness needs compose/inverse
        cache.witness(sequence, min(cache.length_set(sequence)))


def test_product_helper(setup):
    space, catalogue = setup
    cache = space.factorization_cache(catalogue)
    atoms = (catalogue[0], catalogue[-1], catalogue[-1])
    assert cache.product(atoms) == product(space, atoms)
    assert cache.product(()) == space(())
    assert cache.length_set(cache.product(atoms)) == cache.length_set(product(space, atoms))


def test_persistence_round_trip(setup, tmp_path):
    space, catalogue = setup
    group = space.automorphism_group()
    cache = space.factorization_cache(catalogue, group=group)
    rng = random.Random(3)
    queried = [left + right for left, right in pairs(catalogue, rng, 15)]
    expected = {sequence: cache.length_set(sequence) for sequence in queried}
    path = tmp_path / "lengths.jsonl"
    count = cache.to_jsonl(path)
    assert count == cache.statistics()["cached_length_sets"]

    loaded = FactorizationCache.from_jsonl(path, catalogue, group=group)
    for sequence, lengths in expected.items():
        assert loaded.length_set(sequence) == lengths
    assert loaded.statistics()["misses"] == 0 and loaded.statistics()["solver_builds"] == 0
    records = list(loaded.records())
    assert records == list(cache.records())
    assert [r["sequence"] for r in records] == [
        r["sequence"] for r in sorted(records, key=lambda r: (len(r["sequence"]), r["sequence"]))
    ]
    # the file is bound to the catalogue and to the group
    with pytest.raises(ValueError):
        FactorizationCache.from_jsonl(path, catalogue)
    with pytest.raises(ValueError):
        FactorizationCache.from_jsonl(path, catalogue.restrict(max_length=2), group=group)


def test_preload_accepts_external_tables():
    space = space_for(6, bound=6)
    catalogue = space.enumerate_atom_catalogue()
    cache = space.factorization_cache(catalogue)
    sequence = catalogue[-1] + catalogue[-2]
    truth = FactorizationSolver(sequence, atom_catalogue=catalogue).length_set()
    assert cache.preload([(sequence, sorted(truth))]) == 1
    assert cache.length_set(sequence) == frozenset(truth)
    assert cache.statistics()["solver_builds"] == 0
    with pytest.raises(TypeError):
        cache.preload([("not a sequence", [1])])


def test_exhaustive_pair_table_is_consistent_with_orbit_classification():
    # every pair of atoms of C_2 x C_4: the group-aware cache needs one solve
    # per orbit of unordered pairs and agrees with the plain cache
    space = space_for(2, 4, bound=5)
    catalogue = space.enumerate_atom_catalogue()
    group = space.automorphism_group()
    plain = space.factorization_cache(catalogue, maxsize=None)
    aware = space.factorization_cache(catalogue, group=group, maxsize=None)
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


def test_cache_with_anchored_canonicalizer(setup):
    from zero_sum_sequences import AnchoredCanonicalizer

    space, catalogue = setup
    group = space.automorphism_group()
    plain = space.factorization_cache(catalogue)
    cache = space.factorization_cache(
        catalogue, group=group, canonicalize=AnchoredCanonicalizer(catalogue, group)
    )
    rng = random.Random(5)
    for left, right in pairs(catalogue, rng, 10):
        sequence = left + right
        lengths = plain.length_set(sequence)
        for element in group.elements[:: max(1, len(group) // 5)]:
            image = group.apply(element, sequence)
            assert cache.length_set(image) == lengths
            witness = cache.witness(image, min(lengths))
            assert product(space, witness) == image
    # a canonicalizer returning only the image works too (transporter is searched)
    plain_key = space.factorization_cache(catalogue, group=group, canonicalize=lambda s: group.canonical(s))
    sequence = catalogue[-1] + catalogue[-2]
    witness = plain_key.witness(sequence, min(plain.length_set(sequence)))
    assert product(space, witness) == sequence
    with pytest.raises(ValueError):
        space.factorization_cache(catalogue, canonicalize=lambda s: s)
    with pytest.raises(TypeError):
        space.factorization_cache(catalogue, group=group, canonicalize="no")


@pytest.fixture
def problem():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    atom = space([(2,)] * 3)
    return space, atom, AtomCatalogue(space, [atom])


def test_group_aware_cache_rejects_a_nonclosed_catalogue(problem):
    space, atom, catalogue = problem
    group = space.automorphism_group()
    assert atom.length_set(atom_catalogue=catalogue) == {1}
    for factory in (
        lambda: space.factorization_cache(catalogue, group=group),
        lambda: AnchoredCanonicalizer(catalogue, group),
    ):
        with pytest.raises(ValueError, match="closed"):
            factory()
    trivial = group.subgroup([group.identity])
    cache = space.factorization_cache(catalogue, group=trivial)
    assert cache.length_set(atom) == {1}
    assert cache.witness(atom, 1) == (atom,)


def test_closure_validation_is_generator_based_and_reused(monkeypatch):
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(2, 4))
    catalogue = space.enumerate_atom_catalogue()
    group = space.automorphism_group()
    original = AtomCatalogue.permutation
    checks = []

    def counted(self, element, *, group=None):
        checks.append(element)
        return original(self, element, group=group)

    monkeypatch.setattr(AtomCatalogue, "permutation", counted)
    anchored = AnchoredCanonicalizer(catalogue, group)
    cache = space.factorization_cache(catalogue, group=group, canonicalize=anchored)
    for atom in catalogue:
        assert cache.length_set(atom) == {1}
    assert checks == list(group.generators)


def test_closed_catalogue_works_without_composition_callbacks():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    catalogue = space.enumerate_atom_catalogue()
    full = space.automorphism_group()
    group = AutomorphismGroup(full.elements, identity=full.identity, apply_term=full.apply_term)
    cache = space.factorization_cache(catalogue, group=group)
    assert cache.length_set(catalogue[-1]) == {1}


def test_plain_cache_never_discovers_symmetry(problem, monkeypatch):
    space, atom, catalogue = problem

    def forbidden(*args):
        raise AssertionError("plain cache must not discover symmetry")

    monkeypatch.setattr(AdditiveSequenceSpace, "automorphism_group", forbidden)
    assert space.factorization_cache(catalogue).length_set(atom) == {1}


@pytest.mark.parametrize("state", ["cold", "warm", "preloaded"])
@pytest.mark.parametrize("bad", [True, False, 1.0, 1.5, -1, "1", None])
def test_length_validation_does_not_depend_on_cache_state(problem, state, bad):
    space, atom, catalogue = problem
    cache = space.factorization_cache(catalogue)
    if state == "warm":
        cache.length_set(atom)
    elif state == "preloaded":
        cache.preload([(atom, [1])])
    with pytest.raises(ValueError, match="factor count"):
        cache.has_length(atom, bad)


@pytest.mark.parametrize("bad", [True, False, 0, -1, 1.0, 1.5, float("inf"), "1"])
def test_cache_maxsize_requires_an_exact_positive_integer(problem, bad):
    space, _, catalogue = problem
    with pytest.raises(ValueError, match="maxsize"):
        space.factorization_cache(catalogue, maxsize=bad)


@pytest.mark.parametrize("bad", [True, False, 1.5, -1, "1"])
def test_preload_rejects_invalid_lengths_without_inserting_a_record(problem, bad):
    space, atom, catalogue = problem
    cache = space.factorization_cache(catalogue)
    with pytest.raises(ValueError, match="factor count"):
        cache.preload([(atom, [1, bad])])
    assert cache.statistics()["cached_length_sets"] == 0


def test_indexable_counts_and_cache_size_are_supported(problem):
    class Indexable:
        def __index__(self):
            return 1

    space, atom, catalogue = problem
    cache = space.factorization_cache(catalogue, maxsize=Indexable())
    cache.preload([(atom, [Indexable()])])
    assert cache.has_length(atom, Indexable())


@pytest.mark.parametrize("query", ["has_length", "witness"])
def test_factor_counts_are_validated_before_solver_lookup(problem, query, monkeypatch):
    space, atom, catalogue = problem
    cache = space.factorization_cache(catalogue)

    def forbidden(*args):
        raise AssertionError("invalid requests must not look up or build solvers")

    monkeypatch.setattr(FactorizationCache, "_solver_for_key", forbidden)
    with pytest.raises(ValueError, match="factor count"):
        getattr(cache, query)(atom, True)


def test_cache_jsonl_rejects_nonintegral_lengths(problem, tmp_path):
    space, atom, catalogue = problem
    cache = space.factorization_cache(catalogue)
    cache.length_set(atom)
    path = tmp_path / "cache.jsonl"
    cache.to_jsonl(path)
    metadata, record = [json.loads(line) for line in path.read_text().splitlines()]
    record["length_set"] = [1.5]
    path.write_text(json.dumps(metadata) + "\n" + json.dumps(record) + "\n")
    with pytest.raises(ValueError, match="factor count"):
        FactorizationCache.from_jsonl(path, catalogue)
