"""Catalogue symmetry and input validation must not change cached answers."""

import json

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AnchoredCanonicalizer,
    AtomCatalogue,
    AutomorphismGroup,
    FactorizationCache,
    FiniteAdditiveGroup,
)


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
