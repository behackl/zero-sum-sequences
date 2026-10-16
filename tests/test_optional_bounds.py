"""Local factorization and catalogue persistence without an ambient bound."""

import json

import networkx as nx
import pytest

from zero_sum_sequences import (
    AdditiveSequence,
    AdditiveSequenceSpace,
    AtomCatalogue,
    FiniteAdditiveGroup,
)


class IntegerParent:
    def __call__(self, term):
        return int(term)

    def zero(self):
        return 0

    def is_finite(self):
        return False


@pytest.fixture
def integer_problem():
    space = AdditiveSequenceSpace(IntegerParent())
    pair = space([-1, 1])
    other_pair = space([-2, 2])
    positive = space([-2, 1, 1])
    negative = space([-1, -1, 2])
    catalogue = AtomCatalogue(space, [pair, other_pair, positive, negative])
    sequence = positive + negative
    expected = {
        tuple(sorted((positive, negative))),
        tuple(sorted((pair, pair, other_pair))),
    }
    return space, sequence, catalogue, expected


@pytest.mark.parametrize("use_catalogue", [False, True])
def test_boundless_factorization_queries(integer_problem, use_catalogue):
    space, sequence, catalogue, expected = integer_problem
    solver = sequence.factorization_solver(
        atom_catalogue=catalogue if use_catalogue else None
    )

    assert set(solver.factorizations()) == expected
    assert solver.length_set() == {2, 3}
    assert solver.minimum_factorization_length() == 2
    assert solver.maximum_factorization_length() == 3
    assert set(solver.factorization_witnesses()) == {2, 3}
    for length in range(5):
        assert solver.has_factorization_of_length(length) == (length in {2, 3})
        witness = solver.factorization_witness(length)
        assert (witness is not None) == (length in {2, 3})
        if witness is not None:
            assert len(witness) == length
            assert sum(witness, space()) == sequence
    assert solver.factorization_witness(
        3, minimum_matching_factors=3, factor_predicate=lambda atom: len(atom) == 2
    ) in expected
    assert solver.factorization_witness(
        2, minimum_matching_factors=1, factor_predicate=lambda atom: len(atom) == 2
    ) is None

    graph = solver.digraph()
    assert nx.is_directed_acyclic_graph(graph)
    for source, target, data in graph.edges(data=True):
        assert target + data["atom"] == source
    if use_catalogue:
        assert space._arithmetic_context is None
    with pytest.raises(ValueError, match="supply davenport_bound explicitly"):
        _ = solver.davenport_bound


def test_boundless_factorization_cache(integer_problem, tmp_path):
    space, sequence, catalogue, expected = integer_problem
    cache = space.factorization_cache(catalogue)
    assert cache.length_set(sequence) == frozenset({2, 3})
    assert cache.minimum(sequence) == 2
    assert cache.maximum(sequence) == 3
    assert set(cache.factorizations(sequence)) == expected
    assert cache.witness(sequence, 3) in expected
    assert cache.witness(sequence, 4) is None
    assert set(cache.witnesses(sequence)) == {2, 3}
    assert space._arithmetic_context is None
    path = tmp_path / "lengths.jsonl"
    cache.to_jsonl(path)
    restored = type(cache).from_jsonl(path, catalogue)
    assert restored.length_set(sequence) == frozenset({2, 3})


def test_catalogue_backed_solving_never_resolves_a_bound(integer_problem, monkeypatch):
    space, sequence, catalogue, expected = integer_problem

    def forbidden(self):
        raise AssertionError("catalogue-backed solving must not resolve a bound")

    monkeypatch.setattr(AdditiveSequenceSpace, "davenport_bound", property(forbidden))
    assert set(sequence.factorizations(atom_catalogue=catalogue)) == expected
    assert space.factorization_cache(catalogue).length_set(sequence) == frozenset({2, 3})
    assert space._arithmetic_context is None


def test_empty_factorization_never_resolves_a_bound(monkeypatch):
    space = AdditiveSequenceSpace(IntegerParent())

    def forbidden(self):
        raise AssertionError("empty input must not resolve a bound")

    monkeypatch.setattr(AdditiveSequenceSpace, "davenport_bound", property(forbidden))
    solver = space().factorization_solver()
    assert solver.length_set() == {0}
    assert list(solver.factorizations()) == [()]
    assert solver.factorization_witness(0) == ()
    assert solver.minimum_factorization_length() == 0
    assert solver.maximum_factorization_length() == 0
    assert space._arithmetic_context is None


@pytest.mark.parametrize("bound", [None, 2, 3])
def test_large_inputs_keep_the_bound_capped_search(bound, monkeypatch):
    group = FiniteAdditiveGroup.cyclic_product(3)
    space = AdditiveSequenceSpace(group, davenport_bound=bound)
    sequence = space([(1,)] * 15_000 + [(2,)] * 15_000)
    expected_bound = 3 if bound is None else bound
    original_subsequences = AdditiveSequence.subsequences
    original_bound = AdditiveSequenceSpace.davenport_bound.fget
    bound_queries = []
    caps = []

    def checked_subsequences(self, **kwargs):
        if self is sequence:
            caps.append(kwargs.get("max_length"))
            # Check before enumeration: an uncapped search would be enormous.
            assert kwargs.get("max_length") == expected_bound
        yield from original_subsequences(self, **kwargs)

    def counted_bound(self):
        bound_queries.append(self)
        return original_bound(self)

    monkeypatch.setattr(AdditiveSequence, "subsequences", checked_subsequences)
    monkeypatch.setattr(AdditiveSequenceSpace, "davenport_bound", property(counted_bound))
    solver = sequence.factorization_solver()
    assert caps == [expected_bound]
    assert bound_queries == [space]
    assert len(solver._atoms) == (1 if expected_bound == 2 else 3)
    if bound is not None:
        assert space._arithmetic_context is None


def test_local_search_has_no_ambient_cutoff(integer_problem, monkeypatch):
    _, sequence, _, expected = integer_problem
    original = AdditiveSequence.subsequences
    caps = []

    def checked(self, **kwargs):
        if self is sequence:
            caps.append(kwargs.get("max_length"))
        yield from original(self, **kwargs)

    monkeypatch.setattr(AdditiveSequence, "subsequences", checked)
    assert set(sequence.factorizations()) == expected
    assert caps == [None]  # the finite input is exhausted by the product iterator


@pytest.mark.parametrize(
    "terms, lengths", [([], {0}), ([1], set()), ([-1, 2], set()), ([-2, 1, 1], {1})]
)
def test_local_search_handles_empty_atomic_and_nonzero_sums(terms, lengths):
    space = AdditiveSequenceSpace(IntegerParent())
    assert space(terms).length_set() == lengths


def test_local_search_needs_no_finiteness_protocol():
    class MinimalParent:
        def __call__(self, term):
            return int(term)

        def zero(self):
            return 0

    space = AdditiveSequenceSpace(MinimalParent())
    assert space([-1, -1, 1, 1]).length_set() == {2}


def test_local_search_supports_a_finite_parent_without_order_metadata():
    class NonIterableC3(IntegerParent):
        def __call__(self, term):
            return int(term) % 3

        def is_finite(self):
            return True

        def add(self, left, right):
            return (left + right) % 3

    space = AdditiveSequenceSpace(NonIterableC3())
    assert space([1, 1, 1, 2, 2, 2]).length_set() == {2, 3}


@pytest.mark.parametrize("bound", [None, 2])
@pytest.mark.parametrize("empty", [False, True])
def test_catalogue_jsonl_with_optional_bound(tmp_path, bound, empty):
    space = AdditiveSequenceSpace(IntegerParent(), davenport_bound=bound)
    atoms = [] if empty else [space([-1, 1])]
    catalogue = AtomCatalogue(
        space, atoms, annotations=[(atom, {"tag": "pair"}) for atom in atoms]
    )
    path = tmp_path / "atoms.jsonl"
    assert catalogue.to_jsonl(path) == len(atoms)
    lines = path.read_text().splitlines()
    metadata = json.loads(lines[0])
    assert metadata["davenport_bound"] == bound
    assert metadata["digest"] == catalogue.digest()
    for verify in (False, True):
        loaded = AtomCatalogue.from_jsonl(path, space, verify=verify)
        assert loaded.atoms == catalogue.atoms
        assert loaded.annotations == catalogue.annotations
    del metadata["davenport_bound"]
    path.write_text("\n".join([json.dumps(metadata), *lines[1:]]) + "\n")
    assert AtomCatalogue.from_jsonl(path, space).atoms == catalogue.atoms


def test_catalogue_export_preserves_an_automatic_bound(tmp_path, monkeypatch):
    group = FiniteAdditiveGroup.cyclic_product(2, 2)
    space = AdditiveSequenceSpace(group)
    catalogue = AtomCatalogue(space, [space([(0, 1), (0, 1)])])

    def forbidden(*args, **kwargs):
        raise AssertionError("export must not start an exact search")

    monkeypatch.setattr("zero_sum_sequences.invariants.zero_sum_constant", forbidden)
    path = tmp_path / "atoms.jsonl"
    catalogue.to_jsonl(path)
    metadata = json.loads(path.read_text().splitlines()[0])
    assert metadata["davenport_bound"] == 3
    assert "addition_table" not in group._group_arithmetic().__dict__
    assert group._automorphism_group is None


@pytest.mark.parametrize("error_type", [ValueError, RuntimeError])
def test_unrelated_bound_errors_are_not_hidden(
    integer_problem, tmp_path, monkeypatch, error_type
):
    _, sequence, catalogue, _ = integer_problem

    def broken(self):
        raise error_type("broken bound provider")

    monkeypatch.setattr(AdditiveSequenceSpace, "davenport_bound", property(broken))
    with pytest.raises(error_type, match="broken bound provider"):
        sequence.factorization_solver()
    with pytest.raises(error_type, match="broken bound provider"):
        catalogue.to_jsonl(tmp_path / "atoms.jsonl")
