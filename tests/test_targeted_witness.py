"""Single-witness queries must not require the complete remainder DAG."""

from itertools import combinations_with_replacement

import pytest

from zero_sum_sequences import AdditiveSequenceSpace, AtomCatalogue, FactorizationSolver

from groups import cyclic_group


def fail_full_solve(*args, **kwargs):
    pytest.fail("a targeted witness query requested a complete solve")


@pytest.mark.parametrize("warmup", ["none", "existence", "maximum", "length_set"])
def test_witness_queries_use_targeted_search_or_existing_cache(monkeypatch, warmup):
    space = AdditiveSequenceSpace(cyclic_group(4), davenport_bound=4)
    sequence = space([1] * 4 + [3] * 4)
    solver = FactorizationSolver(sequence)
    if warmup == "existence":
        assert solver.has_factorization_of_length(2)
    elif warmup == "maximum":
        assert solver.maximum_factorization_length() == 4
    elif warmup == "length_set":
        assert solver.length_set() == {2, 4}
    if warmup in {"maximum", "length_set"}:
        # Cached transitions should suffice as well as cached length answers.
        monkeypatch.setattr(solver, "_state_transitions", fail_full_solve)

    monkeypatch.setattr(solver, "_solve_lengths", fail_full_solve)
    monkeypatch.setattr(solver, "_build_state_graph", fail_full_solve)

    assert solver.factorization_witness(3) is None  # A hole, not a term bound.
    assert solver.factorization_witness(10**30) is None
    assert solver.factorization_witness(0) is None
    assert solver.factorization_witness(4) == (space([1, 3]),) * 4
    assert solver.factorization_witness(2) == (space([1] * 4), space([3] * 4))

    assert solver.factorization_witness(
        2, minimum_matching_factors=1, factor_predicate=lambda atom: len(atom) == 2
    ) is None
    # Failure for one predicate must not poison an unconstrained or differently
    # constrained query of the same length.
    witness = solver.factorization_witness(
        2, minimum_matching_factors=2, factor_predicate=lambda atom: len(atom) == 4
    )
    assert witness is not None
    assert len(witness) == 2
    assert sum(witness, space()) == sequence
    assert solver.has_factorization_of_length(2)
    assert solver.factorization_witness(2) == witness

    if warmup in {"none", "existence"}:
        assert solver._transitions is None
    if warmup != "length_set":
        assert solver._length_bits is None


@pytest.mark.parametrize("modulus", [3, 4])
def test_targeted_witnesses_match_exhaustive_enumeration(modulus):
    space = AdditiveSequenceSpace(cyclic_group(modulus), davenport_bound=modulus)
    catalogue = space.enumerate_atom_catalogue()
    sequences = {space([1])}  # Also test a sequence with no factorization.
    sequences.update(
        sum(factors, space())
        for count in range(4)
        for factors in combinations_with_replacement(catalogue, count)
    )
    predicates = (
        lambda atom: False,
        lambda atom: len(atom) == 2,
        lambda atom: atom.multiplicity(1) > 0,
        lambda atom: True,
    )
    for sequence in sequences:
        expected = tuple(sequence.factorizations(atom_catalogue=catalogue))
        solver = FactorizationSolver(sequence, atom_catalogue=catalogue)
        saved_witnesses = {}
        for count in range(len(sequence) // 2 + 2):
            witness = solver.factorization_witness(count)
            assert (witness is not None) == any(len(f) == count for f in expected)
            saved_witnesses[count] = witness
            if witness is not None:
                assert len(witness) == count
                assert sum(witness, space()) == sequence
            for predicate in predicates:
                for matches in range(count + 1):
                    witness = solver.factorization_witness(
                        count,
                        minimum_matching_factors=matches,
                        factor_predicate=predicate,
                    )
                    possible = any(
                        len(f) == count and sum(map(predicate, f)) >= matches
                        for f in expected
                    )
                    assert (witness is not None) == possible
                    if witness is not None:
                        assert len(witness) == count
                        assert sum(witness, space()) == sequence
                        assert all(atom.is_atom() for atom in witness)
                        assert sum(map(predicate, witness)) >= matches
        assert solver._transitions is None
        assert solver._length_bits is None
        # Filling the full cache later does not change deterministic witnesses.
        assert solver.length_set() == {len(f) for f in expected}
        for count, witness in saved_witnesses.items():
            assert solver.factorization_witness(count) == witness


def test_targeted_witness_empty_sequence(monkeypatch):
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    solver = FactorizationSolver(space())
    monkeypatch.setattr(solver, "_solve_lengths", fail_full_solve)
    monkeypatch.setattr(solver, "_build_state_graph", fail_full_solve)

    assert solver.factorization_witness(0) == ()
    assert solver.factorization_witness(1) is None


def test_long_targeted_witness_does_not_recurse_or_build_full_graph(monkeypatch):
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    atom = space([1, 1, 1])
    solver = FactorizationSolver(5000 * atom, atom_catalogue=AtomCatalogue(space, [atom]))
    monkeypatch.setattr(solver, "_solve_lengths", fail_full_solve)
    monkeypatch.setattr(solver, "_build_state_graph", fail_full_solve)

    assert solver.factorization_witness(
        5000,
        minimum_matching_factors=5000,
        factor_predicate=lambda factor: factor == atom,
    ) == (atom,) * 5000
    assert solver._transitions is None
    assert solver._length_bits is None
