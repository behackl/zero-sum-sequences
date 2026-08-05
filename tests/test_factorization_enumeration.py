import itertools
from functools import reduce
from operator import add

import networkx as nx
import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AtomCatalogue,
    FactorizationSolver,
)

from groups import cyclic_group

C3 = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)


def c3(*values: int):
    return C3(values)


def product(factors, space=C3):
    return reduce(add, factors, space())


def reference_atoms(space, length):
    zero = space.base_parent.zero()
    for terms in itertools.combinations_with_replacement(
        tuple(space.base_parent), length
    ):
        sequence = space(terms)
        if sequence.is_atom() and zero not in sequence:
            yield sequence


def naive_factorizations(sequence, atoms, minimum_index=0):
    """Reference enumeration with nondecreasing atom indices."""

    if not sequence:
        yield ()
        return
    for index in range(minimum_index, len(atoms)):
        atom = atoms[index]
        if atom.divides(sequence):
            for rest in naive_factorizations(sequence - atom, atoms, index):
                yield (atom, *rest)


def factorization_profile(factorization, atoms):
    return tuple(factorization.count(atom) for atom in atoms)


def test_remainder_digraph_records_states_and_atom_edges():
    sequence = c3(1, 1, 1, 2, 2, 2)
    empty = C3()
    graph = sequence.factorization_digraph()

    assert graph.is_directed()
    assert nx.is_directed_acyclic_graph(graph)
    assert sequence in graph
    assert empty in graph
    assert graph.out_degree(empty) == 0
    assert {
        len(path) - 1 for path in nx.all_simple_paths(graph, sequence, empty)
    } == sequence.length_set()
    assert all(
        source - atom == target
        for source, target, atom in graph.edges(data="atom")
    )


def test_factorizations_length_sets_and_witnesses():
    atom_1 = c3(1, 1, 1)
    atom_2 = c3(2, 2, 2)
    pair = c3(1, 2)
    sequence = atom_1 + atom_2

    assert set(sequence.factorizations()) == {
        (atom_1, atom_2),
        (pair, pair, pair),
    }
    assert sequence.length_set() == {2, 3}

    witnesses = sequence.factorization_witnesses()
    assert set(witnesses) == {2, 3}
    assert all(product(factors) == sequence for factors in witnesses.values())
    assert all(all(atom.is_atom() for atom in factors) for factors in witnesses.values())


def test_empty_atom_nonzero_and_zero_term_cases():
    empty = C3()
    atom = c3(1, 1, 1)
    nonzero = c3(1)
    zero_times_atom = c3(0) + atom

    assert empty.length_set() == {0}
    assert list(empty.factorizations()) == [()]
    assert atom.length_set() == {1}
    assert nonzero.length_set() == set()
    with pytest.raises(ValueError, match="without zero terms"):
        zero_times_atom.length_set()


def test_memoized_enumeration_agrees_with_naive_on_small_products():
    atoms = tuple(
        atom
        for length in range(1, 4)
        for atom in reference_atoms(C3, length)
    )
    sequences = {
        product(factors)
        for number_of_factors in range(4)
        for factors in itertools.combinations_with_replacement(atoms, number_of_factors)
    }

    for sequence in sequences:
        expected = {
            factorization_profile(factorization, atoms)
            for factorization in naive_factorizations(sequence, atoms)
        }
        actual_factorizations = list(sequence.factorizations())
        actual_profiles = [
            factorization_profile(factorization, atoms)
            for factorization in actual_factorizations
        ]
        assert len(actual_profiles) == len(set(actual_profiles))
        assert set(actual_profiles) == expected


def test_precomputed_atom_catalogue_can_be_reused():
    atoms = tuple(
        atom
        for length in range(1, 4)
        for atom in reference_atoms(C3, length)
    )
    catalogue = AtomCatalogue(C3, atoms)
    sequence = c3(1, 1, 1, 2, 2, 2)

    assert sequence.length_set(atom_catalogue=catalogue) == {2, 3}


def test_targeted_witness_search_agrees_with_exhaustive_factorizations():
    atoms = tuple(
        atom
        for length in range(1, 4)
        for atom in reference_atoms(C3, length)
    )
    catalogue = AtomCatalogue(C3, atoms)
    sequences = {
        product(factors)
        for factor_count in range(4)
        for factors in itertools.combinations_with_replacement(atoms, factor_count)
    }

    for sequence in sequences:
        solver = FactorizationSolver(sequence, atom_catalogue=catalogue)
        factorizations = list(solver.factorizations())
        for factor_count in range(4):
            for required_short in range(factor_count + 1):
                expected = any(
                    len(factorization) == factor_count
                    and sum(len(factor) <= 2 for factor in factorization)
                    >= required_short
                    for factorization in factorizations
                )
                witness = solver.factorization_witness(
                    factor_count,
                    minimum_matching_factors=required_short,
                    factor_predicate=(
                        (lambda factor: len(factor) <= 2)
                        if required_short
                        else None
                    ),
                )
                assert (witness is not None) == expected
