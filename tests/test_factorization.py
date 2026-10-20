import math
from fractions import Fraction

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AtomCatalogue,
    FactorizationSolver,
)

from groups import IndexableInteger, cyclic_group


def test_solver_uses_the_space_davenport_bound():
    complete = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    truncated = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=2)

    assert complete(1, 1, 1).length_set() == {1}
    assert truncated(1, 1, 1).length_set() == set()


def test_targeted_factorization_witness_honors_generic_predicate():
    space = AdditiveSequenceSpace(cyclic_group(4), davenport_bound=4)
    short = space(2, 2)
    long = space(1, 1, 1, 1)
    sequence = short ** 2 * long
    catalogue = AtomCatalogue(space, [short, long])
    solver = FactorizationSolver(sequence, atom_catalogue=catalogue)

    witness = solver.factorization_witness(
        3,
        minimum_matching_factors=2,
        factor_predicate=lambda factor: len(factor) <= 3,
    )

    assert witness is not None
    assert math.prod(witness, start=space()) == sequence
    assert sum(len(factor) <= 3 for factor in witness) >= 2
    assert (
        solver.factorization_witness(
            3,
            minimum_matching_factors=3,
            factor_predicate=lambda factor: len(factor) <= 3,
        )
        is None
    )
    assert solver.factorization_witness(2) is None


def test_targeted_factorization_witness_accepts_indexable_integers():
    space = AdditiveSequenceSpace(
        cyclic_group(3),
        davenport_bound=IndexableInteger(3),
    )
    sequence = space(1, 1, 1, 2, 2, 2)
    solver = FactorizationSolver(sequence)

    witness = solver.factorization_witness(
        IndexableInteger(3),
        minimum_matching_factors=IndexableInteger(3),
        factor_predicate=lambda factor: len(factor) == 2,
    )

    assert witness == (space(1, 2),) * 3
    assert solver.has_factorization_of_length(IndexableInteger(3))
    assert list(
        solver.factorizations(factor_count=IndexableInteger(3))
    ) == [(space(1, 2),) * 3]
    with pytest.raises(ValueError, match="factor count"):
        solver.factorization_witness(Fraction(3, 2))
    with pytest.raises(ValueError, match="minimum matching factors"):
        solver.factorization_witness(
            3,
            minimum_matching_factors=Fraction(3, 2),
        )


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"factor_count": -1}, "factor count"),
        ({"factor_count": True}, "factor count"),
        (
            {"factor_count": 2, "minimum_matching_factors": 3},
            "cannot exceed",
        ),
        (
            {"factor_count": 2, "minimum_matching_factors": 1},
            "factor_predicate",
        ),
    ],
)
def test_targeted_factorization_witness_rejects_invalid_constraints(
    arguments, message
):
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    solver = FactorizationSolver(space())

    with pytest.raises(ValueError, match=message):
        solver.factorization_witness(**arguments)


@pytest.mark.parametrize("factor_count", [-1, True, Fraction(3, 2)])
def test_fixed_length_queries_reject_invalid_factor_counts(factor_count):
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    sequence = space(1, 1, 1, 2, 2, 2)
    solver = FactorizationSolver(sequence)

    with pytest.raises(ValueError, match="factor count"):
        solver.has_factorization_of_length(factor_count)
    with pytest.raises(ValueError, match="factor count"):
        list(solver.factorizations(factor_count=factor_count))
    with pytest.raises(ValueError, match="factor count"):
        sequence.has_factorization_of_length(factor_count)
    with pytest.raises(ValueError, match="factor count"):
        list(sequence.factorizations(factor_count=factor_count))


def test_solver_requires_a_catalogue_from_the_same_space():
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    other_space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)

    with pytest.raises(TypeError, match="AtomCatalogue"):
        FactorizationSolver(space(0), atom_catalogue=[space(0)])
    with pytest.raises(TypeError, match="different spaces"):
        FactorizationSolver(
            space(1, 2),
            atom_catalogue=AtomCatalogue(other_space, [other_space(1, 2)]),
        )


@pytest.mark.parametrize(
    ("terms", "expected"),
    [([], 0), ([1], None), ([1, 1, 1, 2, 2, 2], 3)],
    ids=["empty", "unfactorable", "multiple-lengths"],
)
def test_maximum_only_avoids_length_sets_and_caches_its_answer(terms, expected, monkeypatch):
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    solver = FactorizationSolver(space.from_terms(terms))

    def forbidden():
        raise AssertionError("maximum-only queries must not compute length sets")

    monkeypatch.setattr(solver, "_solve_lengths", forbidden)
    assert solver.maximum_factorization_length() == expected
    monkeypatch.setattr(solver, "_build_state_graph", forbidden)
    assert solver.maximum_factorization_length() == expected


@pytest.mark.parametrize(
    "terms", [[], [1], [1, 1, 1, 2, 2, 2]],
    ids=["empty", "unfactorable", "multiple-lengths"],
)
def test_maximum_reuses_existing_length_sets(terms, monkeypatch):
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    solver = FactorizationSolver(space.from_terms(terms))
    lengths = solver.length_set()

    def forbidden():
        raise AssertionError("cached length data must not trigger a graph pass")

    monkeypatch.setattr(solver, "_build_state_graph", forbidden)
    assert solver.maximum_factorization_length() == max(lengths, default=None)


def test_maximum_then_length_set_returns_complete_lengths():
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    solver = FactorizationSolver(space(1, 1, 1, 2, 2, 2))
    assert solver.maximum_factorization_length() == 3
    assert solver.length_set() == {2, 3}
    assert solver.maximum_factorization_length() == 3
    assert solver.minimum_factorization_length() == 2


def test_maximum_of_a_long_atom_power_avoids_length_sets(monkeypatch):
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)
    solver = FactorizationSolver(space.from_multiplicities({1: 60_000}))

    def forbidden():
        raise AssertionError("long maximum-only queries must not compute length sets")

    monkeypatch.setattr(solver, "_solve_lengths", forbidden)
    assert solver.maximum_factorization_length() == 20_000


def test_factorization_rejects_sequences_with_the_term_zero():
    space = AdditiveSequenceSpace(cyclic_group(3), davenport_bound=3)

    with pytest.raises(ValueError, match="without the term 0"):
        FactorizationSolver(space(0, 1, 2))
