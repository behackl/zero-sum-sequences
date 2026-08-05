import pytest

from benchmarks.factorization_cases import factorization_benchmark_cases  # noqa: E402
from zero_sum_sequences import FactorizationSolver


CASES = factorization_benchmark_cases()


@pytest.mark.parametrize("case", CASES, ids=lambda case: case.name)
def test_solver_on_benchmark_corpus(case):
    solver = FactorizationSolver(case.sequence)

    assert solver.length_set() == set(case.expected_lengths)
    assert solver.statistics.states >= 1
    expected_bound = 7 if case.name.startswith("rank-three") else 3
    assert solver.davenport_bound == expected_bound


def test_very_long_factorization_enumeration_is_iterative():
    [case] = [case for case in CASES if case.name == "c3-pure-power-5000"]
    solver = FactorizationSolver(case.sequence)
    [factorization] = solver.factorizations()

    assert len(factorization) == 5000


@pytest.mark.parametrize(
    "case",
    [
        case
        for case in CASES
        if case.expected_factorizations is not None and case.tier != "very-long"
    ],
    ids=lambda case: case.name,
)
def test_safe_benchmark_factorization_counts(case):
    solver = FactorizationSolver(case.sequence)

    assert sum(1 for _ in solver.factorizations()) == case.expected_factorizations
