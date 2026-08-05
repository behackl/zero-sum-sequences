"""Benchmark the memoized factorization implementation."""

from __future__ import annotations

import argparse
import platform
from collections.abc import Sequence
from time import perf_counter

from zero_sum_sequences import FactorizationSolver

from .factorization_cases import factorization_benchmark_cases


def elapsed(callable_):
    start = perf_counter()
    result = callable_()
    return result, perf_counter() - start


def benchmark_case(
    case, *, enumerate_factorizations: bool, repeats: int
) -> dict[str, object]:
    best_run = None
    for _ in range(repeats):
        solver, preparation_seconds = elapsed(
            lambda: FactorizationSolver(case.sequence)
        )
        lengths, solve_seconds = elapsed(solver.length_set)
        if best_run is None or preparation_seconds + solve_seconds < sum(best_run[1:]):
            best_run = (solver, preparation_seconds, solve_seconds)
    assert best_run is not None
    solver, preparation_seconds, solve_seconds = best_run
    lengths = solver.length_set()
    if lengths != set(case.expected_lengths):
        raise AssertionError(
            f"{case.name}: expected {sorted(case.expected_lengths)}, "
            f"got {sorted(lengths)}"
        )
    statistics = solver.statistics

    factorization_count: int | None = None
    enumeration_seconds: float | None = None
    if enumerate_factorizations and case.enumerate_factorizations:
        factorizations, enumeration_seconds = elapsed(
            lambda: list(solver.factorizations())
        )
        factorization_count = len(factorizations)
        if factorization_count != case.expected_factorizations:
            raise AssertionError(
                f"{case.name}: expected {case.expected_factorizations} "
                f"factorizations, got {factorization_count}"
            )

    return {
        "name": case.name,
        "tier": case.tier,
        "terms": len(case.sequence),
        "lengths": len(lengths),
        "atoms": statistics.candidate_atoms,
        "states": statistics.states,
        "transitions": statistics.transitions,
        "prepare_ms": preparation_seconds * 1000,
        "solve_ms": solve_seconds * 1000,
        "factorizations": factorization_count,
        "enumerate_ms": (
            None if enumeration_seconds is None else enumeration_seconds * 1000
        ),
    }


def format_value(value: object) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def print_markdown(rows: list[dict[str, object]]) -> None:
    columns = (
        "name",
        "tier",
        "terms",
        "lengths",
        "atoms",
        "states",
        "transitions",
        "prepare_ms",
        "solve_ms",
        "factorizations",
        "enumerate_ms",
    )
    print("| " + " | ".join(columns) + " |")
    print("| " + " | ".join("---" for _ in columns) + " |")
    for row in rows:
        print("| " + " | ".join(format_value(row[column]) for column in columns) + " |")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--enumerate",
        action="store_true",
        help="also enumerate unique factorizations for cases marked as safe",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="number of timing repetitions; the minimum is reported",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    if arguments.repeats < 1:
        raise SystemExit("--repeats must be positive")
    from sage.version import version as sage_version

    print(f"SageMath {sage_version}; Python {platform.python_version()}")
    print()
    cases = factorization_benchmark_cases()
    rows = [
        benchmark_case(
            case,
            enumerate_factorizations=arguments.enumerate,
            repeats=arguments.repeats,
        )
        for case in cases
    ]
    print_markdown(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
