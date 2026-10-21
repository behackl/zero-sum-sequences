"""Timing runner for catalogue, atom, symmetry and solver operations.

Run from the repository root::

    python -m benchmarks.benchmark_operations          # quick cases
    python -m benchmarks.benchmark_operations --large  # adds the C_3^3 catalogue

The workloads follow the computations of the companion article over
``C_3^3`` (27 911 atoms, 11 232 automorphisms). Times are the best of
``--repeat`` runs, in milliseconds.
"""

from __future__ import annotations

import argparse
import random
import time
from collections.abc import Callable

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    FactorizationSolver,
    FiniteAdditiveGroup,
)


def best_time(function: Callable[[], object], repeat: int) -> float:
    best = float("inf")
    for _ in range(repeat):
        start = time.perf_counter()
        function()
        best = min(best, time.perf_counter() - start)
    return best * 1e3


def space(*moduli: int, bound: int | None = None) -> AdditiveSequenceSpace:
    return AdditiveSequenceSpace(
        FiniteAdditiveGroup.cyclic_product(*moduli), davenport_bound=bound
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--large", action="store_true", help="include the C_3^3 catalogue")
    parser.add_argument("--repeat", type=int, default=3)
    args = parser.parse_args()

    rng = random.Random(2026)
    cases: list[tuple[str, Callable[[], object], int]] = [
        ("catalogue C_2^4 (323 atoms)", lambda: space(2, 2, 2, 2).enumerate_atom_catalogue(), args.repeat),
    ]
    if args.large:
        cases.append(("catalogue C_3 + C_6", lambda: space(3, 6).enumerate_atom_catalogue(), 1))
        cases.append(("catalogue C_3^3 (27 911 atoms)", lambda: space(3, 3, 3).enumerate_atom_catalogue(), 1))

    # Shared C_3^3 data: catalogue, automorphism group, sample inputs.
    c33 = space(3, 3, 3)
    catalogue = c33.enumerate_atom_catalogue()
    group = c33.automorphism_group()
    atoms = list(catalogue)
    long_atoms = [atom for atom in atoms if len(atom) >= 5]
    sample_atoms = rng.sample(atoms, 2000)
    products = [rng.choice(long_atoms) * rng.choice(long_atoms) * rng.choice(atoms) for _ in range(100)]
    pairs = [rng.choice(long_atoms) * rng.choice(long_atoms) for _ in range(5)]
    seven = [atom for atom in atoms if len(atom) == 7][:10]
    automorphisms = group.elements[:2000]
    sequence = pairs[0]

    cases += [
        ("is_atom, 2000 C_3^3 atoms", lambda: [atom.is_atom() for atom in sample_atoms], args.repeat),
        ("catalogue-backed length sets, 100 triple products",
         lambda: [FactorizationSolver(p, atom_catalogue=catalogue).length_set() for p in products], args.repeat),
        ("apply 2000 automorphisms to a pair product",
         lambda: [group.apply(element, sequence) for element in automorphisms], args.repeat),
        ("stabilizer of 10 length-7 atoms", lambda: [group.stabilizer(atom) for atom in seven], args.repeat),
        ("canonical form of 5 pair products", lambda: [group.canonical(pair) for pair in pairs], 1),
    ]

    width = max(len(name) for name, _, _ in cases)
    print(f"{'case':<{width}}  {'ms':>10}")
    for name, function, repeat in cases:
        print(f"{name:<{width}}  {best_time(function, repeat):10.1f}", flush=True)


if __name__ == "__main__":
    main()
