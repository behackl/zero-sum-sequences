"""Timing runner for catalogue, atom, factorization, symmetry and constant operations.

Run from the repository root::

    python -m benchmarks.benchmark_operations          # quick cases
    python -m benchmarks.benchmark_operations --large  # adds larger catalogues

The cases cover cyclic groups, groups of rank two, elementary abelian groups
and ``C_3^3`` (27 911 atoms, 11 232 automorphisms), with symmetric and
generic inputs. Times are the best of ``--repeat`` runs, in milliseconds.
"""

from __future__ import annotations

import argparse
import math
import random
import time
from collections.abc import Callable

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    FactorizationSolver,
    FiniteAdditiveGroup,
    erdos_ginzburg_ziv_constant,
    eta_constant,
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


def nonzero(sequence_space: AdditiveSequenceSpace) -> list:
    zero = sequence_space.base_parent.zero()
    return [term for term in sequence_space.base_parent if term != zero]


def group_name(moduli: tuple[int, ...]) -> str:
    if len(moduli) > 1 and len(set(moduli)) == 1:
        return f"C_{moduli[0]}^{len(moduli)}"
    return " + ".join(f"C_{m}" for m in moduli)


def fresh_group(*moduli: int):
    sequence_space = space(*moduli)
    return sequence_space, sequence_space.automorphism_group()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--large", action="store_true", help="include larger catalogues")
    parser.add_argument("--repeat", type=int, default=3)
    args = parser.parse_args()
    rng = random.Random(2026)
    repeat = args.repeat
    cases: list[tuple[str, Callable[[], object], int]] = []

    # Catalogues: complete enumeration of all atoms.
    for moduli in [(12,), (2, 8), (4, 4), (2, 2, 2, 2)]:
        name = group_name(moduli)
        cases.append((f"catalogue {name}", lambda m=moduli: space(*m).enumerate_atom_catalogue(), repeat))
    if args.large:
        for moduli in [(5, 5), (3, 6), (3, 3, 3)]:
            name = group_name(moduli)
            cases.append((f"catalogue {name}", lambda m=moduli: space(*m).enumerate_atom_catalogue(), 1))

    # Shared C_3^3 data.
    c33 = space(3, 3, 3)
    catalogue = c33.enumerate_atom_catalogue()
    group = c33.automorphism_group()
    atoms = list(catalogue)
    long_atoms = [atom for atom in atoms if len(atom) >= 5]

    # Atom tests.
    c55 = space(5, 5)
    random_c55 = [c55.from_terms(rng.choices(list(c55.base_parent), k=rng.randrange(1, 12))) for _ in range(3000)]
    c200 = space(200)
    powers = [c200((1,)) ** 200, c200((3,)) ** 200, c200((1,)) ** 100 * c200((199,)) ** 100]
    sample_atoms = rng.sample(atoms, 2000)
    cases += [
        ("is_atom, 3000 random sequences over C_5 + C_5", lambda: [s.is_atom() for s in random_c55], repeat),
        ("is_atom, high multiplicities over C_200", lambda: [s.is_atom() for s in powers], repeat),
        ("is_atom, 2000 atoms over C_3^3", lambda: [atom.is_atom() for atom in sample_atoms], repeat),
    ]

    # Length sets.
    c55_atoms = list(c55.enumerate_atom_catalogue())
    c55_products = [math.prod(rng.sample(c55_atoms, 3), start=c55()) for _ in range(15)]
    c18 = space(18)
    c18_sequence = c18((1,)) ** 18 * c18((5,)) ** 18 * c18((7,)) ** 18
    products = [rng.choice(long_atoms) * rng.choice(long_atoms) * rng.choice(atoms) for _ in range(100)]
    cases += [
        ("length sets, 15 atom triples over C_5 + C_5 (no catalogue)",
         lambda: [FactorizationSolver(p).length_set() for p in c55_products], repeat),
        ("length set of 1^18 5^18 7^18 over C_18", lambda: FactorizationSolver(c18_sequence).length_set(), repeat),
        ("length sets, 100 atom triples over C_3^3 (catalogue)",
         lambda: [FactorizationSolver(p, atom_catalogue=catalogue).length_set() for p in products], repeat),
    ]

    # Symmetry: generic inputs have large orbits, symmetric inputs small ones.
    c2222, group_2222 = fresh_group(2, 2, 2, 2)
    c55_group = c55.automorphism_group()
    generic_2222 = [c2222.from_terms(rng.choices(nonzero(c2222), k=8)) for _ in range(10)]
    everything = c2222.from_terms(nonzero(c2222))
    symmetric_2222 = [everything, everything ** 2]
    generic_c55 = [c55.from_terms(rng.choices(nonzero(c55), k=8)) for _ in range(10)]
    classified_c55 = [c55.from_terms(rng.choices(nonzero(c55), k=rng.randrange(2, 7))) for _ in range(300)]
    pairs = [rng.choice(long_atoms) * rng.choice(long_atoms) for _ in range(5)]
    seven = [atom for atom in atoms if len(atom) == 7][:10]
    automorphisms = group.elements[:2000]
    cases += [
        ("apply 2000 automorphisms over C_3^3", lambda: [group.apply(e, pairs[0]) for e in automorphisms], repeat),
        ("canonical forms, 10 generic sequences over C_2^4",
         lambda: [group_2222.canonical(s) for s in generic_2222], repeat),
        ("canonical forms, symmetric sequences over C_2^4",
         lambda: [group_2222.canonical(s) for s in symmetric_2222], repeat),
        ("canonical forms, 10 generic sequences over C_5 + C_5",
         lambda: [c55_group.canonical(s) for s in generic_c55], repeat),
        ("canonical forms, 5 atom pairs over C_3^3", lambda: [group.canonical(p) for p in pairs], repeat),
        ("stabilizers of 10 length-7 atoms over C_3^3", lambda: [group.stabilizer(a) for a in seven], repeat),
        ("orbit classes of 300 sequences over C_5 + C_5",
         lambda: c55_group.orbit_representatives(classified_c55), repeat),
    ]

    # Zero-sum constants (searches with symmetry reduction).
    cases += [
        ("eta(C_3^3)", lambda: eta_constant(space(3, 3, 3)), 1),
        ("s(C_2^4)", lambda: erdos_ginzburg_ziv_constant(space(2, 2, 2, 2)), 1),
        ("s(C_4 + C_4)", lambda: erdos_ginzburg_ziv_constant(space(4, 4)), 1),
    ]

    width = max(len(name) for name, _, _ in cases)
    print(f"{'case':<{width}}  {'ms':>10}")
    for name, function, runs in cases:
        print(f"{name:<{width}}  {best_time(function, runs):10.2f}", flush=True)


if __name__ == "__main__":
    main()
