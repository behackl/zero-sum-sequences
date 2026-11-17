# Getting started

## Installation

The package requires Python 3.12 or newer:

```console
python -m pip install zero-sum-sequences
```

## Groups and sequences

We work over $G=C_2\oplus C_4$, whose elements are coordinate tuples. A
sequence space constructs sequences over $G$; as in the literature, they
are written multiplicatively.

```{doctest}
>>> from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup
>>> G = FiniteAdditiveGroup.cyclic_product(2, 4)
>>> S = AdditiveSequenceSpace(G)
>>> A = S((0, 1))**3 * S((1, 0), (1, 1))
>>> A
(0, 1)^3 · (1, 0) · (1, 1)
>>> len(A), A.sum(), A.is_atom()
(5, (0, 0), True)
```

$A$ is an atom, a minimal zero-sum sequence: its sum is zero, but no
proper nontrivial subsequence sums to zero.

## Factorizations

Every zero-sum sequence factors into atoms, usually in several ways. The
**length** of a factorization is its number of atoms. Take $B=-A$:

```{doctest}
>>> B = -A
>>> B
(0, 3)^3 · (1, 0) · (1, 3)
>>> sorted((A * B).length_set())
[2, 4, 5]
```

So $AB$ factors into two, four or five atoms, but not three. For several
questions about one sequence, keep a solver:

```{doctest}
>>> from zero_sum_sequences import FactorizationSolver
>>> solver = FactorizationSolver(A * B)
>>> solver.factorization_witness(3) is None
True
>>> for atom in solver.factorization_witness(4):
...     print(atom)
(0, 1) · (0, 3)
(0, 1) · (0, 3)
(0, 1) · (1, 0) · (1, 3)
(0, 3) · (1, 0) · (1, 1)
```

## All atoms and their orbits

$G$ has 38 atoms besides the sequence $0$. The longest have length
$D(G)=5$, the Davenport constant.

```{doctest}
>>> atoms = S.enumerate_atom_catalogue()
>>> len(atoms), max(map(len, atoms))
(38, 5)
```

Automorphisms of $G$ act on sequences term by term. The orbit of $A$
consists of the eight atoms of length five, and all 38 atoms fall into
eleven orbits.

```{doctest}
>>> len(A.orbit()), A.is_in_same_orbit(B)
(8, True)
>>> len(atoms.orbits())
11
```

## Zero-sum constants

The Davenport, eta and Erdős–Ginzburg–Ziv constants are computed exactly:

```{doctest}
>>> from zero_sum_sequences import (
...     davenport_constant, erdos_ginzburg_ziv_constant, eta_constant)
>>> davenport_constant(S)
5
>>> eta_constant(S)
6
>>> erdos_ginzburg_ziv_constant(S)
9
```

{py:func}`~zero_sum_sequences.zero_sum_constant` handles arbitrary sets of
allowed subsequence lengths and also returns the extremal sequences.

## Next steps

- The [tutorial](tutorial.md) continues these examples interactively in
  the browser, including the factorization graph and atoms on a prescribed
  support.
- [Definitions](definitions.md) fixes the notation.
- The [API reference](api/index.rst) documents every class and function.

All results are exact. Atom searches are limited by an upper bound for
$D(G)$, which the space derives automatically. If you supply your own
`davenport_bound=` or atom catalogue, results are complete only if it
suffices; see {py:class}`~zero_sum_sequences.FactorizationSolver`.

## SageMath

Compatible Sage parents, such as `GF(3)**2`, can be used directly in place
of `FiniteAdditiveGroup`. Sage is optional; on supported platforms, the
`sage` extra installs its runtime:

```console
python -m pip install 'zero-sum-sequences[sage]'
```

## Development

From a checkout, install the locked environment and run the tests:

```console
uv sync --locked
uv run --locked python -m pytest
```

Add `--extra sage` to `uv sync` for the Sage integration. The
[benchmark guide](https://github.com/behackl/zero-sum-sequences/blob/main/benchmarks/README.md)
describes the performance measurements.
