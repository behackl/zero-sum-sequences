# Zero-sum sequences

A Python library for computing with sequences over finite additive groups:
minimal zero-sum sequences (atoms), factorizations and their length sets,
automorphism orbits, and zero-sum constants.

The package provides ordinary-Python group adapters, including products of
cyclic groups. Compatible parents from SageMath can also be used directly;
Sage is optional.

## A first calculation

Over $C_2\oplus C_4$, take the atom $A=(0,1)^3\cdot(1,0)\cdot(1,1)$ of
length five and its negative $B=-A$. Sequences are written
multiplicatively, as in the literature:

```{doctest}
>>> from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup
>>> G = FiniteAdditiveGroup.cyclic_product(2, 4)
>>> S = AdditiveSequenceSpace(G)
>>> A = S((0, 1))**3 * S((1, 0), (1, 1))
>>> B = -A
>>> B
(0, 3)^3 · (1, 0) · (1, 3)
>>> A.is_atom(), B.is_atom()
(True, True)
>>> sorted((A * B).length_set())
[2, 4, 5]
```

The product $AB$ has length ten, but its factorizations consist of two,
four, or five atoms. There is no factorization into three atoms: a length
set need not be an interval.

**New here?** [Getting started](getting-started.md) continues this example
step by step. The [tutorial](tutorial.md) is an interactive notebook that
runs in your browser, no installation needed.

## What you can compute

- **Sequences and group arithmetic:** construct and multiply sequences, test
  for zero-sum and minimal zero-sum sequences, and work with additive orders
  and inverses.
  See [Groups and sequences](api/sequences.rst).
- **Factorizations and relations:** compute length sets and extrema, find
  targeted witnesses, enumerate unordered factorizations, compare them by
  cancellation and distance, and visualize remainder graphs with NetworkX.
  See [Atoms and factorizations](api/factorizations.rst).
- **Catalogues and repeated queries:** enumerate atoms on a full or restricted
  support, keep stable indices, and exchange JSONL data with
  {py:class}`~zero_sum_sequences.AtomCatalogue`. Reuse factorization results
  through {py:class}`~zero_sum_sequences.FactorizationCache`.
- **Symmetry:** compute automorphism orbits, canonical representatives,
  stabilizers, and transporters. See [Automorphisms and orbits](api/symmetry.rst).
- **Zero-sum constants:** compute Davenport, Erdős–Ginzburg–Ziv and eta
  constants, or search with a chosen finite set of allowed subsequence
  lengths and inspect the extremal sequences. See
  [Zero-sum constants](api/constants.rst).

These are exact computations for small groups, not a replacement for
structural theorems. [Scope and performance](performance.md) describes
what is feasible.

## Explore the documentation

- **[Getting started](getting-started.md)** — installation and a first
  walk through the main features.
- **[Tutorial](tutorial.md)** — an interactive notebook in the browser,
  covering every main feature.
- **[How-to guides](how-to/index.md)** — recipes for common research tasks.
- **[Definitions](definitions.md)** — the notation and conventions used here.
- **[Scope and performance](performance.md)** — completeness of results and
  what is feasible.
- **[API reference](api/index.rst)** — parameters, examples, and computational
  notes, generated from the code's docstrings.
- **[References](references.md)** — mathematical background, attribution, and
  how to cite the package.
- **[Changelog](changelog.md)** — changes in each release.

```{toctree}
:hidden:
:maxdepth: 2

getting-started
tutorial
how-to/index
definitions
performance
api/index
references
changelog
```
