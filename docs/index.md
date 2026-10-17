# Zero-sum sequences

A Python library for computing with finite multisets over additive groups:
zero-sum atoms, factorizations and their lengths, automorphism orbits, and
zero-sum constants.

The package provides ordinary-Python group adapters, including products of
cyclic groups. Compatible parents from SageMath can also be used directly;
Sage is optional.

## A first calculation

Over $C_3$, the sequence $1^3 2^3$ has factorizations $(1^3)(2^3)$ and
$(1\cdot2)^3$, with two and three atom factors respectively:

```{doctest}
>>> from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup
>>> G = FiniteAdditiveGroup.cyclic_product(3)
>>> S = AdditiveSequenceSpace(G, davenport_bound=3)
>>> sequence = S([(1,)] * 3 + [(2,)] * 3)
>>> sequence.is_zero_sum(), sequence.is_atom()
(True, False)
>>> sorted(sequence.length_set())
[2, 3]
```

Sequences are immutable multisets: input order does not matter, and their
term count is distinct from the number of factors in a factorization.
The [getting-started guide](getting-started.md) explains the construction
and installation options.

## What you can compute

- **Sequences and group arithmetic:** construct and combine multisets, test
  zero-sumness and minimality, and work with additive orders and inverses.
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

These are explicit computations intended for small examples, not a
replacement for structural theorems or a large-scale computation service.
Exhaustive searches can grow quickly; the API documents their completeness
conditions and which queries are targeted.

## Explore the documentation

- **[Getting started](getting-started.md)** — installation, a small calculation,
  optional Sage support, and development setup.
- **[Tutorial](tutorial.md)** — a longer tour with examples and a graph
  visualization.
- **[Definitions](definitions.md)** — the notation and conventions used here.
- **[API reference](api/index.rst)** — parameters, examples, and computational
  notes, generated from the code's docstrings.
- **[References](references.md)** — mathematical background and attribution.

```{toctree}
:hidden:
:maxdepth: 2

getting-started
tutorial
definitions
api/index
references
```
