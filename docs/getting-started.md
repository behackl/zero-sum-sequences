# Getting started

## Installation

The package requires Python 3.12 or newer. These examples describe the
development version, including invariant functions not yet part of the
published 0.2.0 release. To run all of them, install from a checkout of
`zero-sum-constants` or `docs/shibuya`:

```console
python -m pip install --editable .
```

The published version is available separately:

```console
python -m pip install zero-sum-sequences
```

## A small sequence space

Build the cyclic group $C_3$ and a constructor for its sequences. Cyclic
products use coordinate tuples, even in rank one.

```{doctest}
>>> from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup
>>> S = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
>>> sample = S([(1,)] * 3)
>>> len(sample), sample.is_zero_sum(), sample.is_atom()
(3, True, True)
```

An atom is a minimal nonempty zero-sum sequence. An input's term order is
irrelevant: these are finite multisets, not ordered lists.

## Factorizations and their lengths

The sequence $1^3 2^3$ has factorizations $(1^3)(2^3)$ and $(1\cdot2)^3$,
of lengths 2 and 3.

```{doctest}
>>> sequence = S([(1,)] * 3 + [(2,)] * 3)
>>> sorted(sequence.length_set())
[2, 3]
```

The sequence has six terms; the numbers 2 and 3 count **atom factors**, not
terms. See [Definitions](definitions.md) for this distinction and the
identity-term convention.

## A zero-sum constant

For allowed subsequence lengths $L=\{3\}$, compute $s_L(C_3)$:

```{doctest}
>>> from zero_sum_sequences import zero_sum_constant
>>> result = zero_sum_constant(S, {3}, symmetry=False)
>>> result.value
5
>>> all(len(extremal) == 4 for extremal in result.extremal)
True
```

The result also contains the extremal sequences. Here `symmetry=False`
keeps all of them; symmetry reduction instead keeps one representative
per orbit.

## Bounds and completeness

The space derives a certified Davenport upper bound lazily when needed:

```{doctest}
>>> S.davenport_bound
3
```

An explicit `davenport_bound=...` remains authoritative for enumeration
and local atom discovery. If too small, it can make catalogues and
factorization results incomplete; an incomplete supplied catalogue has
the same limitation. The atom predicate and exact zero-sum constants
ignore enumeration bounds. Automatic bound resolution does not launch an
exhaustive exact search.

For the precise bound policy, see
{py:attr}`zero_sum_sequences.AdditiveSequenceSpace.davenport_bound`;
for the invariant search, see {py:func}`zero_sum_sequences.zero_sum_constant`.

## Optional Sage support

Compatible Sage parents can be used directly. If you do not already have a
Sage environment, the `sage` extra installs its runtime on supported platforms:

```console
python -m pip install --editable '.[sage]'
```

Sage is not needed for the ordinary-Python examples or the tutorial.

## Development

From a checkout, install the locked development environment and run the tests:

```console
uv sync --locked
uv run --locked python -m pytest
```

The default `dev` group includes the test tools, tutorial dependencies, and
JupyterLab. Add `--extra sage` to `uv sync` to include the Sage integration.
Without uv, use `python -m pip install --editable . --group dev` (pip 25.1+).

For performance measurements, run
`uv run --locked python -m benchmarks.benchmark_factorization --enumerate`.
The [benchmark guide](https://github.com/behackl/zero-sum-sequences/blob/main/benchmarks/README.md)
describes the corpus and its mathematical expectations.
