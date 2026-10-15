# Getting started

## Installation

The package requires Python 3.12 or newer:

```console
python -m pip install zero-sum-sequences
```

For the development version in this checkout:

```console
python -m pip install --editable .
```

The examples below describe the development version on the
`zero-sum-constants` branch; its new invariant functions are not yet part of
the published 0.2.0 release.

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
reduced-atom convention.

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

An explicit `davenport_bound=...` is also accepted and remains authoritative.
An insufficient supplied bound or incomplete supplied atom catalogue can
make atom and factorization results incomplete. Automatic bound resolution
does not launch an exhaustive exact search.

For the precise bound policy, see
{py:attr}`zero_sum_sequences.AdditiveSequenceSpace.davenport_bound`;
for the invariant search, see {py:func}`zero_sum_sequences.zero_sum_constant`.
