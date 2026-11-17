# Work on a restricted support

**Goal:** factorizations of sequences whose terms lie in a subset
$G_0\subseteq G$, in a group too large for a full atom catalogue.

Over $G=C_7\oplus C_7$ there are about 9.5 million atoms. On the support
$G_0=\{\pm e_1,\pm e_2,\pm(e_1+e_2)\}$ there are 47.

```{doctest}
>>> from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup
>>> S = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(7, 7))
>>> G0 = [(1, 0), (0, 1), (1, 1), (6, 0), (0, 6), (6, 6)]
>>> atoms = S.enumerate_atom_catalogue(support=G0)
>>> atoms
AtomCatalogue(47 atoms)
```

$G_0$ need not be a subgroup: sums are taken in $G$.

## Use the catalogue

Every atom dividing a sequence over $G_0$ is itself a sequence over $G_0$,
so the restricted catalogue is complete for such sequences:

```{doctest}
>>> U = S((1, 0))**3 * S((0, 1))**3 * S((1, 1))**4
>>> U.is_atom()
True
>>> sorted((U * -U).length_set(atom_catalogue=atoms))
[2, 7, 8, 9, 10]
```

## Sequences with other terms

For sequences with terms outside $G_0$, the catalogue is incomplete. When
this makes factorization impossible, queries raise an error instead of
returning a wrong answer:

```{doctest}
>>> V = U * S((2, 3))
>>> (V * -V).length_set(atom_catalogue=atoms)
Traceback (most recent call last):
...
ValueError: the zero-sum sequence ... has no factorization into the available atoms: the atom catalogue lacks some of its atom divisors
```

Missing lengths cannot be detected this way, so only pass catalogues that
contain every atom divisor. Without a catalogue, the atoms dividing each
sequence are found directly.

See also {py:meth}`~zero_sum_sequences.AdditiveSequenceSpace.enumerate_atom_catalogue`
and {py:meth}`~zero_sum_sequences.AtomCatalogue.restrict`.
