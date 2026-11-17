# Find factorizations with given properties

**Goal:** decide whether a factorization with a prescribed length, and
possibly a prescribed kind of factors, exists, and find one.

Over $C_3\oplus C_3$, take the atom $A=(0,1)^2(1,0)^2(1,1)$ and
$S=A^2(-A)^2$, of length 20:

```{doctest}
>>> from zero_sum_sequences import (
...     AdditiveSequenceSpace, FactorizationSolver, FiniteAdditiveGroup)
>>> G = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3, 3))
>>> A = G((0, 1))**2 * G((1, 0))**2 * G((1, 1))
>>> solver = FactorizationSolver(A**2 * (-A)**2)
```

## A given length

Targeted queries search only as far as needed; they do not compute the
full length set first.

```{doctest}
>>> solver.has_factorization_of_length(3), solver.minimum_factorization_length()
(False, 4)
>>> for atom in solver.factorization_witness(5):
...     print(atom)
(0, 1) · (0, 2)
(0, 1) · (2, 0)^2 · (2, 2)
(0, 1)^2 · (1, 0)^2 · (1, 1)
(0, 2) · (1, 0)^2 · (1, 1)
(0, 2)^2 · (2, 0)^2 · (2, 2)
```

## A given length and kind of factors

`factor_predicate` selects factors; `minimum_matching_factors` requires
at least that many of them. Here: six factors, at least two of length two.

```{doctest}
>>> def short(atom):
...     return len(atom) == 2
>>> solver.factorization_witness(
...     6, minimum_matching_factors=2, factor_predicate=short) is not None
True
>>> solver.factorization_witness(
...     6, minimum_matching_factors=3, factor_predicate=short) is None
True
```

So every factorization of length six has at most two factors of length two.
`None` means that no factorization satisfies both conditions together; the
length itself can still occur.

## All factorizations of one length

```{doctest}
>>> sum(1 for _ in solver.factorizations(factor_count=6))
53
```

Enumerating all factorizations grows quickly with the sequence; prefer
witnesses and length sets where they suffice.

See also {py:class}`~zero_sum_sequences.FactorizationSolver` and
{py:class}`~zero_sum_sequences.FactorizationRelation`.
