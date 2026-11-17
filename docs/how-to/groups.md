# Use your own group

**Goal:** compute over a group given by your own elements and addition,
or by a SageMath parent.

## Elements and addition

A {py:class}`~zero_sum_sequences.FiniteAdditiveGroup` needs the elements, the
zero, and the addition. The Klein four-group on $\{0,1,2,3\}$, with
addition by bitwise XOR:

```{doctest}
>>> from operator import xor
>>> from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup
>>> klein = AdditiveSequenceSpace(FiniteAdditiveGroup(range(4), zero=0, add=xor))
>>> list(klein.enumerate_atom_catalogue())
[(1)^2, (2)^2, (3)^2, 1 · 2 · 3]
```

The group axioms are not checked; the addition must define a finite
abelian group.

## Automorphisms

Orbit and symmetry methods need automorphisms. Supply generators as
functions on elements; the full group is generated from them. For
$C_9$, multiplication by 2 generates all automorphisms:

```{doctest}
>>> C9 = FiniteAdditiveGroup(
...     range(9), zero=0, add=lambda a, b: (a + b) % 9,
...     automorphism_generators=[lambda x: 2 * x % 9])
>>> S = AdditiveSequenceSpace(C9)
>>> S.automorphism_group()
AutomorphismGroup(order=6)
>>> len(S(1, 1, 7).orbit())
6
```

Groups from {py:meth}`~zero_sum_sequences.FiniteAdditiveGroup.cyclic_product`
come with their full automorphism group, and are usually the simplest
choice.

## SageMath parents

Finite additive parents from Sage, such as `Zmod(n)`, `GF(p)**k` or
`AdditiveAbelianGroup`, can be passed directly:

```python
from sage.all import GF
from zero_sum_sequences import AdditiveSequenceSpace

S = AdditiveSequenceSpace(GF(3)**2)
```

Sage is optional; see [Getting started](../getting-started.md#sagemath).
