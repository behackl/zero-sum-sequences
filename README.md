# Zero-sum sequences

`zero-sum-sequences` computes with finite multisets over additive groups:
zero-sum atoms, factorizations and their lengths, automorphism orbits, and
zero-sum constants. It runs in ordinary Python; SageMath is optional.

## Installation

Requires Python 3.12 or newer:

```console
python -m pip install zero-sum-sequences
```

## A small example

Over `C₂ × C₄`, take a five-term atom and its termwise negative.
Their product has factorization lengths `2`, `4`, and `5`—but not `3`:

```python
from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup

G = FiniteAdditiveGroup.cyclic_product(2, 4)
S = AdditiveSequenceSpace(G, davenport_bound=5)
A = S([(0, 1)] * 3 + [(1, 0), (1, 1)])
B = S([(0, 3)] * 3 + [(1, 0), (1, 3)])

A.is_atom(), B.is_atom()        # (True, True)
sorted((A + B).length_set())    # [2, 4, 5]
```

**[Documentation][docs]** — installation options, tutorial, definitions, and
API reference.

[MIT License](https://github.com/behackl/zero-sum-sequences/blob/main/LICENSE).

[docs]: https://zero-sum-sequences.readthedocs.io/
