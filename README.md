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

Over the cyclic group `C₃`, the six-term sequence `1³ 2³` has
factorizations into either two or three atoms:

```python
from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup

G = FiniteAdditiveGroup.cyclic_product(3)
S = AdditiveSequenceSpace(G, davenport_bound=3)
sequence = S([(1,)] * 3 + [(2,)] * 3)

sequence.is_zero_sum()          # True
sequence.is_atom()              # False
sorted(sequence.length_set())   # [2, 3]
```

**[Documentation sources][docs]** — installation options, tutorial, definitions,
and API reference.

[MIT License](https://github.com/behackl/zero-sum-sequences/blob/main/LICENSE).

[docs]: https://github.com/behackl/zero-sum-sequences/blob/docs/shibuya/docs/index.md
