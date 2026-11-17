# Length sets for a family of sequences

**Goal:** the length sets of all products $AB$ of two atoms over
$G=C_3\oplus C_3$, counted with multiplicity, without solving each of the
$68^2$ products separately.

## Enumerate the atoms once

```{doctest}
>>> from collections import Counter
>>> from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup
>>> S = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3, 3))
>>> atoms = S.enumerate_atom_catalogue()
>>> len(atoms)
68
```

## Iterate over orbit representatives

Automorphisms do not change length sets. Choose $A$ up to automorphism,
then $B$ up to the automorphisms fixing $A$: every product $AB$ is then
the image of exactly one solved pair. Orbit sizes recover the counts.

```{doctest}
>>> group = S.automorphism_group()
>>> cache = S.factorization_cache(atoms)
>>> counts = Counter()
>>> pairs = 0
>>> for first in atoms.orbits():
...     A = first.representative
...     for second in atoms.orbits(group=group.stabilizer(A)):
...         B = second.representative
...         lengths = tuple(sorted(cache.length_set(A * B)))
...         counts[lengths] += len(first.indices) * len(second.indices)
...         pairs += 1
>>> pairs
126
>>> sorted(counts.items())
[((2,), 1800), ((2, 3), 2248), ((2, 3, 4), 552), ((2, 3, 4, 5), 24)]
>>> sum(counts.values()) == len(atoms)**2
True
```

126 pairs cover all 4624 ordered pairs; some of their products coincide,
so the cache solves even fewer. Passing the catalogue to the cache means
atoms are looked up rather than searched for each product.

## Save the results

Long computations can be saved and continued later. The cache file is
bound to the catalogue by a digest, so it can only be loaded with the same
catalogue.

```python
atoms.to_jsonl("atoms.jsonl")
cache.to_jsonl("lengths.jsonl")

atoms = AtomCatalogue.from_jsonl("atoms.jsonl", S)
cache = FactorizationCache.from_jsonl("lengths.jsonl", atoms)
```

## Notes

- For sequences that arrive in no particular order, the cache can use
  symmetry itself, with `symmetry=True`. This pays off only when single
  solves are expensive; see [Scope and performance](../performance.md).
- See also {py:meth}`~zero_sum_sequences.AtomCatalogue.orbits`,
  {py:meth}`~zero_sum_sequences.AutomorphismGroup.stabilizer` and
  {py:class}`~zero_sum_sequences.FactorizationCache`.
