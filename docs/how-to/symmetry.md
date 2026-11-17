# Classify sequences up to automorphism

**Goal:** sort sequences into automorphism classes and find the
automorphism relating two of them.

## Test and witness equivalence

```{doctest}
>>> from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup
>>> S = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(2, 4))
>>> A = S((0, 1))**3 * S((1, 0), (1, 1))
>>> B = S((0, 1), (1, 0)) * S((1, 1))**3
>>> A.is_in_same_orbit(B)
True
>>> group = S.automorphism_group()
>>> group.transporter(A, B)
Automorphism(generator_images=((1, 0), (1, 1)))
```

An automorphism is shown by the images of the generators $(1,0)$ and
$(0,1)$; this one maps $(0,1)\mapsto(1,1)$.

## Classify a list

`orbit_representatives` groups sequences by their canonical
representative, the smallest sequence in the orbit:

```{doctest}
>>> sequences = [S((0, 1), (0, 3)), S((1, 1), (1, 3)), S((0, 2))**2, S((1, 0))**2]
>>> for representative, members in group.orbit_representatives(sequences).items():
...     print(representative, "<-", members)
(0, 1) · (0, 3) <- [(0, 1) · (0, 3), (1, 1) · (1, 3)]
(0, 2)^2 <- [(0, 2)^2]
(1, 0)^2 <- [(1, 0)^2]
```

## Classify all atoms

For atoms, the catalogue precomputes the classes:

```{doctest}
>>> atoms = S.enumerate_atom_catalogue()
>>> for orbit in atoms.orbits()[-3:]:
...     print(orbit.representative, len(orbit.indices))
(0, 1)^2 · (1, 1)^2 4
(0, 1) · (0, 2) · (1, 0) · (1, 1) 4
(0, 1)^3 · (1, 0) · (1, 1) 8
```

## Notes

- Groups from `cyclic_product` know their automorphisms. For other groups,
  see [Use your own group](groups.md).
- A canonical form costs one pass over the orbit. For large groups, see
  [Scope and performance](../performance.md).
- See also {py:class}`~zero_sum_sequences.AutomorphismGroup` and
  {py:class}`~zero_sum_sequences.CanonicalForm`.
