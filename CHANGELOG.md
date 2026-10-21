# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[semantic versioning](https://semver.org/). While the version is below 1.0,
minor releases may contain breaking changes; they are listed first.

## [Unreleased]

## [0.3.1] - 2026-10-21

A performance release with no API changes. Results are unchanged, except
that symmetry-reduced factorization caches may return different (equally
valid) witnesses.

### Changed

- Faster sequence comparisons, which no longer expand the term lists, and
  faster catalogue divisor queries.
- Atom catalogue enumeration extends only zero-sum free prefixes and needs no
  separate atom tests; catalogues are tens to thousands of times faster, for
  instance under a second instead of about half a minute for C₃³.
  `is_atom()` checks a single set of subsequence sums instead of enumerating
  subsequences.
- Automorphism groups of image-table automorphisms, such as those built by
  `AutomorphismGroup.closure()`, compute images, orbits, orbit minima,
  stabilizers and transporters with integer keys instead of building image
  sequences. Small orbits are traversed through the generators, larger ones
  found by scanning the group. Canonical forms of generic sequences are 25
  to 70 times faster, stabilizers and orbit classification 3 to 9 times.
  `FactorizationCache` with symmetry benefits accordingly; for such groups it
  transports witnesses with the smallest automorphism realizing the
  canonical form, so returned witnesses may differ from 0.3.0.
- The factorization solver packs remainder multiplicities into single integers
  with guard bits, so testing whether an atom fits and computing the remainder
  take one subtraction. Together with the faster atom tests, length sets are
  about four to five times faster.
- Zero-sum constant searches with symmetry reduction find canonical forms
  without scanning the whole group at every step; η(C₃³) takes under a second
  instead of about four.
- `sum()` and `subsequences()` no longer coerce terms that are already
  canonical.

### Added

- `benchmarks/benchmark_operations.py` times catalogue enumeration, atom
  tests, length sets, automorphism operations, and zero-sum constants over
  several kinds of groups.

## [0.3.0] - 2026-10-20

Sequences now use the multiplicative notation of the literature, and the
documentation follows its terminology. Breaking changes take effect without a
deprecation period.

### Breaking changes

- Sequences are written multiplicatively. `S(g, h, ...)` constructs the
  sequence `g·h·…`, `S * T` is the product, `S ** k` a power, `S / T` is
  `S·T⁻¹` for a divisor `T` of `S`, and `-S` negates every term.
  `S.from_terms(iterable)` builds a sequence from an iterable, and
  `math.prod()` multiplies sequences. Removed: `S([...])`, `+`, `-`, integer
  `*`, and `sum()` support.
- `AdditiveSequence.total()` is renamed to `sum()`.
- Symmetry keywords follow one rule. Queries whose answer depends on the group
  take `group=`, which defaults to the space's automorphism group: the
  catalogue orbit methods, `AnchoredCanonicalizer` (whose group is now
  optional), and the sequence orbit methods, where `group=` replaces `action=`
  and also accepts an `AutomorphismGroup`. Queries that use symmetry only to
  save work take `symmetry=True|False|AutomorphismGroup`: the zero-sum
  constants (default `True`) and `FactorizationCache`,
  `AdditiveSequenceSpace.factorization_cache()` and
  `FactorizationCache.from_jsonl()` (default `False`, replacing `group=`).
- `AdditiveSequence.orbit()` returns the orbit sorted, not in breadth-first
  order.
- Factorization queries raise `ValueError` when an explicit Davenport bound or
  a supplied atom catalogue provably cannot factor a zero-sum sequence, instead
  of returning empty results. Queries for a single length are unchanged.
- `FactorizationCache.statistics` is a property returning a
  `FactorizationCacheStatistics` dataclass instead of a method returning a
  dict. `clear()` now also resets the counters.
- `AutomorphismGroup` finds the identity among `Automorphism` elements or, when
  `compose` is given, as the idempotent element; otherwise `identity=` is
  required and its absence raises `TypeError`.
- Short `repr`s: sequences display as their product, groups and spaces by
  their construction, e.g. `FiniteAdditiveGroup.cyclic_product(2, 4)`.
- Length caps and exponents reject booleans and other non-integers.
- The tutorial is a marimo notebook (`notebooks/tutorial.py`);
  `notebooks/tutorial.ipynb` and the Binder setup are removed, and the
  `tutorial` extra installs marimo and Matplotlib. Development tools moved
  from the `dev` extra to a `dev` dependency group
  (`pip install --group dev`).

### Added

- Exact zero-sum constants: `zero_sum_constant` (with extremal sequences),
  `davenport_constant`, `erdos_ginzburg_ziv_constant`, `eta_constant`, and
  `exponent`, with symmetry reduction.
- The Davenport bound of a sequence space is optional. Without one, a
  certified bound is derived without an exhaustive search, using exact
  formulas where they apply; exact values are cached and shared between
  spaces over the same group.
- `FiniteAdditiveGroup.order()`, `element_order()`, `exponent()`, and
  `negate()`.
- Groups, spaces, sequences, and catalogues can be pickled; sequences sent to
  worker processes compare equal to the caller's.
- Orbit methods use an automorphism group supplied by the parent when it has
  no explicit generators.
- `FactorizationCache.from_jsonl()` accepts `canonicalize=` and loads files
  saved with any symmetry setting.
- `FactorizationStatistics` and `FactorizationCacheStatistics` are exported.
- Documentation at <https://zero-sum-sequences.readthedocs.io/>, with
  definitions, an API reference, and the tutorial, runnable in the browser.
- Python 3.15 is tested.

### Changed

- Maximum-only factorization queries use less memory.
- NetworkX is imported only when a factorization graph is requested.
- Sage is probed only for Sage parents, and unavailable symmetry is resolved
  once per space.

### Fixed

- Parents whose `zero()` returns a non-canonical label give correct sums and
  zero-sum tests.

## [0.2.0] - 2026-09-17

### Added

- **Materialized automorphism groups.** For a finite parent, the whole
  automorphism group can now be built explicitly
  (`space.automorphism_group()`), giving orbits, stabilizers, canonical orbit
  representatives (smallest element + witness automorphism), transporters,
  and orbit classification of atom catalogues in one pass. `AdditiveSequence`
  gained a total order (length, then terms) so canonical forms are
  well-defined regardless of how the group was found. Cross-checked against
  brute force for C_42, C_2×C_4, C_2×C_6, C_4×C_4, C_3×C_9, C_2^3, C_2^4,
  C_3^3.
- **Indexed, serializable atom catalogues.** `AtomCatalogue` atoms now have
  stable indices (`catalogue.index(atom)`, `catalogue[i]`, `in`), can be
  filtered/restricted (`filter`, `restrict(max_length=..., support=...)`),
  classified into automorphism orbits (`AtomOrbit`), and round-tripped
  through `to_jsonl`/`from_jsonl` with digest-based re-verification on load.
  Parents gained `encode_term`/`decode_term`; spaces gained `decode()`,
  sequences gained `encode()`.
- **`FactorizationCache`**: a memoized, optionally group-aware front end for
  length-set queries (`length_set`, `minimum`, `maximum`, `has_length`,
  `witness`, `witnesses`, `factorizations`) over a catalogue, with an
  unbounded length-set cache and LRU-bounded solvers. Supports
  `preload`/`to_jsonl`/`from_jsonl` scoped to the catalogue digest and group
  order. Ships alongside `AnchoredCanonicalizer`, which transports the
  largest atom divisor to its orbit representative and reduces under its
  stabilizer.
- New public API: `AutomorphismGroup`, `Automorphism`, `CanonicalForm`,
  `AtomOrbit`, `AnchoredCanonicalizer`, `FactorizationCache`.

### Breaking changes

- `FactorizationOracle` is renamed to `FactorizationCache` (`space.oracle()`
  → `space.factorization_cache()`).

## [0.1.2] - 2026-09-05

### Added

- Targeted extraction of individual factorization witnesses, including
  profile constraints, without computing the full length set.
- Atom catalogue generation on prescribed supports via
  `enumerate_atom_catalogue(support=...)`.

### Changed

- Faster exact-length and minimum/maximum factorization queries, with
  improved cache reuse.
- Expanded documentation, tutorial examples, and regression tests.

## [0.1.1] - 2026-09-04

### Added

- Fixed-length factorization queries and enumeration via
  `has_factorization_of_length()` and `factorizations(factor_count=...)`.
- Immutable `FactorizationRelation` objects with canonicalization,
  common-factor cancellation, reduced relations, and standard factorization
  distance.

### Changed

- Expanded tests and tutorial documentation for the new APIs.
- Removed runtime assertions from the factorization engine.

## [0.1.0] - 2026-08-20

First release as a standalone package: sequences over configurable additive
parents, factorization queries, exhaustive atom catalogue enumeration,
automorphism orbits, and optional SageMath integration.

[Unreleased]: https://github.com/behackl/zero-sum-sequences/compare/v0.3.1...HEAD
[0.3.1]: https://github.com/behackl/zero-sum-sequences/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/behackl/zero-sum-sequences/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/behackl/zero-sum-sequences/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/behackl/zero-sum-sequences/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/behackl/zero-sum-sequences/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/behackl/zero-sum-sequences/releases/tag/v0.1.0
