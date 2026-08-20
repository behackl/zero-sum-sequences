# Zero-sum sequences

`zero-sum-sequences` provides immutable finite additive sequences and tools for
enumerating their factorizations into minimal zero-sum sequences. The runtime
is ordinary Python with NetworkX; SageMath is supported as an optional source
of additive parents, but is not required.

The package keeps group-specific mathematics explicit. Callers configure the
ambient parent and provide an upper bound for its Davenport constant; the
package does not infer structural invariants from the parent.

## Installation

Install a tagged source release into Python 3.12 or newer with:

```console
python -m pip install \
  "zero-sum-sequences @ git+https://github.com/behackl/zero-sum-sequences.git@v0.1.0"
```

For a reproducible development environment using the committed `uv.lock`,
clone the repository and run:

```console
uv sync --extra dev
uv run python -m pytest
```

Alternatively, install the package and its test tools with `pip`:

```console
python -m pip install -e '.[dev]'
python -m pytest
```

The optional `sage` extra installs the `sagelite` runtime on supported
platforms and enables the Sage integration tests:

```console
python -m pip install -e '.[dev,sage]'
```

With `uv`, use `uv sync --extra dev --extra sage` instead.

## Tutorial

[![Launch the tutorial on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/behackl/zero-sum-sequences/main?urlpath=lab/tree/notebooks/tutorial.ipynb)

The executable [tutorial](https://github.com/behackl/zero-sum-sequences/blob/main/notebooks/tutorial.ipynb) introduces the public API
with small hand-checkable examples in ordinary Python.
Its execution dependencies are available through
`uv sync --extra tutorial` or the corresponding `pip` extra.

## Additive sequences

Configure an ambient parent and a Davenport upper bound once, then use the
resulting callable space to construct sequences. `FiniteAdditiveGroup` is a
small convenience adapter for groups represented by ordinary Python values:

```python
from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup

group = FiniteAdditiveGroup(
    range(3),
    zero=0,
    add=lambda left, right: (left + right) % 3,
    coerce=lambda value: int(value) % 3,
)
Sequences = AdditiveSequenceSpace(group, davenport_bound=3)

sequence = Sequences([1, 1, 2, 2])
sequence.is_zero_sum()  # True
sequence.is_atom()      # False: 1 * 2 is a proper zero-sum subsequence
sequence.multiplicities
```

Existing algebra systems can be used directly. A compatible parent is
callable for coercion and provides `zero()`; its elements must be hashable and
mutually orderable. Elements may implement `+`, or the parent may provide
`add(left, right)`. Exhaustive catalogue enumeration additionally requires the
parent to be finite, iterable, and to provide `is_finite()`.

A sequence is immutable and hashable. Addition combines multisets,
subtraction removes a subsequence, and multiplication by a non-negative
integer repeats a sequence. Arithmetic preserves the sequence space, and the
empty sequence retains the base parent and its zero element.
`sequence.map_terms(mapping)` applies a map to every term and reconstructs the
result as a canonical multiset; pass `target_space=` when the image belongs to
a different sequence space.

The configured bound must not be smaller than the actual Davenport constant
when complete atom or factorization results are required.

## Factorizations

The factorization engine indexes relevant atom divisors as sparse
multiplicity vectors and merges equal remainders in a directed acyclic graph.
Attained lengths are represented internally as integer bitsets.

```python
lengths = sequence.length_set()
witnesses = sequence.factorization_witnesses()
factorizations = list(sequence.factorizations())
graph = sequence.factorization_digraph()
```

`factorization_witnesses()` retains one factorization for every attained
length. Exhaustive `factorizations()` is necessarily output-sensitive, but it
emits each unordered factorization once. `factorization_digraph()` returns a
NetworkX `DiGraph` whose vertices are remainder sequences and whose edges
store the removed atom in their `"atom"` attribute.

For several queries against the same remainder DAG, use the public solver:

```python
from zero_sum_sequences import FactorizationSolver

solver = FactorizationSolver(sequence)
solver.length_set()
solver.factorization_witnesses()
solver.statistics
solver.digraph()
```

A complete precomputed catalogue can avoid rediscovering atoms:

```python
from zero_sum_sequences import AtomCatalogue, FactorizationSolver

catalogue = AtomCatalogue(Sequences, atoms)
solver = FactorizationSolver(sequence, atom_catalogue=catalogue)
```

For a small finite parent, a complete reduced catalogue can instead be
generated exhaustively through the configured Davenport bound:

```python
catalogue = Sequences.enumerate_atom_catalogue()
```

The parent must be a finite iterable additive group. Enumeration completes
each sorted prefix with its uniquely determined final term, rather than
testing multisets whose sum is nonzero. Completeness depends on the configured
Davenport bound being valid.

The caller is responsible for catalogue completeness. A catalogue used for a
complete result must contain every atom divisor relevant to the input.

Factorizations are computed in the reduced block monoid: identity terms are
not accepted by the solver and are not stored in an `AtomCatalogue`.

## Automorphism orbits

Group automorphisms act on sequences term by term. Because an
`AdditiveSequence` is a multiset, the induced action automatically disregards
the order of its terms. The package can materialize an orbit, test orbit
membership, and return a shortest deterministic word in the configured
automorphism generators.

The convenience constructor for a product of cyclic groups configures the
coordinate group and generators of its full automorphism group together. For
example, the eight maximal-length atoms over $C_2\oplus C_4$ form one orbit:

```python
G = FiniteAdditiveGroup.cyclic_product(2, 4)
C2xC4 = AdditiveSequenceSpace(G, davenport_bound=5)

atom = C2xC4([(0, 1)] * 3 + [(1, 0), (1, 1)])
other = C2xC4([(0, 1), (1, 0)] + [(1, 1)] * 3)

len(atom.orbit())                         # 8
atom.is_in_same_orbit(other)              # True
witness = atom.orbit_witness(other)
witness.show()
# (1, 0) ↦ (1, 0)
# (0, 1) ↦ (1, 1)
```

Automorphism data is resolved only on the first orbit query and then cached on
the sequence space. A returned witness retains this context, so `show()` can
display the induced homomorphism on the base parent's distinguished additive
generators. `FiniteAdditiveGroup.cyclic_product(...)` supplies both the
standard additive generators and elementary coordinate scalings and shears.
Finite-dimensional Sage vector spaces over finite fields are recognized
automatically and use their basis and generators of their general linear
group. A custom `FiniteAdditiveGroup` can receive `additive_generators=` and
callable `automorphism_generators=` at construction; callers can also pass an
`AutomorphismAction` explicitly through the `action=` keyword.

Orbit traversal is breadth-first and therefore requires a finite orbit.
Automatic discovery raises `AutomorphismActionUnavailable` when the parent
does not expose suitable generators, in which case an explicit action is
required.

## Benchmarks

Run the short-to-very-long performance corpus with:

```console
python -m benchmarks.benchmark_factorization --enumerate
```

The cases and their mathematical expectations are documented in
[`benchmarks/README.md`](benchmarks/README.md).

## License

The source code is available under the MIT License.
