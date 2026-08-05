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

For development, clone the repository and install it with its test tools:

```console
python -m pip install -e '.[dev]'
python -m pytest
```

The optional `sage` extra installs the `sagelite` runtime on supported
platforms and enables the Sage integration tests:

```console
python -m pip install -e '.[dev,sage]'
```

## Tutorial

[![Launch the tutorial on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/behackl/zero-sum-sequences/main?urlpath=lab/tree/notebooks/tutorial.ipynb)

The executable [tutorial](notebooks/tutorial.ipynb) introduces the public API
with small hand-checkable examples in ordinary Python.

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

The parent must be finite and iterable. Enumeration is combinatorial, and its
completeness depends on the configured Davenport bound being valid.

The caller is responsible for catalogue completeness. A catalogue used for a
complete result must contain every atom divisor relevant to the input.

Factorizations are computed in the reduced block monoid: identity terms are
not accepted by the solver and are not stored in an `AtomCatalogue`.

## Benchmarks

Run the short-to-very-long performance corpus with:

```console
python -m benchmarks.benchmark_factorization --enumerate
```

The cases and their mathematical expectations are documented in
[`benchmarks/README.md`](benchmarks/README.md).

## License

The source code is available under the MIT License.
