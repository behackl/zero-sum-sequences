# Zero-sum sequences

`zero-sum-sequences` provides immutable finite additive sequences and tools for
enumerating their factorizations into minimal zero-sum sequences. It is built
for SageMath and works with compatible additive Sage parents.

The package keeps group-specific mathematics explicit. Callers configure the
ambient parent and provide an upper bound for its Davenport constant; the
package does not infer structural invariants from the parent.

## Installation

SageMath 10.8 is the supported runtime. Install a tagged source release into
Sage's Python environment with:

```console
sage -pip install \
  "zero-sum-sequences @ git+https://github.com/behackl/zero-sum-sequences.git@v0.1.0"
```

For development, clone the repository and install it with its test tools:

```console
sage -pip install -e '.[dev]'
sage -python -m pytest
```

There is no separate ordinary-CPython execution mode.

## Tutorial

[![Launch the tutorial on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/behackl/zero-sum-sequences/main?urlpath=lab/tree/notebooks/tutorial.ipynb)

The executable [tutorial](notebooks/tutorial.ipynb) introduces the public API
with small hand-checkable examples. It uses ordinary Python syntax rather than
Sage's preparser syntax.

## Additive sequences

Configure an ambient Sage parent and a Davenport upper bound once, then use the
resulting callable space to construct sequences:

```python
from sage.all import Zmod
from zero_sum_sequences import AdditiveSequenceSpace

group = Zmod(3)
Sequences = AdditiveSequenceSpace(group, davenport_bound=3)

sequence = Sequences([1, 1, 2, 2])
sequence.is_zero_sum()  # True
sequence.is_atom()      # False: 1 * 2 is a proper zero-sum subsequence
sequence.multiplicities
```

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
Sage `DiGraph` whose vertices are remainder sequences and whose edge labels
are the removed atoms.

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
sage -python -m benchmarks.benchmark_factorization --enumerate
```

The cases and their mathematical expectations are documented in
[`benchmarks/README.md`](benchmarks/README.md).

## License

The source code is available under the MIT License.
