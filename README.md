# Zero-sum sequences

`zero-sum-sequences` provides immutable finite additive sequences and tools for
enumerating their factorizations into minimal zero-sum sequences. The runtime
is ordinary Python with NetworkX; SageMath is supported as an optional source
of additive parents, but is not required.

Callers configure the ambient additive parent. A certified Davenport upper
bound is derived when needed, or can be supplied explicitly; deriving it
never launches an exhaustive search.

## Installation

Install the latest release from PyPI into Python 3.12 or newer:

```console
python -m pip install zero-sum-sequences
```

The optional `sage` extra installs the `sagelite` runtime on supported
platforms:

```console
python -m pip install "zero-sum-sequences[sage]"
```

For a reproducible development environment using the committed `uv.lock`,
clone the repository and run:

```console
uv sync
uv run python -m pytest
```

The `dev` dependency group, which `uv sync` installs by default, contains the
test tools, the tutorial dependencies, and JupyterLab. Alternatively, install it with `pip`
(version 25.1 or newer):

```console
python -m pip install -e . --group dev
python -m pytest
```

To include the Sage integration tests, use `uv sync --extra sage` or
`python -m pip install -e '.[sage]' --group dev`.

## Tutorial

[![Launch the tutorial on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/behackl/zero-sum-sequences/main?urlpath=lab/tree/notebooks/tutorial.ipynb)

The executable [tutorial](https://github.com/behackl/zero-sum-sequences/blob/main/notebooks/tutorial.ipynb) introduces the public API
with small hand-checkable examples in ordinary Python.
Its execution dependencies are included in the development environment above,
and are also available on their own as the `tutorial` extra
(`python -m pip install "zero-sum-sequences[tutorial]"`). In the development
environment, open it with `uv run jupyter lab`.

## Additive sequences

Configure an ambient parent, then use the resulting callable space to
construct sequences. `FiniteAdditiveGroup` is a small convenience adapter
for groups represented by ordinary Python values:

```python
from zero_sum_sequences import AdditiveSequenceSpace, FiniteAdditiveGroup

group = FiniteAdditiveGroup(
    range(3),
    zero=0,
    add=lambda left, right: (left + right) % 3,
    coerce=lambda value: int(value) % 3,
)
Sequences = AdditiveSequenceSpace(group)

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

An explicit Davenport bound is used as supplied, without being checked:

```python
Sequences = AdditiveSequenceSpace(group, davenport_bound=3)
```

Without one, `Sequences.davenport_bound` derives and caches a certified upper
bound on first use: the exact value when known (computed earlier, or the
formula for cyclic groups, groups of rank at most two and p-groups with
known invariant factors), else `min(n, e * (1 + ceil(log2(n / e))))` from
the order `n` and a cheaply known exponent `e`, else `n`. (For other groups
`1 + sum(n_i - 1)` is only a lower bound.) Exact results tighten a derived
bound but never replace an explicit one. If no bound can be derived,
`davenport_bound` raises `ValueError`.

Complete atom and factorization results need a bound covering the relevant
atom lengths; a valid bound for `D(G)` suffices. Without any bound,
factorizing a sequence searches all its distinct subsequences, which can be
expensive. A complete atom catalogue removes the need for a bound.

### Group arithmetic

The dependency-free adapter exposes cardinality, additive orders, exponent
and inverses directly:

```python
G = FiniteAdditiveGroup.cyclic_product(2, 4)
G.order()                      # 8
G.element_order((1, 1))         # 4
G.exponent()                   # 4
G.negate((1, 1))               # (1, 3)
```

This data is computed lazily and cached, and spaces over the same
`FiniteAdditiveGroup` share it. Other parents are queried through their
native methods (`cardinality()`, `exponent()`, `invariants()`,
`additive_order()`, ...) when available; `exponent(Sequences)` works for any
supported finite parent.

## Factorizations

The factorization engine indexes relevant atom divisors as sparse
multiplicity vectors and merges equal remainders in a directed acyclic graph.
Attained lengths are represented internally as integer bitsets.

```python
lengths = sequence.length_set()
minimum = sequence.minimum_factorization_length()
maximum = sequence.maximum_factorization_length()
witnesses = sequence.factorization_witnesses()
factorizations = list(sequence.factorizations())
length_three = list(sequence.factorizations(factor_count=3))
has_length_three = sequence.has_factorization_of_length(3)
graph = sequence.factorization_digraph()
```

`factorization_witnesses()` retains one factorization for every attained
length. Exhaustive `factorizations()` is necessarily output-sensitive, but it
emits each unordered factorization once. Pass `factor_count=` to enumerate
only factorizations of one length; `has_factorization_of_length()` tests that
length with a targeted search unless complete length data is already cached.
`minimum_factorization_length()` uses a shortest-path search, while
`maximum_factorization_length()` solves the length bitsets on the remainder
DAG (so a later `length_set()` is free). Both return `None` when no
factorization exists.
`factorization_digraph()` returns a NetworkX `DiGraph` whose vertices are
remainder sequences and whose edges store the removed atom in their `"atom"`
attribute.

For several queries against the same remainder DAG, use the public solver:

```python
from zero_sum_sequences import FactorizationSolver

solver = FactorizationSolver(sequence)
solver.length_set()
solver.minimum_factorization_length()
solver.maximum_factorization_length()
solver.has_factorization_of_length(3)
solver.factorization_witness(3)
list(solver.factorizations(factor_count=3))
solver.factorization_witnesses()
solver.statistics
solver.digraph()
```

To find just one factorization, use `solver.factorization_witness(k)`. It
returns a deterministic tuple of atoms, or `None` if the requested length is
absent, without first constructing the full remainder DAG or computing the
complete length set. Existing fixed-length answers and complete length data
are reused when available. An optional predicate can constrain the witness:

```python
solver.factorization_witness(
    3,
    minimum_matching_factors=2,
    factor_predicate=lambda atom: len(atom) == 2,
)
```

This asks for three factors, at least two of which have term length two.
`None` means that no witness satisfies both constraints; it does not imply
that length three itself is absent. Profile-constrained searches are also
targeted, but may need to explore more alternatives. For the empty sequence,
the length-zero witness is `()` rather than `None`.

Inspecting `solver.statistics` or calling `solver.digraph()` materializes the
full remainder graph, even if preceding queries were targeted.

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
each sorted prefix with its uniquely determined final term; completeness
depends on the Davenport bound being valid.

To enumerate only atoms supported on a subset of the ambient group, pass
`support=`:

```python
# Here Sequences is the C_3 space from the first example.
restricted = Sequences.enumerate_atom_catalogue(support=[1])
list(restricted)  # Only the atom 1^3; 2 need not belong to the support.
Sequences([1] * 6).length_set(atom_catalogue=restricted)  # {2}
```

The support can be any iterable of terms, not necessarily a subgroup or
closed under negation; duplicates and zero are ignored. Sums and inverses
are computed in the ambient group. The result is complete for sequences
supported there, not for sequences using other terms.

In general the caller is responsible for catalogue completeness: a complete
result needs every atom divisor of the input in the catalogue.

Factorizations use only nonzero group elements: identity terms are not
accepted by the solver or stored in an `AtomCatalogue`. The singleton zero
is still an atom in the full block monoid.

## Factorization relations

A `FactorizationRelation` represents an oriented equality between two
unordered factorizations. It canonicalizes the factors on each side, verifies
that all factors use one sequence space and have equal products, and provides
common-factor cancellation and the standard factorization distance.

For example, over $C_3$ let $A=1^3$, $B=2^3$, and $P=1\cdot2$. Since
$AB=P^3$, adjoining one common factor gives a relation between $ABP$ and
$P^4$:

```python
from zero_sum_sequences import FactorizationRelation

A = Sequences([1, 1, 1])
B = Sequences([2, 2, 2])
P = Sequences([1, 2])

relation = FactorizationRelation(
    source=(A, B, P),
    target=(P, P, P, P),
)
relation.common_factor       # (P,)
relation.reduced_source      # (A, B)
relation.reduced_target      # (P, P, P)
relation.distance            # 3

reduced = relation.reduced()
reduced.is_reduced           # True
reduced.product == A + B     # True
```

Relations are immutable and keep the orientation of source and target.
Factors are assumed to be atoms (as solver results and catalogue entries
are) and are not re-tested. For a relation with two empty sides, pass
`space=Sequences`.

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

## Materialized automorphism groups

For a finite parent the whole automorphism group can be materialized, which
gives stabilizers, canonical orbit representatives, transporters and
classification of many objects at once:

```python
G = FiniteAdditiveGroup.cyclic_product(2, 4)
C2xC4 = AdditiveSequenceSpace(G, davenport_bound=5)
group = C2xC4.automorphism_group()        # AutomorphismGroup(order=8)

x = C2xC4([(1, 1), (1, 1), (0, 2)])
group.orbit(x)                            # the four images, sorted
len(group.stabilizer(x))                  # 2
form = group.canonical_form(x)            # orbit minimum + smallest automorphism reaching it
form.automorphism.generator_images()      # ((1, 0), (1, 1)): the "matrix", column by column
group.transporter(x, group.apply(form.automorphism, x))

atoms = tuple(C2xC4.enumerate_atom_catalogue())
group.orbit_representatives(atoms)        # 11 orbits: representative -> members
```

The group acts on sequences and on tuples of sequences; tuples are treated as
multisets unless `ordered=True` (the canonical form of a tuple also records
the permutation that sorts the image). Canonical forms are minima with
respect to the total order on sequences (length first, then the sorted term
list, exposed through `<`), so representatives do not depend on how the group
was found.

The parent decides how the group is obtained: a `FiniteAdditiveGroup` may
receive a callable `automorphism_group=` returning an `AutomorphismGroup`
built from any element type with `apply_term`, `compose` and `inverse`
callbacks (for example units acting on a cyclic group, or CAS matrices);
otherwise its `automorphism_generators` are closed under composition.
Orbit traversals cost proportionally to the orbit size; stabilizers and
smallest transporters scan the whole group once, so they are cheap for
groups of a few thousand elements and expensive beyond that.
`AutomorphismGroupUnavailable` is raised when neither a provider nor
generators exist.

## Indexed and serialized catalogues

An `AtomCatalogue` stores its atoms in the sequence order (length first,
then terms), so `catalogue.index(atom)` and `catalogue[i]` form a stable
identifier namespace. It can be filtered, classified into automorphism
orbits, and written to or read from JSON lines:

```python
catalogue = C2xC4.enumerate_atom_catalogue()
catalogue.index(atom), catalogue[7], atom in catalogue
short = catalogue.restrict(max_length=3)            # also restrict(support=...)
for orbit in catalogue.orbits():                    # AtomOrbit(representative, indices, stabilizer_order)
    ...
catalogue.representative(atom)                      # orbit minimum of an atom
catalogue.transporter(atom)                         # an automorphism mapping atom to it
catalogue.permutation(group.elements[1])            # an automorphism as an index permutation

catalogue.to_jsonl("atoms.jsonl", annotate=lambda atom: {"label": ...})
same = AtomCatalogue.from_jsonl("atoms.jsonl", C2xC4)
same.annotations[atom]["label"]
catalogue.digest()                                  # sha256 of the atoms, annotation-independent
```

`from_jsonl` checks the order, the indices and the digest, and re-verifies
every record as an atom unless `verify=False`. Terms are serialized through
the parent: `FiniteAdditiveGroup` accepts `encode_term=`/`decode_term=`,
coordinate groups and Sage vector spaces use lists of integers by default,
and `sequence.encode()` / `space.decode(data)` apply the codec to whole
sequences.

## Memoized factorization queries

A `FactorizationCache` answers length-set questions over one catalogue with
caching, and can exploit that length sets are automorphism invariants:

```python
cache = C2xC4.factorization_cache(catalogue)
cache.length_set(x), cache.minimum(x), cache.maximum(x)
cache.has_length(x, 3), cache.witness(x, 3), cache.witnesses(x)
cache.statistics()                                      # hits, misses, solver builds

aware = C2xC4.factorization_cache(catalogue, group=group)
aware.length_set(group.apply(a, x)) == aware.length_set(x)   # one solve per orbit
aware.witness(group.apply(a, x), 3)                     # transported back through a^-1

cache.to_jsonl("lengths.jsonl")                         # bound to catalogue.digest()
FactorizationCache.from_jsonl("lengths.jsonl", catalogue)
```

Length sets are kept without bound; solvers are kept in a bounded LRU
(`maxsize=`) so that follow-up questions about a recent sequence reuse its
remainder graph. With `group=` every query is keyed by a canonical form;
the catalogue must then be closed under the group, which is checked on
construction. The default is the full canonical form, whose orbit traversal is worthwhile
when solving is the expensive part; for large groups pass
`canonicalize=AnchoredCanonicalizer(catalogue, group)`, which only
transports the largest atom dividing the sequence to its orbit
representative and reduces under that representative's stabilizer. Any
callable returning an automorphic image of its argument (optionally with the
automorphism) is accepted as `canonicalize=`.

## Zero-sum constants

For a set `L` of positive integers, `s_L(G)` is the smallest length that
forces a nonempty zero-sum subsequence with length in `L`. The Davenport
constant `D(G) = s_N(G)`, the Erdős–Ginzburg–Ziv constant
`s(G) = s_{exp(G)}(G)` and `η(G) = s_[1, exp(G)](G)` are special cases:

```python
from zero_sum_sequences import (
    davenport_constant,
    erdos_ginzburg_ziv_constant,
    eta_constant,
    exponent,
    zero_sum_constant,
)

C3xC3 = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3, 3))
exponent(C3xC3)                      # 3
davenport_constant(C3xC3)            # 5
eta_constant(C3xC3)                  # 7
erdos_ginzburg_ziv_constant(C3xC3)   # 9

result = zero_sum_constant(C3xC3, {2, 3})
result.value                         # 8
str(result.extremal[0])              # '(0, 0) · (0, 1)^2 · (1, 0)^2 · (1, 1)^2', the only orbit
zero_sum_constant(C3xC3, {1, 2}).value   # math.inf: no multiple of exp(G)
```

`extremal` holds the longest sequences without such a subsequence, one
minimum per orbit of the chosen symmetry group. Because the property is
inherited by subsequences, these are enumerated length by length, tracking
subsequence-sum lengths as bitsets.

Every constant function accepts the same `symmetry=` keyword:

```python
zero_sum_constant(C3xC3, {2, 3}, symmetry=True)   # default: discover automorphisms
zero_sum_constant(C3xC3, {2, 3}, symmetry=False)  # no discovery or orbit reduction
zero_sum_constant(C3xC3, {2, 3}, symmetry=C3xC3.automorphism_group())
```

`True` falls back to no reduction if automorphisms are unavailable. An
explicit `AutomorphismGroup` uses exactly that group, possibly a proper
subgroup. With `False`, `extremal` contains every extremal sequence.
Symmetry affects running time and representatives, never the value.

`lengths` must be a finite nonempty iterable of positive integers. An
infinite `L` is passed as an initial segment: a result that does not exceed
`max(L)` is valid for every set agreeing with `L` up to `max(L)`. The
configured Davenport bound is not used for these constants. The package needs no
computer algebra system; for example `s(C_2^3 ⊕ C_4) = 13` takes about a minute.

## References

Mathematical citations in the API docstrings use the shared
[`docs/references.bib`](docs/references.bib). The
[references page](docs/references.md) distinguishes theorems used by the
package from related computational work.

## Benchmarks

Run the short-to-very-long performance corpus with:

```console
python -m benchmarks.benchmark_factorization --enumerate
```

The cases and their mathematical expectations are documented in
[`benchmarks/README.md`](https://github.com/behackl/zero-sum-sequences/blob/main/benchmarks/README.md).

## License

The source code is available under the MIT License.
