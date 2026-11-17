# Scope and performance

The library is for explicit computations in small groups: examples,
counterexamples, and checking conjectures on small cases. Every result is
exact, and every search is exhaustive. This page explains when results are
complete and what is feasible.

## When results are complete

- **Zero-sum constants** are always exact; they need no bounds.
- **Atom catalogues and factorizations** search atoms up to an upper bound
  for the Davenport constant $D(G)$. The bound the space derives itself is
  always sufficient. A bound you supply with `davenport_bound=` is used as
  given; if it is smaller than $D(G)$, atoms can be missed.
- **Catalogues you pass** to factorization queries must contain every atom
  dividing the sequence, as the full catalogue and catalogues of a support
  $G_0$ for sequences over $G_0$ do.

If the available atoms provably cannot factor a zero-sum sequence,
queries raise an error rather than answer. Missing lengths in a length
set cannot be detected this way.

## What is feasible

Single runs of version 0.3.1 on one core of a 2024 laptop (Apple M4).
Other machines differ; a desktop with an AMD Ryzen 7 7700X took 10–30%
longer. Memory is given where it exceeds a few hundred megabytes.

| Computation | Time | Memory |
|---|---:|---:|
| all atoms of $C_5^2$ (31 028) or $C_3^3$ (27 911) | 0.3 s | |
| all atoms of $C_6^2$ (390 860) | 6 s | 0.4 GB |
| all atoms of $C_7^2$ (9.5 million) | 4–5 min | 8 GB |
| length set of a product of two atoms, $C_3^3$, with catalogue | < 1 ms | |
| length set of $(AB)^6$, $\lvert AB\rvert=14$, $C_3^3$ | 0.6 s | |
| $\eta(C_3^3)$, $\eta(C_7^2)$ | 1 s, 7 s | |
| $\mathsf s(C_3^3)$, $\mathsf s(C_6^2)$ | 5 s, 2.5 min | 1 GB |
| $\eta$ or $\mathsf s$ of $C_2^5$ | > 10 min | |

The atoms are usually the limit: their number grows rapidly with the
group, and a catalogue must fit in memory. A
[restricted support](how-to/support.md) avoids the full catalogue.
Factorizations get expensive with the length and multiplicities of the
sequence, not with the group. Zero-sum constants search all sequences
avoiding the condition, up to symmetry, and grow quickly with the order.

## Symmetry

Symmetry saves work when one representative stands for many sequences,
but finding the representative of a sequence has a cost of its own. Over
$C_3^3$, with 11 232 automorphisms, it takes about 5 ms per sequence,
while the length set of a product of two atoms takes 0.2 ms. Hence:

- **For a family of sequences**, iterate over orbit representatives once,
  as in [Length sets for a family](how-to/families.md). Over $C_3^3$, the
  13 104 products $A\cdot(-A)$ of atoms of maximal length come down to
  three solves.
- **For single long sequences**, `symmetry=True` in a
  {py:class}`~zero_sum_sequences.FactorizationCache` pays off: the length
  set of $(AB)^3$ over $C_3^3$ already takes as long as a representative.
- **For many short sequences in no particular order**, leave symmetry
  off. Over $C_3^3$, the 13 104 products above take 2 s without symmetry
  and 1.5–2 min with full canonical forms;
  {py:class}`~zero_sum_sequences.AnchoredCanonicalizer` is a cheaper
  compromise at 3 s.

Zero-sum constants use symmetry by default; it changes the cost, not the
result.

## Reproducing timings

The repository's
[benchmark suite](https://github.com/behackl/zero-sum-sequences/blob/main/benchmarks/README.md)
measures related workloads and records reference timings.
