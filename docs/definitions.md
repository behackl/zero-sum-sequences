# Definitions

This is a short notation guide, not an introduction to zero-sum theory.
For background, see {cite:p}`gao-geroldinger-2006,geroldinger-halter-koch-2006`.

## Sequences and subsequences

A sequence over an additive abelian group $G$ is a finite multiset, written

$$
S = \prod_{g\in G} g^{v_g(S)}.
$$

Its **length** $|S|=\sum_g v_g(S)$ counts terms, its **support** consists
of the terms with positive multiplicity, and its **sum** is
$\sigma(S)=\sum_g v_g(S)g$. A subsequence is obtained by decreasing
multiplicities. The empty sequence has length zero and sum zero.

The implementation uses `S + T` for multiset union, despite the
multiplicative notation customary in the literature.

## Zero-sum sequences and atoms

A sequence is **zero-sum** if its sum is zero. An **atom** is a nonempty
zero-sum sequence with no proper nonempty zero-sum subsequence. In
particular, the empty sequence is not an atom.

The singleton containing the group identity is an atom. The package's
atom catalogues omit that singleton, and its factorization solver accepts
sequences without identity terms. Zero-sum constants include identity terms.

In standard monoid terminology, the full monoid of zero-sum sequences is
already **reduced**: its only unit is the empty sequence.

Restricting the permitted terms to a subset $G_0\subseteq G$ does not
change the ambient addition. The support need not be a subgroup.

## Factorizations and sets of lengths

A factorization of a zero-sum sequence expresses it as a product of
atoms. Its **factorization length** is the number of atom factors, not
the number of terms in the original sequence. The **set of lengths**
$\mathsf L(S)$ contains the lengths of all such factorizations.

For example, over $C_3$,

$$
\begin{gathered}
1^3 2^3 = (1^3)(2^3) = (1\cdot 2)^3,\\
\mathsf L(1^3 2^3)=\{2,3\}.
\end{gathered}
$$

See {py:class}`~zero_sum_sequences.FactorizationRelation` for cancellation
and the factorization distance.

## Zero-sum constants

For a finite abelian group $G$, the **exponent** $\exp(G)$ is the least
common multiple of its element orders.

For a set $L$ of positive integers, $s_L(G)$ is the least $\ell$ such
that every sequence of length $\ell$ over $G$ has a nonempty zero-sum
subsequence with length in $L$. If no such $\ell$ exists, the constant
is infinite. An **$L$-free** sequence has no such subsequence.

Common special cases are

$$
\begin{aligned}
D(G)&=s_{\mathbb N_{>0}}(G),\\
s(G)&=s_{\{\exp(G)\}}(G),\\
\eta(G)&=s_{\{1,\ldots,\exp(G)\}}(G).
\end{aligned}
$$

Here $D(G)$ is the **Davenport constant**, also the maximum length of an
atom. The notation $s(G)$ is the Erdős–Ginzburg–Ziv constant. The general
API accepts finite nonempty sets of positive allowed lengths; see
{py:func}`~zero_sum_sequences.zero_sum_constant`.

## Symmetry

Group automorphisms act on sequences term by term. They preserve
zero-sumness, atoms, and factorization lengths. A **canonical
representative** selects one element from an orbit.

With symmetry reduction, the extremal sequences returned for a zero-sum
constant are representatives, not necessarily the complete list. The
chosen action and representative convention are documented in the
[API reference](api/index.rst).
