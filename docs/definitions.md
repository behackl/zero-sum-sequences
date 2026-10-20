# Definitions

This is a short notation guide, not an introduction to zero-sum theory.
For background, see {cite:p}`gao-geroldinger-2006,geroldinger-halter-koch-2006`.

## Sequences and subsequences

Let $G$ be an additive abelian group. By a **sequence** over $G$ we mean a
finite unordered sequence of terms from $G$, with repetition of terms
allowed. As usual, sequences are elements of the free abelian monoid
$\mathcal F(G)$ with basis $G$ and are written multiplicatively:

$$
S = g_1\cdot\ldots\cdot g_\ell = \prod_{g\in G} g^{\mathsf v_g(S)}.
$$

Here $\mathsf v_g(S)$ is the **multiplicity** of $g$ in $S$. The
**length** $|S|=\ell$ counts terms, the **support** consists of the terms
with positive multiplicity, and the **sum** is
$\sigma(S)=g_1+\ldots+g_\ell$. A sequence $T$ is a **subsequence** of $S$,
or **divides** $S$, if $\mathsf v_g(T)\le\mathsf v_g(S)$ for every $g$;
then $ST^{-1}$ is the sequence of the remaining terms. The empty sequence
$1$ has length zero and sum zero. We write
$-S=(-g_1)\cdot\ldots\cdot(-g_\ell)$.

In code, `S(g, h)` is the sequence $g\cdot h$, `S * T` is the product,
`S ** k` the $k$-th power, `S / T` is $ST^{-1}$, `-S` is the negative, and
`S.sum()` is $\sigma(S)$.

## Zero-sum sequences and atoms

A sequence $S$ is a **zero-sum sequence** if $\sigma(S)=0$, and a
**minimal zero-sum sequence** if it is a nontrivial zero-sum sequence
without proper nontrivial zero-sum subsequences. The zero-sum sequences
form a submonoid $\mathcal B(G)\subseteq\mathcal F(G)$, the **monoid of
zero-sum sequences**, whose only unit is the empty sequence. Its atoms are
precisely the minimal zero-sum sequences; the package calls them **atoms**.

The sequence $0$ of length one is an atom. The package's atom catalogues
omit it, and factorization queries take sequences over
$G^\bullet=G\setminus\{0\}$, that is, without the term $0$; this loses
nothing, since $\mathsf L(0^kS)=k+\mathsf L(S)$. Zero-sum constants allow
the term $0$.

For a subset $G_0\subseteq G$, sequences over $G_0$ are still added in $G$;
$G_0$ need not be a subgroup.

## Factorizations and length sets

A factorization of a zero-sum sequence writes it as a product of atoms.
Its **length** is the number of atoms, not the length $|S|$ of the
sequence. The **length set** $\mathsf L(S)$ contains the lengths of all
factorizations of $S$.

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

Group automorphisms act on sequences term by term, $\varphi(S)=\varphi(g_1)\cdot\ldots\cdot\varphi(g_\ell)$.
They map zero-sum sequences to zero-sum sequences and atoms to atoms, and
preserve length sets. A **canonical
representative** selects one element from an orbit.

With symmetry reduction, the extremal sequences returned for a zero-sum
constant are representatives, not necessarily the complete list. The
chosen action and representative convention are documented in the
[API reference](api/index.rst).
