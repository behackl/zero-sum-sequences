# Attribution review before release

The constants branch owns `docs/references.bib` and `docs/references.md`.
Citations live beside the corresponding claims in the Python docstrings.
This is a maintainer review checklist, not an algorithms chapter.

**Status:** reviewed and accepted by the maintainer. Per-value test citations
are deferred and are not a blocker for finalizing the constants work.

## 1. High priority: exact formulas and historical credit

Review `davenport_constant` in `src/zero_sum_sequences/invariants.py`.

- The formula `D(G) = 1 + sum(n_i - 1)` is used only with known invariant
  factors `1 < n_1 | ... | n_r` for p-groups or rank at most two.
- Gao–Geroldinger, Theorem 3.1, states `d(G) = d*(G)` and explicitly credits
  **D. Kruyswijk and J. E. Olson independently**. The docstring converts
  this using `D(G) = d(G) + 1`.
- Olson I and II are included as primary references. Their identities,
  DOIs and page ranges were checked through bibliographic records and the
  survey, **not by reading both original papers**. No original Olson theorem
  numbers are asserted. Please review the primary papers and whether you
  want an additional original citation for Kruyswijk's independent work.
- Survey theorem/definition numbers were checked in the author preprint
  linked below. Its pagination differs from the published article; no
  published survey page locators have been guessed.

## 2. High priority: the upper bound and its integer relaxation

Review `AdditiveSequenceSpace.davenport_bound` and the arithmetic comment.

- van Emde Boas–Kruyswijk, report ZW 1969-008, Theorem 7.1, **printed p. 21
  (PDF p. 23)**, states `mu(G) <= m * (1 + log(|G| / m))`.
  The scan and the definition of `mu(G)` on that page were inspected:
  `mu(G)` is our `D(G)`, and maximal element order `m` equals `exp(G)`
  for finite abelian groups.
- The survey's Theorem 3.6(1) states the equivalent inequality for `d(G)`.
  The 2012 Geroldinger–Liebmann–Philipp paper repeats it; it is **not the
  original source**. The unverified published-page-214 locator was removed.
- For integer `q = |G| / exp(G) >= 1`,
  `ceil(log2(q)) >= ln(q)`. Thus the package's bound is conservative.
  Taking the minimum with `|G|` is valid because `D(G) <= |G|`.
  Please check this conversion as well as the distinction between the
  classical estimate and the package's rounding choice.
- Group metadata and the abelian group axioms are trusted, not certified
  by the adapter. Explicit bounds are also not checked for sufficiency.

## 3. High priority: algorithmic influence and provenance

Review the related-work wording in `invariants.py` and the references page.

- Zerdoum's thesis, Chapter IV, Section 2, printed pp. 56–57, propagates
  sequences containing a target zero sum; its efficient representation is
  interval-compressed. Chapter V adapts that approach to other constants.
- Geroldinger–Liebmann–Philipp, Section 4, describes a different extension
  search using negative subsums.
- The package extends avoiding sequences, retaining orbit representatives
  and their subsequence sums/lengths. **No source has been verified as the
  provenance of this particular bitset-and-stabilizer implementation.**
- Neither a claim of algorithmic novelty nor a claim that either source
  describes this implementation has been made. Please confirm whether
  either source or previously inspected software directly influenced the
  implementation and should receive a more specific acknowledgement.

## 4. Mathematical conventions to check carefully

Review `zero_sum_constant`, `ZeroSumConstant`, and the scalar wrappers.

- Finiteness criterion for finite abelian groups: if no allowed length is
  divisible by the exponent, repeat an element of maximal order. Conversely,
  if an allowed length `m` is divisible by the exponent, a sequence of
  length `|G| * (m - 1) + 1` contains `m` equal terms by pigeonhole,
  giving a zero sum. This elementary justification is not attributed to
  an unverified theorem or paper.
- A finite truncation determines the infinite allowed-length problem only
  when the returned value does not exceed the cutoff. An infinite iterable
  is not an accepted way to specify that problem.
- Constants include the identity element. Factorization queries and
  catalogues exclude identity terms. In particular, the singleton zero
  **is an atom of the full block monoid**. Calling its omission “reduction
  of the block monoid” is misleading: the full block monoid is already
  reduced in the usual monoid-theoretic sense. That wording was corrected.
- Finite `.extremal` results contain exactly one minimum per orbit of the
  chosen group, or every extremal sequence without symmetry reduction.
  An infinite result has no stored extremal sequences.

## 5. Deferred: literature values in the tests

`tests/test_invariants.py` labels several expectations as published values
without per-value sources. The Guillot–Schmid–Zerdoum two-length formula
also lacks a complete bibliographic entry. I did not recover the previously
consulted manuscript during this pass; no unrelated paper has been
substituted for it.

The maintainer can add the missing manuscript details and per-value sources
later. This is non-blocking documentation follow-up and does not affect the
exhaustive cross-check tests. The whole test corpus should not be described
as independently literature-verified until those sources have been checked.

## Sources inspected

- [Gao–Geroldinger author preprint](https://cfc.nankai.edu.cn/_upload/article/files/c6/e1/a2c52bf04b1896f59003b5993582/5c9e49ea-af5b-44ac-b153-5dbb6d8ae9a3.pdf):
  Definition 2.1, Theorems 3.1 and 3.6, and the corresponding references.
- [van Emde Boas–Kruyswijk original report](https://ir.cwi.nl/pub/7217/7217D.pdf):
  Theorem 7.1 and its notation, including visual inspection of printed p. 21.
- [Geroldinger–Liebmann–Philipp preprint](https://arxiv.org/abs/1009.5835):
  introduction and Section 4. Published metadata was checked separately.
- [Zerdoum's complete thesis](https://theses.fr/2021PA080086/document):
  title page, Chapter IV's algorithm description, and Chapter V.

When rebasing `docs/shibuya`, retain this branch's reviewed bibliography and
references page rather than maintaining a second copy. The Shibuya source
already uses these paths; reconcile the add/add changes during the rebase.
