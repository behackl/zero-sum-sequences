# Working on the documentation

This is the initial Shibuya scaffold on `docs/shibuya`, branched from a
checkpoint on `zero-sum-constants`. Finish the constants work first, then
rebase this branch onto its final tip before expanding the docs.

The repository's other worktree remains on `zero-sum-constants`. The
experimental Thebe and marimo sites are not part of this docs branch.

## Build locally

From the **docs worktree root**:

```console
uv sync --group docs --no-default-groups
uv run --group docs --no-default-groups sphinx-build -E -W --keep-going -b html docs docs/_build/html
uv run --group docs --no-default-groups sphinx-build -E -W --keep-going -b doctest docs docs/_build/doctest
uv run --group docs --no-default-groups python -m http.server 8768 --bind 127.0.0.1 --directory docs/_build/html
```

Open http://127.0.0.1:8768/. `-W` makes Sphinx warnings fail the build;
`-E` checks a fresh environment, including references. The doctest builder
checks the quickstart examples. No notebook runtime or Sage is needed to
build this scaffold. Without uv, install the package editable plus
`python -m pip install --group docs` (pip 25.1+) in a virtual environment,
then use `python -m sphinx` for the two builds.

The published-site deployment, marimo export and browser tests are deferred
until the constants branch is ready. The `docs` dependency group is
separate from the package's runtime dependencies and normal test tools.

## Keep the site small

- **Getting started:** a few executable examples, not a second tutorial.
- **Tutorial:** the existing notebook for now; convert to a canonical
  marimo `.py` source later, generating any compatibility exports.
- **Definitions:** only terminology needed to read the API. No literature
  review or full mathematical introduction.
- **API reference:** generated from docstrings, grouped by subject.
- **References:** one shared bibliography and software acknowledgements.

There is deliberately **no algorithms chapter**, no inline execution in
API pages, and no separate prose copy of each algorithm's explanation.

## What belongs in the docstrings

Use the existing NumPy-style sections (rendered by `sphinx.ext.napoleon`):

1. A short description of the result and its mathematical meaning.
2. Parameters, return values, and relevant failure modes.
3. A `Notes` section when needed: the core recurrence/search idea,
   completeness assumptions, expensive operations, and cache behaviour.
4. Citations beside the theorem or bound actually being used, with a
   theorem/page locator where verified. Add entries to `references.bib`.
5. A small example when it clarifies a convention or an important edge case.

Do not attribute an implementation to a paper merely because the paper
contains related mathematics. Distinguish a mathematical result, a borrowed
algorithm, and our implementation choices. A short explanation is usually
more useful than a proof or an unsubstantiated complexity claim.

For example, citations in Python docstrings use reStructuredText syntax:

```rst
Notes
-----
Describe the computational idea and its assumptions here. Attribute a
specific mathematical result with :cite:p:`gao-geroldinger-2006`.
```

In MyST pages the same citation is written
`` {cite:p}`gao-geroldinger-2006` ``. Put a precise locator in the surrounding
sentence when appropriate. Sphinx resolves both forms into the same
bibliography. Code docstrings should not acquire marimo/islands directives.

### Constants-branch docstring checklist

Use the API pages to inspect these as the constants branch is finalized:

- `zero_sum_constant`: finite allowed lengths, infinity criterion,
  heredity/bitset recurrence, symmetry semantics, representative coverage,
  and practical search cost. Do not imply that an arbitrary infinite
  iterable is an accepted input.
- `davenport_constant`: exact-formula hypotheses, formula citations versus
  exhaustive fallback, caching, and the distinction between computing
  an exact value and resolving an automatic bound.
- `AdditiveSequenceSpace.davenport_bound`: authoritative explicit bounds,
  exact/cache/formula resolution, the order/exponent fallback, the
  conservative integer rounding, and the failure policy. The existing
  source comment cites the 2012 Geroldinger–Liebmann–Philipp paper, p. 214;
  verify the statement and `d(G)` versus `D(G)` convention before promoting
  that locator into the API docs.
- Scalar wrappers and `ZeroSumConstant`: consistent notation and a clear
  definition of `.extremal`, including its meaning with no symmetry and
  when the value is infinite.
- `FactorizationSolver` and `AtomCatalogue`: sufficient bounds/catalogue
  completeness, the reduced convention, and which queries build the
  full remainder graph. Much of this is already documented; avoid
  duplicating it on a new narrative page.

Bibliographic metadata in the initial `.bib` was checked against the DOI
records. This is **not** a completed theorem-by-theorem attribution audit;
original sources and precise theorem locators still need review.

## Deferred until after constants

- Finish the citation/docstring audit on the constants branch.
- Rebase `docs/shibuya` onto the finished `zero-sum-constants` tip.
- Convert the full tutorial, adapting repeated names and stateful solver
  examples to marimo's reactive model; check plotting in Pyodide.
- Add static reading/download/browser-editing exports from one source.
- Add documentation CI and publication, update the development-release
  notice in the quickstart, and review navigation on small screens.
