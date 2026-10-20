# Working on the documentation

The site uses Sphinx and Shibuya, with a short notation guide and an API
reference generated from the library's docstrings. It is published at
https://zero-sum-sequences.readthedocs.io/ from `main`, configured by
`.readthedocs.yaml` (Read the Docs runs `uv sync --group docs`). The experimental Thebe
and marimo sites are not part of this documentation source.

## Build locally

From the checkout root:

```console
uv sync --locked --group docs --no-default-groups
uv run --locked --group docs --no-default-groups sphinx-build -E -W --keep-going -b html docs docs/_build/html
uv run --locked --group docs --no-default-groups sphinx-build -E -W --keep-going -b doctest docs docs/_build/doctest
uv run --locked --group docs --no-default-groups python -m http.server 8768 --bind 127.0.0.1 --directory docs/_build/html
```

Open http://127.0.0.1:8768/. `-W` makes Sphinx warnings fail the build;
`-E` checks a fresh environment, including references. The doctest builder
checks examples in the quickstart and API docstrings. No notebook runtime
or Sage is needed to build the site.

Without uv, install the package editable plus
`python -m pip install --group docs` (pip 25.1+) in a virtual environment,
then use `python -m sphinx` for the two builds. The `docs` dependency group
is separate from runtime dependencies and normal test tools.

## Keep the site small

- **Getting started:** a few executable examples, not a second tutorial.
- **Tutorial:** the existing notebook for now; replace it with one canonical
  marimo `.py` source, generating any compatibility exports.
- **Definitions:** terminology needed to read the API, not a literature review.
- **API reference:** docstrings grouped by subject, without a duplicate prose
  explanation of each algorithm or inline execution widgets.
- **References:** the shared, reviewed bibliography and acknowledgements.

Use NumPy-style docstrings. Keep mathematical, completeness, and caching
qualifications where they help users understand behavior. Cite the actual
result being used, with a theorem/page locator only where verified; related
mathematics alone does not establish algorithmic provenance.

Citations in docstrings use reStructuredText, for example
`` :cite:p:`gao-geroldinger-2006` ``. In MyST pages use
`` {cite:p}`gao-geroldinger-2006` ``. Both resolve against `references.bib`.
Do not add marimo/islands directives to library docstrings.

## Remaining documentation work

- Convert the full tutorial to marimo, adapting stateful examples and plotting
  for its reactive model and Pyodide.
- Add static, download, and browser-editing exports from that single source.
- Add a documentation build and doctest job to CI.
- Update the quickstart's development-version notice when the new APIs are
  released.
