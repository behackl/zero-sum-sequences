# Tutorial

The tutorial is a [marimo](https://marimo.io) notebook. It covers sequence
construction, factorizations and targeted witnesses, a NetworkX/Matplotlib
view of the factorization graph, atom catalogues, automorphism orbits, and
zero-sum constants. It uses ordinary Python, without Sage.

- <a href="tutorial/index.html">Read the tutorial</a>, with precomputed outputs.
- <a href="tutorial/edit/index.html">Open it in the browser</a> to edit and
  rerun the examples. Python runs on your device via Pyodide; the first start
  downloads the runtime and can take a minute. The installed package is built
  from the same commit as this documentation.
- Download it as a <a href="tutorial/tutorial.py" download>marimo notebook</a>
  or a <a href="tutorial/tutorial.ipynb" download>Jupyter notebook</a>.

The downloaded marimo notebook declares its dependencies, including this
version of the package, so it runs in an isolated environment with
[uv](https://docs.astral.sh/uv/):

```console
uvx marimo edit --sandbox tutorial.py
```

From a checkout, install the tutorial dependencies and open the source:

```console
python -m pip install --editable '.[tutorial]'
marimo edit notebooks/tutorial.py
```
