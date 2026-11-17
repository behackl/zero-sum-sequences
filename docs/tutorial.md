# Tutorial

The tutorial is an interactive [marimo](https://marimo.io) notebook that
runs in your browser. Edit any cell and the results update.

**<a href="tutorial/index.html">Open the tutorial →</a>** The first start
downloads a Python runtime and can take a minute.

It takes you from a first sequence to the main computations:

1. Sequences over $C_3$: construction, products, atoms
2. Factorizations and length sets, relations between factorizations
3. Targeted queries: a factorization of a given length, with constraints
4. The factorization graph, drawn with NetworkX
5. All atoms over $C_2\oplus C_4$, and a catalogue for reuse
6. Automorphism orbits
7. Atoms on a prescribed support in $C_3^3$
8. Zero-sum constants for a group of your choice

## Run it locally

Download the <a href="tutorial/tutorial.py" download>marimo notebook</a>
and open it with [uv](https://docs.astral.sh/uv/); it installs its own
dependencies:

```console
uvx marimo edit --sandbox tutorial.py
```

A <a href="tutorial/tutorial.ipynb" download>Jupyter version</a> is also
available.
