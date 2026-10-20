"""Export the marimo tutorial next to the HTML documentation.

The canonical source is ``notebooks/tutorial.py``. After an HTML build this
writes, under ``<outdir>/tutorial/``:

``index.html``
    A static, executed reading view.
``edit/``
    A browser-editable WebAssembly notebook with precomputed outputs. It
    installs a wheel built from the same checkout, so the browser runtime
    matches the documented code. Precomputing runs the notebook in an
    isolated environment created by ``uv``, which must be on ``PATH``.
``tutorial.py``
    A download whose script header installs this commit from GitHub.
``tutorial.ipynb``
    A generated Jupyter notebook with outputs.

Set ``ZSS_SKIP_TUTORIAL_EXPORT=1`` to skip the exports in quick local builds.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK = ROOT / "notebooks" / "tutorial.py"
REPOSITORY = "https://github.com/behackl/zero-sum-sequences"
HEADER_END = "# ///\n"


def _run(*args: str, cwd: Path = ROOT) -> str:
    result = subprocess.run(
        [sys.executable, *args], cwd=cwd, capture_output=True, text=True
    )
    if result.returncode:
        raise RuntimeError(
            f"tutorial export failed: {' '.join(args)}\n{result.stdout}{result.stderr}"
        )
    return result.stdout


def _commit() -> str | None:
    commit = os.environ.get("READTHEDOCS_GIT_COMMIT_HASH")
    if commit:
        return commit
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True
        )
    except OSError:
        return None
    return result.stdout.strip() or None


def _with_package_source(source: str, specification: str) -> str:
    """Pin the notebook's ``zero-sum-sequences`` dependency in its script header."""

    if source.count(HEADER_END) != 1:
        raise RuntimeError("the tutorial must have one inline script metadata block")
    block = f"# [tool.uv.sources]\n# zero-sum-sequences = {specification}\n"
    return source.replace(HEADER_END, block + HEADER_END)


def _build_wheel(directory: Path) -> str:
    code = (
        "import sys; from hatchling.build import build_wheel; "
        "print(build_wheel(sys.argv[1]))"
    )
    return _run("-c", code, str(directory)).strip().splitlines()[-1]


def export_tutorial(destination: Path) -> None:
    source = NOTEBOOK.read_text(encoding="utf-8")
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True)

    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        wheels = work / "wheels"
        wheels.mkdir()
        wheel = _build_wheel(wheels)
        browser_notebook = work / "tutorial.py"
        browser_notebook.write_text(
            _with_package_source(source, f'{{ path = "wheels/{wheel}" }}'),
            encoding="utf-8",
        )
        _run(
            "-m", "marimo", "export", "html-wasm", str(browser_notebook),
            "-o", str(destination / "edit"), "--mode", "edit", "--execute",
            "--force",
        )
    # Instructions for marimo's AI assistant; not part of the published site.
    (destination / "edit" / "CLAUDE.md").unlink(missing_ok=True)

    _run(
        "-m", "marimo", "export", "html", str(NOTEBOOK),
        "-o", str(destination / "index.html"), "--force",
    )
    _run(
        "-m", "marimo", "export", "ipynb", str(NOTEBOOK),
        "-o", str(destination / "tutorial.ipynb"), "--include-outputs", "--force",
    )

    commit = _commit()
    download = source
    if commit:
        download = _with_package_source(
            source, f'{{ git = "{REPOSITORY}", rev = "{commit}" }}'
        )
    (destination / "tutorial.py").write_text(download, encoding="utf-8")


def _on_build_finished(app, exception) -> None:
    if exception is not None or app.builder.format != "html":
        return
    if os.environ.get("ZSS_SKIP_TUTORIAL_EXPORT"):
        return
    export_tutorial(Path(app.outdir) / "tutorial")


def setup(app):
    app.connect("build-finished", _on_build_finished)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
