"""Small, static documentation site; executable tutorials are a separate export."""

import sys
from importlib.metadata import version as _distribution_version
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "_ext"))

project = "zero-sum-sequences"
author = "Benjamin Hackl"
copyright = "Benjamin Hackl"
release = _distribution_version("zero-sum-sequences")
version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinxcontrib.bibtex",
    "marimo_tutorial",
]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
root_doc = "index"
exclude_patterns = ["_build", "README.md"]
myst_enable_extensions = ["colon_fence", "dollarmath"]
myst_heading_anchors = 3

# API pages render the source docstrings, including their Notes and citations.
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
python_use_unqualified_type_names = True

bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "alpha"

html_theme = "shibuya"
html_title = project
html_theme_options = {
    "accent_color": "teal",
    "globaltoc_expand_depth": 1,
    "github_url": "https://github.com/behackl/zero-sum-sequences",
}
html_show_sourcelink = False
