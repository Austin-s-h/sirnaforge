"""Sphinx configuration for siRNAforge v2 documentation."""

import sys
from pathlib import Path

# Add source to path
sys.path.insert(0, str(Path("../src").resolve()))

from sirnaforge import __version__

# -- Project information -----------------------------------------------------
project = "siRNAforge"
copyright = "2025, Austin S. Hovland"
author = "Austin S. Hovland"
release = __version__

# -- General configuration ---------------------------------------------------
extensions = [
    # Core Sphinx
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    # Markdown support
    "myst_parser",
    # UI enhancements
    "sphinx_design",
    "sphinxcontrib.mermaid",
    # Live command output
    "sphinxcontrib.programoutput",
    # Type hints
    "sphinx_autodoc_typehints",
]

templates_path = ["_templates"]
exclude_patterns = ["_build"]

# -- HTML output -------------------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "../docs/branding/sirnaforge_logo_3.svg"
html_show_sphinx = False
html_copy_source = False

html_theme_options = {
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_depth": 3,
    "prev_next_buttons_location": "both",
}

# -- Extension configuration -------------------------------------------------

# MyST parser (Markdown)
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "tasklist",
    "substitution",
]

# Napoleon (docstrings)
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# Autodoc
autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
}

# Autosummary
autosummary_generate = True

# Program output (live CLI examples)
programoutput_prompt_template = "$ {command}\n{output}"
# Allow commands to fail gracefully during docs build (e.g., network unavailable)
# Error output will be shown in the docs instead of failing the build
programoutput_use_ansi = True

# Intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}
intersphinx_timeout = 5

# Type hints
typehints_fully_qualified = False
