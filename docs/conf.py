"""Sphinx configuration file for siRNAforge documentation."""

import sys
from pathlib import Path

from sirnaforge.__init__ import __version__

# Add the source directory to the path
sys.path.insert(0, str(Path("../src").resolve()))

# -- Project information -----------------------------------------------------
project = "siRNAforge"
copyright = "2025, Austin S. Hovland"
author = "Austin S. Hovland"
release = __version__

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.todo",
    "sphinx.ext.ifconfig",
    "sphinx_autodoc_typehints",
    "myst_parser",
    "sphinx_design",
    "sphinxcontrib.mermaid",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "branding/sirnaforge_logo_3.svg"

# Enhanced HTML theme options
html_theme_options = {
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
    "prev_next_buttons_location": "both",
    "style_external_links": True,
}

# HTML context for better navigation
html_context = {
    "display_github": True,
    "github_user": "austin-s-h",
    "github_repo": "sirnaforge",
    "github_version": "master",
    "conf_py_path": "/docs/",
}

# Additional HTML options
html_use_index = True
html_split_index = True
html_show_sourcelink = True
html_show_sphinx = False
html_copy_source = False
html_show_copyright = True

# Custom CSS and JS files
html_css_files = [
    "custom.css",
]

html_js_files = [
    "custom.js",
]

# -- Extension configuration -------------------------------------------------

# Napoleon settings (Google/NumPy docstring support)
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__",
    "show-inheritance": True,
}

# Autosummary settings
autosummary_generate = True
autosummary_imported_members = True

# Todo extension settings
todo_include_todos = True
todo_emit_warnings = False

# MyST parser settings (Markdown support)
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

myst_fence_as_directive: list[str] = []

# Mermaid configuration for GitHub Pages compatibility
mermaid_version = "10.6.1"  # Use a specific stable version
mermaid_init_js = """
mermaid.initialize({
    startOnLoad: true,
    theme: 'default',
    themeVariables: {
        primaryColor: '#007acc',
        primaryTextColor: '#ffffff',
        primaryBorderColor: '#005a9e',
        lineColor: '#333333',
        sectionBkgColor: '#f0f0f0'
    },
    flowchart: {
        useMaxWidth: true,
        htmlLabels: true
    }
});
"""

# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "biopython": ("https://biopython.org/docs/1.81/api/", None),
}

# Type hints configuration
typehints_fully_qualified = False
always_document_param_types = True
typehints_document_rtype = True
