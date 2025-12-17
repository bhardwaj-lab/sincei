# -*- coding: utf-8 -*-

# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import re
import sys

# Ensure the repository root is on sys.path so `import sincei` works
# When building from docs/, the parent directory is the repo root
sys.path.insert(0, os.path.abspath(".."))
import sincei

VPATH = "../sincei/_version.py"


def get_version(path=VPATH):
    try:
        f = open(path)
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match('__version__ = "([^//"]+)"', line)
        if mo:
            ver = mo.group(1)
            return ver
    return None


# -- Project information -----------------------------------------------------

project = "sincei"
copyright = "2025, Vivek Bhardwaj"
author = "Vivek Bhardwaj"

# The full version, including alpha/beta/rc tags
release = get_version()


# -- General configuration ---------------------------------------------------

# To avoid issues with smart quotes in code examples (e.g., '--' turns to an em dash)
smartquotes = False

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinxarg.ext",
    "sphinx_toolbox.collapse",
    "nbsphinx",
]
#    'numpydoc'

# Do not execute tutorial notebooks
nbsphinx_execute = "never"

# This is needed to suppress autosummary reordering
numpydoc_show_class_members = False

# Order members by source order instead of alphabetically
autodoc_member_order = "bysource"

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

# If optional heavy dependencies aren't present locally, mock them so autodoc doesn't fail.
# ReadTheDocs installs extras=['doc'], but local builds might not have everything.
# autodoc_mock_imports = [
#     "numpy", "pandas", "scipy", "scanpy", "anndata", "mudata", "umap", "leidenalg",
#     "matplotlib", "networkx", "igraph", "torch", "mctorch", "deeptools", "gensim",
#     "tqdm", "deeptoolsintervals", "sklearn", "pysam", "pyBigWig", "joblib", "py2bit",
# ]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
# html_theme = "sphinx_rtd_theme"
# on_rtd = os.environ.get("READTHEDOCS", None) == "True"

# if not on_rtd:  # only import and set the theme if we're building docs locally
#     import sphinx_rtd_theme

#     html_theme = "sphinx_rtd_theme"

html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/bhardwaj-lab/sincei",
    "use_repository_button": True,
    "pygments_light_style": "tango",
    "pygments_dark_style": "monokai",
}

html_logo = "./content/images/sincei-logo-transparent.png"
html_show_sphinx = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
