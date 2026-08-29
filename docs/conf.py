# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from importlib import metadata

# -- Project information -----------------------------------------------------

project = "sincei"
copyright = "2026, Vivek Bhardwaj"
author = metadata.metadata(project).get("Author")

# The full version, including alpha/beta/rc tags
release = metadata.version(project)


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

# Do not execute tutorial notebooks
nbsphinx_execute = "never"

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

# -- Options for HTML output -------------------------------------------------

# The theme to use for the HTML help pages.
html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/bhardwaj-lab/sincei",
    "use_repository_button": True,
    "show_toc_level": 3,
    "pygments_light_style": "tango",
    "pygments_dark_style": "monokai",
    "logo": {
        "image_light": "content/images/sincei-logo-light.png",
        "image_dark": "content/images/sincei-logo-dark.png",
    },
}
html_show_sphinx = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["content/_static"]
html_css_files = ["custom.css"]
