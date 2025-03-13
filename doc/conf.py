# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import warnings

import setuptools_scm

# for hoverxref<=1.5.0
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message="The '_Opt' object tuple interface is deprecated, use attribute access instead for 'default', 'rebuild', and 'valid_types'",
)

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "C-COMPASS"
copyright = "2024, Daniel Haas"
author = "Daniel Haas"
version = release = setuptools_scm.get_version(root="..", relative_to=__file__)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["hoverxref.extension"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# hoverxref configuration
hoverxref_auto_ref = True
hoverxref_roles = [
    "term",
]
hoverxref_role_types = {
    "hoverxref": "tooltip",
    "ref": "tooltip",
    "term": "tooltip",
}

# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-the-linkcheck-builder
linkcheck_anchors_ignore_for_url = [
    # https://github.com/sphinx-doc/sphinx/issues/9016
    r"https://github.com/ICB-DCM/C-COMPASS?",
]
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
# html_static_path = ["_static"]
