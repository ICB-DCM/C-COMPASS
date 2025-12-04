# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import setuptools_scm

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "C-COMPASS"
copyright = "2024, Daniel Haas"
author = "Daniel Haas"
version = release = setuptools_scm.get_version(root="..", relative_to=__file__)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-the-linkcheck-builder
linkcheck_anchors_ignore_for_url = [
    # https://github.com/sphinx-doc/sphinx/issues/9016
    r"https://github.com/ICB-DCM/C-COMPASS?",
]
linkcheck_ignore = [
    # 403 from https://www.biorxiv.org/ during linkcheck
    "https://doi.org/10.1101/2024.08.05.606647",
]
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
# html_static_path = ["_static"]
