# Sphinx configuration

import os
import sys

from setuptools_scm import get_version

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath(os.pardir))

extensions = [
    "sphinx.ext.mathjax",
    "sphinx_issues",
    "sphinx_better_subsection",
]

templates_path = ["_templates"]

source_suffix = ".rst"

master_doc = "index"

authors = "Arda Soylev, Samarendra Pani, Fawaz Dabbaghie and Tobias Marschall"
project = "gaftools"
copyright = "2022"
version = get_version(root="..", relative_to=__file__)

# Read The Docs modifies the conf.py script and we therefore get
# version numbers like 0.12+0.g27d0d31
if os.environ.get("READTHEDOCS") == "True":
    version = ".".join(version.split(".")[:2])

html_theme = "sphinx_rtd_theme"

release = version

issues_github_path = "marschall-lab/gaftools"

exclude_patterns = ["_build"]

# The reST default role (used for this markup: `text`) to use for all
# documents.
# default_role = None

pygments_style = "sphinx"

html_static_path = ["_static"]
