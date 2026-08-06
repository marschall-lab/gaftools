.. image:: https://img.shields.io/pypi/v/gaftools.svg?branch=main
    :target: https://pypi.python.org/pypi/gaftools

.. image:: https://github.com/marschall-lab/gaftools/workflows/CI/badge.svg


gaftools
========

gaftools is a fast and comprehensive toolkit designed for processing pangenome alignments. It provides various functionalities such as
indexing, sorting, realignment, viewing and statistical analysis of rGFA-based GAF files.

`Link to GitHub <https://github.com/marschall-lab/gaftools/tree/main>`_

`Documentation <https://gaftools.readthedocs.io/>`_ is available on Read The Docs. Documentation can be locally built under :code:`docs` folder
using :code:`sphinx`. :code:`sphinx` is automatically installed when installing gaftools in `dev` mode.

Publication is available `here <https://doi.org/10.1093/bioinformatics/btag406>`_.

Features
--------

* Viewing GAF files based on user-defined regions or node IDs.
* Conversion of GAF format from unstable segment coordinates to stable coordinates and vice-versa.
* New tags for GFA files for ordering of bubbles and sorting GAF files.
* Post-processing GAF alignments using Wavefront Alignment and generating basic alignment statistics.
* Easy-to-install. Available in PyPI and bioconda.
* Open Source (MIT License)

Citation
--------

Pani, S., Dabbaghie, F., Marschall, T., & Söylev, A. (2026). gaftools: a toolkit for analyzing and manipulating pangenome alignments. Bioinformatics, 42(8). https://doi.org/10.1093/bioinformatics/btag406
