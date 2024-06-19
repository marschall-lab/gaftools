GAFtools
========

GAFtools is a collection of programs for interacting and working with GAF files and their underlying GFA files.

`Link to GitHub <https://github.com/marschall-lab/gaftools/tree/main>`_
`Link to Preprint <>`_ TODO


Features
--------

* Viewing GAF files based on user-defined regions or node IDs.
* Conversion of GAF format from unstable segment coordinates to stable coordinates and vice-versa.
* New tags for GFA files for ordering of bubbles and sorting GAF files.
* Post-processing GAF alignments using Wavefront Alignment and generating basic alignment statistics.
* Easy-to-install
* Open Source (MIT License)


Installation
------------

GAFtools can be installed by building from source::

    pip install git+https://github.com/marschall-lab/gaftools

We recommend installing inside a conda environment to allow easy removal::

    conda create -n gaftools-env python=3.10
    conda activate gaftools-env
    pip install git+https://github.com/marschall-lab/gaftools

To remove GAFtools, the conda environment needs to be removed using::

    conda env remove -n gaftools-env

To install GAFtools in editable mode, the following steps can be followed::

    git clone git@github.com:marschall-lab/gaftools.git
    cd gaftools
    conda create -n gaftools-env python=3.10
    conda activate gaftools-env
    pip install -e .[dev]


Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   guide
   develop
   changes
