gaftools
========

gaftools is a collection of programs for interacting and working with GAF files and their underlying GFA files.

`Link to GitHub <https://github.com/marschall-lab/gaftools/tree/main>`_


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

gaftools can be installed by building from source::

    pip install git+https://github.com/marschall-lab/gaftools

We recommend installing inside a conda environment to allow easy removal::

    conda create -n gaftools-env python=3.10
    conda activate gaftools-env
    pip install git+https://github.com/marschall-lab/gaftools

If you already have a clone locally, then run::

    git clone https://github.com/marschall-lab/gaftools
    conda create -n gaftools-env python=3.10
    conda activate gaftools-env
    cd gaftools
    pip install .

To remove gaftools, the conda environment needs to be removed using::

    conda env remove -n gaftools-env

gaftools can be used with python>=3.8


Requirements
------------

gaftools has the following requirement:

* python>=3.8
* pysam
* pywfa

The documentation is built using :code:`sphinx`.

The users do not need to explicitly install the requirements. The requirements are installed as part of the pip installation step.


Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   guide
   develop
   changes
