Installation
============

Requirements
------------

gaftools has the following requirement:

* python>=3.9
* pysam
* pywfa

The users do not need to explicitly install the requirements. The requirements are installed as part of the pip installation step.

The documentation is built using :code:`sphinx` (only installed installed automatically in development mode installation of gaftools).


Building from Source
--------------------

Unreleased development versions of :code:`gaftools` can be installed by building from source::

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

To remove :code:`gaftools`, the conda environment needs to be removed using::

    conda env remove -n gaftools-env

:code:`gaftools` can be used with python>=3.8


Installing from PyPI
--------------------

:code:`gaftools` is now available on PyPI and can be installed with::

    pip install gaftools

We recommend installation inside a clean conda environment to allow easy removal and avoid dependency issues. For example, to install :code:`gaftools` in a conda environment named :code:`gaftools-env`, run the following commands::

    conda create -n gaftools-env python=3.10
    conda activate gaftools-env
    pip install gaftools

Installing from Conda
---------------------

:code:`gaftools` is now available on bioconda and can be installed with::

    conda install bioconda::gaftools

We recommend installation inside a clean conda environment to allow easy removal and avoid dependency issues. For example, to install :code:`gaftools` in a conda environment named :code:`gaftools-env`, run the following commands::

    conda create -n gaftools-env python=3.10
    conda activate gaftools-env
    conda install bioconda::gaftools
