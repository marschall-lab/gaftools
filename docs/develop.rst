.. _developing

Developing
==========

Developer's Installation
------------------------

For installation GAFtools for the purpose of development,
please use a conda environment.

Conda environment can be used using these commands::

    git clone git@github.com:marschall-lab/gaftools.git
    cd gaftools
    conda create -n gaftools-dev python=3.10 pytest coverage
    conda activate gaftools-dev
    pip install -e .

Executing Test Cases
--------------------

To execute all the tests (collected in `tests`) run::

    pytest

To generate an html report on the test coverage::

    coverage run -m pytest
    coverage report -m
    coverage html
    firefox  htmlcov/index.html

Writing Documentation
---------------------

The documentation for GAFtools is written in
`reStructuredText format <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_
and is translated by `Sphinx <http://www.sphinx-doc.org/>`_ into HTML format.
The documentation is found under `docs`

The documentation is hosted on `Read the Docs <https://readthedocs.org/>`_.


Adding a new subcommand
-----------------------

For creating a new subcommand under gaftools, add a new script under `gaftools/cli/`.
Make sure to follow the same format and add test cases for the subcommand under `tests`.

Since GAFtools is purely written in Python, adding a script is enough. 
If the new script adds external dependencies, then add the dependencies to `requirements.txt` and `setup.py`.


Adding C++ and Cython scripts
-----------------------------

Currently, the setup script does not support building of C++ scripts and wrapping them into Python using Cython.
If it is required for any future work, refer to the WhatsHap setup structure and its documentation to add C++ and Cython scripts.

Accordingly this documentation for developing will have to be changed.