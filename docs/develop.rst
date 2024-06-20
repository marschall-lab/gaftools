.. _developing:

Developing
==========

Documentation is inspired from `WhatsHap <https://whatshap.readthedocs.io/en/latest/>`_

Developer's Installation
------------------------

For installation gaftools for the purpose of development,
please use a conda environment.

Conda environment can be used using these commands::

    git clone git@github.com:marschall-lab/gaftools.git
    cd gaftools
    conda create -n gaftools-dev python=3.10
    conda activate gaftools-dev
    pip install -e .[dev]


Adding a new subcommand
-----------------------

For creating a new subcommand under gaftools, add a new script under :code:`gaftools/cli/`.
Make sure to follow the same format and add test cases for the subcommand under :code:`tests`.

Since gaftools is purely written in Python, adding a script is enough. 
If the new script adds external dependencies, then add the dependencies to :code:`requirements.txt`.


Executing Test Cases
--------------------

To execute all the tests (collected in :code:`tests`) run::

    pytest

To generate an html report on the test coverage::

    coverage run -m pytest
    coverage report -m
    coverage html
    firefox  htmlcov/index.html

Before creating a pull request
    * write test cases for the new code.
    * run :code:`pytest`.


Writing Documentation
---------------------

The documentation for gaftools is written in
`reStructuredText format <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_
and is translated by `Sphinx <http://www.sphinx-doc.org/>`_ into HTML format.
The documentation is found under :code:`docs`

The documentation is hosted on `Read the Docs <https://readthedocs.org/>`_.

For testing the documentation, run developers installation and under :code:`docs`, run :code:`make html`. This creates the htmls for the
files under :code:`_build`.


Adding C++ and Cython scripts
-----------------------------

Currently, the setup script does not support building of C++ scripts and wrapping them into Python using Cython.
If it is required for any future work, refer to the WhatsHap setup structure and its documentation to add C++ and Cython scripts.

Accordingly this documentation for developing will have to be changed.