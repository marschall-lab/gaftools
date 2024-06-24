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

To format using :code:`ruff`, run::

    ruff format gaftools/ test/ setup.py

To check for syntax and style errors using :code:`ruff`, run::

    ruff check gaftools/ test/ setup.py

To check for proper documentation build, proper syntax and styling, and pytest on multiple python versions (installed locally), run::

    tox

Before creating a pull request

#. write test cases for the new code
#. run :code:`tox` after fixing formatting and syntax with :code:`ruff`
#. create the pull request, if the tox check passes with your local python version


Using Pre-Commit
----------------

The tox pipeline has been integrated into GitHub as a CI/CD pipeline.
Instead of doing all the checks locally, the checks can be done online.

There is also pre-commit support which also allows developers to skip the
formatting and syntax checking using :code:`ruff`.

To install :code:`pre-commit`, run::

    pip install pre-commit
    cd <directory of git repository>
    pre-commit install

Now you :code:`pre-commit` will run everytime you commit to the repository. But the pre-commit
run is restricted to the staged files. As an optional step, you can run :code:`pre-commit` on
all the file using::

    pre-commit run --all-files

When the staged files fail the conditions defined in the pre-commit conditions, :code:`ruff` formatting
will attempt to fix it. After adding the files updated by :code:`ruff`, attempt to commit again. If it fails,
then the errors have to be manually fixed.

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


Making a Release
----------------

#. Update ``CHANGES.rst``: Set the correct version number and ensure that
   all nontrivial, user-visible changes are listed.

#. Ensure you have no uncommitted changes in the working copy.

#. Run ``tox``.

#. Tag the current commit with the version number (there must be a ``v`` prefix)::

       git tag -a -m "Version 0.1" v0.1

#. Push the tag::

       git push --tags

#. Wait for the GitHub Action to finish. It will deploy the sdist and wheels to
   PyPI if everything worked correctly.

If something went wrong, fix the problem and follow the above instructions again,
but with an incremented revision in the version number. That is, go from version
x.y to x.y.1. PyPI will not allow you to change a version that has already been
uploaded.
