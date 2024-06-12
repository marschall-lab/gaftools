"""
The code has been taken from WhatsHap.

Setup has been modified to not use Cython and reduced installation dependencies.

Link to WhatsHap: https://github.com/whatshap/whatshap
"""

import os
from setuptools import setup


extensions = []

# Avoid compilation if we are being installed within Read The Docs
if os.environ.get("READTHEDOCS") == "True":
    install_requires = []
else:
    install_requires = []
    if os.path.exists("requirements.txt"):
        with open("requirements.txt") as f:
            for line in f:
                install_requires.append(line.strip())
    else:
        install_requires = [
            "pysam",
            "pywfa"
        ]

setup(
    use_scm_version={"write_to": "gaftools/_version.py"},
    install_requires=install_requires,
)