import os
from setuptools import setup

# Avoid compilation if we are being installed within Read The Docs
if os.environ.get("READTHEDOCS") == "True":
    install_requires = []
else:
    install_requires = []
    with open("requirements.txt") as f:
        for line in f:
            install_requires.append(line.strip())

setup(
    use_scm_version={"write_to": "gaftools/_version.py"},
    install_requires=install_requires,
)
