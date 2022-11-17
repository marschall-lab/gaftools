# gaftools
General purpose utility related to GAF files

## Developer Installation

```
git clone git@github.com:marschall-lab/gaftools.git
cd gaftools
python -m venv venv
source venv/bin/activate
pip install -e .
```

## Developer Editing

Since no C++ code has been implemented yet (the .cpp scripts in src directory are not utilised anywhere), you can simply edit the python scripts and then run the command. When C++ scripts get involved, to make edits and see the changes, the `pip install -e .` command has to be run again.

## Usage

Use the command line function `gaftools` inside the python virutal environment.
