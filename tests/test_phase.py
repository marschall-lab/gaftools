"""
Tests for gaftools phase
"""

from gaftools.cli.phase import (
    run,
    HaplotagInformation,
    IncorrectHaplotypeError,
    IncorrectPhaseSetError,
    DuplicateHaplotagError,
    IncompatibleHaplotagError,
)

import pytest


def parse_line(line):
    return tuple(s.strip() for s in line.split("\t"))


def parse_output(filename):
    return [parse_line(l) for l in open(filename)]


def test_HaplotagInformation():
    try:
        HaplotagInformation("chr1", "none", "10x")
        assert False  # not supposed to be executed
    except IncorrectPhaseSetError:
        # expected
        pass
    try:
        HaplotagInformation("chr1", "none", "10")
        assert False  # not supposed to be executed
    except IncompatibleHaplotagError:
        # expected
        pass
    try:
        HaplotagInformation("chr1", "H10X", "10")
        assert False  # not supposed to be executed
    except IncorrectHaplotypeError:
        # expected
        pass

    HaplotagInformation("chr1", "H10", "00")


def test_duplicate_entries():
    tsv = "tests/data/phase/duplicate-entries.tsv"
    gaf = "tests/data/index_and_view/graphaligner.gaf"
    try:
        run(gaf, tsv)
        assert False
    except DuplicateHaplotagError:
        # expected
        pass


@pytest.mark.parametrize("gaf_file", ["graphaligner.gaf", "graphaligner.gaf.gz"])
def test_customgraph_alignments(tmp_path, gaf_file):
    tsv = "tests/data/phase/customgraph-haplotags.tsv"
    gaf = f"tests/data/index_and_view/{gaf_file}"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/phase/tagged.gaf"
    run(gaf, tsv, output)
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    assert len(output_lines) == len(truth_lines)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
