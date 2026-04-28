"""
Test for gaftools find_path
"""

from gaftools.cli.find_path import run
from gaftools.cli.order_gfa import run_order_gfa
from gaftools.errors import CommandLineError
import pytest


def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split("\t"))

    return [parse_line(l) for l in open(filename)]


def test_smallgraph(tmp_path):
    input_gfa = "tests/data/smallgraph.gfa"
    input_path_list = "tests/data/find_path/smallgraph-input.txt"
    output = str(tmp_path) + "/output.fasta"
    truth = "tests/data/find_path/smallgraph-output.fasta"
    run(
        input_gfa, path=None, paths_file=input_path_list, keep_going=True, output=output, fasta=True
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


@pytest.mark.parametrize(
    "gfa_file",
    [
        "customgraph.gfa",
        "customgraph.gfa.gz",
        "gfa2rgfa/reference-graph.gfa",
    ],
)
def test_customgraph(tmp_path, gfa_file):
    input_gfa = f"tests/data/{gfa_file}"
    input_path_list = "tests/data/find_path/customgraph-input.txt"
    output = str(tmp_path) + "/output.fasta"
    truth = "tests/data/find_path/customgraph-output.fasta"
    run(
        input_gfa, path=None, paths_file=input_path_list, keep_going=True, output=output, fasta=True
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


@pytest.mark.parametrize("gfa_file", ["smallgraph-ordered.gfa", "customgraph-ordered.gfa"])
def test_orderedgraph1(tmp_path, gfa_file):
    input_gfa = f"tests/data/order_gfa/{gfa_file}"
    input_path_list = (
        "tests/data/find_path/customgraph-input.txt"
        if "custom" in gfa_file
        else "tests/data/find_path/smallgraph-input.txt"
    )
    output = str(tmp_path) + "/output.fasta"
    truth = (
        "tests/data/find_path/customgraph-output.fasta"
        if "custom" in gfa_file
        else "tests/data/find_path/smallgraph-output.fasta"
    )
    run(
        input_gfa, path=None, paths_file=input_path_list, keep_going=True, output=output, fasta=True
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


def test_orderedgraph2(tmp_path):
    input_gfa = "tests/data/gfa2rgfa/reference-graph.gfa"
    run_order_gfa(
        gfa_filename=input_gfa,
        outdir=str(tmp_path),
        by_chrom=False,
        without_sequence=True,
    )
    ordered_gfa = str(tmp_path) + "/reference-graph-complete.gfa"
    input_path_list = "tests/data/find_path/customgraph-input.txt"
    output = str(tmp_path) + "/output.fasta"
    truth = "tests/data/find_path/customgraph-output.fasta"
    run(
        ordered_gfa,
        path=None,
        paths_file=input_path_list,
        keep_going=True,
        output=output,
        fasta=True,
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# the path extraction is tested in test_gfa_class.
# this test exists to see if the code around that extraction function works
@pytest.mark.parametrize(
    "gfa_file",
    [
        "customgraph.gfa",
        "customgraph.gfa.gz",
        "gfa2rgfa/reference-graph.gfa",
    ],
)
def test_path_as_text(tmp_path, gfa_file):
    input_gfa = f"tests/data/{gfa_file}"
    output = str(tmp_path) + "/output.fasta"
    run(input_gfa, path=">s1", keep_going=True, output=output, fasta=False)
    output_lines = parse_output(output)
    assert len(output_lines) == 1
    assert output_lines[0][0] == "ACTTATCCCC"

    run(input_gfa, path="<s1", keep_going=True, output=output, fasta=False)
    output_lines = parse_output(output)
    assert len(output_lines) == 1
    assert output_lines[0][0] == "GGGGATAAGT"

    run(input_gfa, path=">s1>s4", keep_going=True, output=output, fasta=False)
    output_lines = parse_output(output)
    assert len(output_lines) == 0

    try:
        run(input_gfa, path=">s1>s4", keep_going=False, output=output, fasta=False)
        assert False
    except CommandLineError:
        # expected
        pass
