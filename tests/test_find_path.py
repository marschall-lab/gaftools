"""
Test for gaftools find_path
"""

from gaftools.cli.find_path import run


def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split("\t"))

    return [parse_line(l) for l in open(filename)]


# testing for a file input for path list
def test_file_input(tmp_path):
    input_gfa = "tests/data/smallgraph.gfa"
    input_path_list = "tests/data/find_path-input.txt"
    output = str(tmp_path) + "/output.fasta"
    truth = "tests/data/find_path-output.fasta"
    run(input_gfa, input_path_list, output=output, fasta=True)
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
