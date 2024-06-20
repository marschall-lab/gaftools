"""
Tests for 'gaftools sort'
"""

from gaftools.cli.sort import run_sort


def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split("\t"))

    return [parse_line(l) for l in open(filename)]


def test_sort_gzipped_gaf(tmp_path):
    """
    Testing if a gzipped GAF file works as an input
    """
    input_gaf = "tests/data/alignments-more-graphaligner.gaf.gz"
    input_gfa = "tests/data/smallgraph-ordered.gfa"
    output = str(tmp_path) + "/output.log"
    run_sort(gfa=input_gfa, gaf=input_gaf, outgaf=output)
    output_lines = parse_output(output)
    assert output_lines[0][0] == "read_s3"
    assert output_lines[1][0] == "read_s464827_revcomp"
    assert output_lines[2][0] == "read_s4_s5_100_100"
    assert output_lines[3][0] == "read_s4_s6_100_100"
    assert output_lines[4][0] == "read_s4_s5_95_105_revcomp"
    assert output_lines[5][0] == "read_s4_s5_90_110"
    assert output_lines[6][0] == "read_s7"
    assert output_lines[7][0] == "read_s8_s9"
