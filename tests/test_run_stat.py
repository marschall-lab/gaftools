"""
Tests for 'gaftools stat'
"""

from gaftools.cli.stat import run_stat


def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split(":"))

    return [parse_line(l) for l in open(filename)]


def test_stat_graphaligner(tmp_path):
    input_gaf = "tests/data/alignments-graphaligner.gaf"
    output = str(tmp_path) + "/output.log"
    run_stat(
        gaf_path=input_gaf,
        cigar_stat=False,
        output=output,
    )
    output_lines = parse_output(output)
    assert output_lines[0] == ("Total alignments", "2")
    assert output_lines[1] == ("Primary", "2")
    assert output_lines[2] == ("Secondary", "0")
    assert output_lines[3] == ("Reads with at least one alignment", "2")
    assert output_lines[4] == ("Total aligned bases", "781")
    assert output_lines[5] == ("Average mapping quality", "60.0")


def test_stat_minigraph_stable(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable.gaf"
    output = str(tmp_path) + "/output.log"
    run_stat(
        gaf_path=input_gaf,
        cigar_stat=False,
        output=output,
    )
    output_lines = parse_output(output)
    assert output_lines[0] == ("Total alignments", "2")
    assert output_lines[1] == ("Primary", "2")
    assert output_lines[2] == ("Secondary", "0")
    assert output_lines[3] == ("Reads with at least one alignment", "2")
    assert output_lines[4] == ("Total aligned bases", "744")
    assert output_lines[5] == ("Average mapping quality", "60.0")


def test_stat_minigraph_unstable(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable.gaf"
    output = str(tmp_path) + "/output.log"
    run_stat(
        gaf_path=input_gaf,
        cigar_stat=False,
        output=output,
    )
    output_lines = parse_output(output)
    assert output_lines[0] == ("Total alignments", "2")
    assert output_lines[1] == ("Primary", "2")
    assert output_lines[2] == ("Secondary", "0")
    assert output_lines[3] == ("Reads with at least one alignment", "2")
    assert output_lines[4] == ("Total aligned bases", "744")
    assert output_lines[5] == ("Average mapping quality", "60.0")
