"""
Tests for 'gaftools stat'
"""

from gaftools.cli.stat import run_stat
import pytest


def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split(":"))

    return [parse_line(l) for l in open(filename)]


@pytest.mark.parametrize("gaf_file", ["graphaligner.gaf", "graphaligner.gaf.gz"])
def test_stat_graphaligner(tmp_path, gaf_file):
    input_gaf = f"tests/data/stat/{gaf_file}"
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


@pytest.mark.parametrize("gaf_file", ["minigraph-stable.gaf", "minigraph-stable.gaf.gz"])
def test_stat_minigraph_stable(tmp_path, gaf_file):
    input_gaf = f"tests/data/stat/{gaf_file}"
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


@pytest.mark.parametrize("gaf_file", ["minigraph-unstable.gaf", "minigraph-unstable.gaf.gz"])
def test_stat_minigraph_unstable(tmp_path, gaf_file):
    input_gaf = f"tests/data/stat/{gaf_file}"
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
