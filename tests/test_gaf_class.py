"""
Tests for gaf.py.
"""

from gaftools.gaf import Alignment, GAF
import pytest


@pytest.mark.parametrize("gaf_file", ["smallexample.gaf", "smallexample.gaf.gz"])
def test_parse_gaf(gaf_file):
    gaf = GAF(f"tests/data/gaf_class/{gaf_file}")
    gaf_lines = list(gaf.read_file())
    assert len(gaf_lines) == 3
    for gaf_line in gaf_lines:
        assert isinstance(gaf_line, Alignment)

    assert gaf_lines[0].query_name == "read1"
    assert gaf_lines[0].query_length == 0
    assert gaf_lines[0].query_start == 1
    assert gaf_lines[0].query_end == 2
    assert gaf_lines[0].strand == "+"
    assert gaf_lines[0].path == ">1"
    assert gaf_lines[0].path_length == 3
    assert gaf_lines[0].path_start == 4
    assert gaf_lines[0].path_end == 5
    assert gaf_lines[0].residue_matches == 6
    assert gaf_lines[0].alignment_block_length == 7
    assert gaf_lines[0].mapping_quality == 8
    assert gaf_lines[0].cigar == "11="
    assert gaf_lines[0].is_primary
    assert len(gaf_lines[0].tags.keys()) == 7
    assert "NM:i:" in gaf_lines[0].tags
    assert "AS:f:" in gaf_lines[0].tags
    assert "dv:f:" in gaf_lines[0].tags
    assert "id:f:" in gaf_lines[0].tags
    assert "cm:i:" in gaf_lines[0].tags
    assert "tp:A:" in gaf_lines[0].tags
    assert "cg:Z:" in gaf_lines[0].tags
    assert gaf_lines[0].tags["NM:i:"] == "9"
    assert gaf_lines[0].tags["AS:f:"] == "10"
    assert gaf_lines[0].tags["dv:f:"] == "0.5"
    assert gaf_lines[0].tags["id:f:"] == "0.1"
    assert gaf_lines[0].tags["cm:i:"] == "127"
    assert gaf_lines[0].tags["tp:A:"] == "P"
    assert gaf_lines[0].tags["cg:Z:"] == "11="

    assert gaf_lines[1].query_name == "read2"
    assert gaf_lines[1].query_length == 0
    assert gaf_lines[1].query_start == 1
    assert gaf_lines[1].query_end == 2
    assert gaf_lines[1].strand == "+"
    assert gaf_lines[1].path == ">s1"
    assert gaf_lines[1].path_length == 3
    assert gaf_lines[1].path_start == 4
    assert gaf_lines[1].path_end == 5
    assert gaf_lines[1].residue_matches == 6
    assert gaf_lines[1].alignment_block_length == 7
    assert gaf_lines[1].mapping_quality == 8
    assert gaf_lines[1].cigar == "11="
    assert not gaf_lines[1].is_primary
    assert len(gaf_lines[1].tags.keys()) == 7
    assert "NM:i:" in gaf_lines[1].tags
    assert "AS:f:" in gaf_lines[1].tags
    assert "dv:f:" in gaf_lines[1].tags
    assert "id:f:" in gaf_lines[1].tags
    assert "ds:Z::" in gaf_lines[1].tags
    assert "tp:A:" in gaf_lines[1].tags
    assert "cg:Z:" in gaf_lines[1].tags
    assert gaf_lines[1].tags["NM:i:"] == "9"
    assert gaf_lines[1].tags["AS:f:"] == "10"
    assert gaf_lines[1].tags["dv:f:"] == "0.5"
    assert gaf_lines[1].tags["id:f:"] == "0.1"
    assert (
        gaf_lines[1].tags["ds:Z::"]
        == "22*ga:2*ca:18*at:30*ta+[tc]tctttgtttcac[ccca]:13*ct:8*gt*ga:10*ga:14*ag:1*ag:3*ct:1*ct:13*ag*ct:11*at:5*ga:4*at:9*tg:20*ga:1*cg:2*ct:17"
    )
    assert gaf_lines[1].tags["tp:A:"] == "S"
    assert gaf_lines[1].tags["cg:Z:"] == "11="

    assert gaf_lines[2].query_name == "3"
    assert gaf_lines[2].query_length == 0
    assert gaf_lines[2].query_start == 1
    assert gaf_lines[2].query_end == 2
    assert gaf_lines[2].strand == "+"
    assert gaf_lines[2].path == ">s1"
    assert gaf_lines[2].path_length == 3
    assert gaf_lines[2].path_start == 4
    assert gaf_lines[2].path_end == 5
    assert gaf_lines[2].residue_matches == 6
    assert gaf_lines[2].alignment_block_length == 7
    assert gaf_lines[2].mapping_quality == 8
    assert gaf_lines[2].cigar == "boo="
    assert not gaf_lines[2].is_primary
    assert len(gaf_lines[2].tags.keys()) == 6
    assert "NM:i:" in gaf_lines[2].tags
    assert "AS:f:" in gaf_lines[2].tags
    assert "dv:f:" in gaf_lines[2].tags
    assert "id:f:" in gaf_lines[2].tags
    assert "tp:A:" in gaf_lines[2].tags
    assert "cg:Z:" in gaf_lines[2].tags
    assert gaf_lines[2].tags["NM:i:"] == "9"
    assert gaf_lines[2].tags["AS:f:"] == "10"
    assert gaf_lines[2].tags["dv:f:"] == "0.5"
    assert gaf_lines[2].tags["id:f:"] == "0.1"
    assert gaf_lines[2].tags["tp:A:"] == "I"
    assert gaf_lines[2].tags["cg:Z:"] == "boo="
