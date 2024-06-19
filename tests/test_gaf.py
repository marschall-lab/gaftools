"""
Tests for gaf.py.
"""

from gaftools.gaf import Alignment, GAF


def test_parse_gaf():
    # gaftools realign alignments.gaf smallgraph.gfa reads.fa
    gaf = GAF("tests/data/alignments-graphaligner.gaf")
    gaf_lines = list(gaf.read_file())
    assert len(gaf_lines) == 2
    for gaf_line in gaf_lines:
        assert isinstance(gaf_line, Alignment)

    assert gaf_lines[0].query_name == "read_s8_s9"
    assert gaf_lines[0].query_length == 398
    assert gaf_lines[0].query_start == 0
    assert gaf_lines[0].query_end == 398
    assert gaf_lines[0].strand == "+"
    assert gaf_lines[0].path == ">s8>s9"
    assert gaf_lines[0].path_length == 10109
    assert gaf_lines[0].path_start == 6166
    assert gaf_lines[0].path_end == 6564
    assert gaf_lines[0].residue_matches == 398
    assert gaf_lines[0].alignment_block_length == 398
    assert gaf_lines[0].mapping_quality == 60
    assert gaf_lines[0].cigar == "398="
    assert gaf_lines[0].is_primary

    assert gaf_lines[1].query_name == "read_s8_s9_deletion15"
    assert gaf_lines[1].query_length == 383
    assert gaf_lines[1].query_start == 0
    assert gaf_lines[1].query_end == 383
    assert gaf_lines[1].strand == "+"
    assert gaf_lines[1].path == ">s8>s9"
    assert gaf_lines[1].path_length == 10109
    assert gaf_lines[1].path_start == 6166
    assert gaf_lines[1].path_end == 6564
    assert gaf_lines[1].residue_matches == 383
    assert gaf_lines[1].alignment_block_length == 398
    assert gaf_lines[1].mapping_quality == 60
    assert gaf_lines[1].cigar == "189=2D1=13D193="
    assert gaf_lines[1].is_primary
