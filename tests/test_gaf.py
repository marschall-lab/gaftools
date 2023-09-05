"""
Tests for gaf.py.
"""

from collections import namedtuple
from gaftools.gaf import parse_gaf, Alignment


def test_parse_gaf():
    #gaftools realign alignments.gaf smallgraph.gfa reads.fa
    gaf_lines = list(parse_gaf('tests/data/alignments-graphaligner.gaf'))
    assert len(gaf_lines) == 2
    for gaf_line in gaf_lines:
        assert isinstance(gaf_line, Alignment)

    assert gaf_lines[0].query_name == 'read_s8_s9'
    assert gaf_lines[0].query_length == 398
    assert gaf_lines[0].query_start == 0
    assert gaf_lines[0].query_end == 398
    assert gaf_lines[0].strand == '+'
    assert gaf_lines[0].path == '>s8>s9'
    assert gaf_lines[0].path_length == 10109
    assert gaf_lines[0].path_start == 6166
    assert gaf_lines[0].path_end == 6564
    assert gaf_lines[0].residue_matches == 398
    assert gaf_lines[0].alignment_block_length == 398
    assert gaf_lines[0].mapping_quality == 60
    assert gaf_lines[0].is_primary == False
