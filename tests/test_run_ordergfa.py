"""
Tests for 'gaftools order_gfa'
"""

from collections import namedtuple
from gaftools.cli.order_gfa import run_ordergfa


def test_ordergfa(tmp_path):
    input_gfa = 'tests/data/smallgraph.gfa'
    run_ordergfa(
        gfa_filename = 'tests/data/smallgraph.gfa',
        outdir = str(tmp_path),
        chromosome_order='chr1',
        with_sequence=False,
    )
    output_gfa = str(tmp_path) + '/smallgraph-chr1.gfa'

    input_lines =  [l.split("\t") for l in open(input_gfa)]
    output_lines = [l.split("\t") for l in open(output_gfa)]
    assert len(input_lines) == len(output_lines)

    for input_line, output_line in zip(input_lines, output_lines):
        assert input_line[0] == output_line[0]

