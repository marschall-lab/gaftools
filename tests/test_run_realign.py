"""
Tests for 'gaftools realign'
"""

from collections import namedtuple
from gaftools.cli.realign import run_realign


def test_order_gfa(tmp_path):
    #gaftools realign alignments.gaf smallgraph.gfa reads.fa
    input_gaf = 'tests/data/alignments-graphaligner.gaf'
    output_gaf = str(tmp_path) + '/output.gaf'
    run_realign(
        gaf  = input_gaf,
        graph = 'tests/data/smallgraph.gfa',
        fasta = 'tests/data/reads.fa',
        output = output_gaf,
        ext=False
    )

    input_lines =  [l.split("\t") for l in open(input_gaf)]
    output_lines = [l.split("\t") for l in open(output_gaf)]
    assert len(input_lines) == len(output_lines)

