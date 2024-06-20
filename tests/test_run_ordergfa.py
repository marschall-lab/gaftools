"""
Tests for 'gaftools order_gfa'
"""

from gaftools.gfa import GFA
from gaftools.cli.order_gfa import run_order_gfa


def test_order_gfa(tmp_path):
    input_gfa = "tests/data/smallgraph.gfa"
    run_order_gfa(
        gfa_filename=input_gfa,
        outdir=str(tmp_path),
        chromosome_order="chr1",
        with_sequence=False,
    )
    output_gfa = str(tmp_path) + "/smallgraph-chr1.gfa"
    graph1 = GFA("tests/data/smallgraph-ordered.gfa", low_memory=True)
    graph2 = GFA(output_gfa, low_memory=True)
    assert graph1.is_equal_to(graph2)

    # testing the output CSV file that it makes sense
    with open(str(tmp_path) + "/smallgraph-chr1.csv", "r") as infile:
        for l in infile:
            if l.startswith("Name"):
                assert l == "Name,Color,SN,SO,BO,NO\n"
            if l.startswith("s9"):
                assert l == "s9,blue,chr1,18093,8,2\n"
