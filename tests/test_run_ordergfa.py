"""
Tests for 'gaftools order-gfa'
"""

from gaftools.gfa import GFA
from gaftools.cli.order_gfa import run_order_gfa
import pytest


# compare two ordered GFAs based on the nodes BO and NO tags only
# this comparison should only consider graph topology (ideally)
def compare_bo_no_tags(graph1: GFA, graph2: GFA):
    assert len(graph1.nodes.keys()) == len(graph2.nodes.keys())
    for id in graph1.nodes.keys():
        assert id in graph2.nodes
        bo1 = graph1.nodes[id].tags["BO"]
        no1 = graph1.nodes[id].tags["NO"]
        bo2 = graph2.nodes[id].tags["BO"]
        no2 = graph2.nodes[id].tags["NO"]

        assert bo1 == bo2
        assert no1 == no2


def test_order_gfa_smallgraph(tmp_path):
    input_gfa = "tests/data/smallgraph.gfa"
    run_order_gfa(
        gfa_filename=input_gfa,
        outdir=str(tmp_path),
        by_chrom=False,
        without_sequence=True,
    )
    output_gfa = str(tmp_path) + "/smallgraph-complete.gfa"
    graph1 = GFA("tests/data/order_gfa/smallgraph-ordered.gfa", low_memory=True)
    graph2 = GFA(output_gfa, low_memory=True)
    assert graph1.is_equal_to(graph2)
    compare_bo_no_tags(graph1, graph2)

    # testing the output CSV file that it makes sense
    with open(str(tmp_path) + "/smallgraph-complete.csv", "r") as infile:
        for l in infile:
            if l.startswith("Name"):
                assert l == "Name,Color,SN,SO,BO,NO\n"
            if l.startswith("s9"):
                assert l == "s9,orange,chr1,18093,10,0\n"


@pytest.mark.parametrize(
    "gfa_file",
    [
        "reference-graph.gfa",
        "reference-graph-seqfile.gfa",
        "reference-graph-partial-seqfile.gfa",
        #'reference-graph-FOO#2.gfa'   <-- this is not working since the rank 0 nodes come from 2 different assemblies and hence have two SN tags
    ],
)
def test_order_gfa_customgraph(tmp_path, gfa_file):
    input_gfa = f"tests/data/gfa2rgfa/{gfa_file}"
    run_order_gfa(
        gfa_filename=input_gfa,
        outdir=str(tmp_path),
        by_chrom=False,
        without_sequence=True,
    )
    outname = gfa_file.split(".")[0] + "-complete.gfa"
    output_gfa = f"{str(tmp_path)}/{outname}"
    graph1 = GFA("tests/data/order_gfa/customgraph-ordered.gfa", low_memory=True)
    graph2 = GFA(output_gfa, low_memory=True)
    compare_bo_no_tags(graph1, graph2)
