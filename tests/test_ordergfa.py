"""
Tests for 'gaftools order-gfa'
"""

from gaftools.gfa import GFA
from gaftools.cli.order_gfa import decompose_and_order, run_order_gfa, name_comps
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
        "gfa2rgfa/reference-graph.gfa",
        "gfa2rgfa/reference-graph-seqfile.gfa",
        "gfa2rgfa/reference-graph-partial-seqfile.gfa",
        #'reference-graph-FOO#2.gfa'   <-- this is not working since the rank 0 nodes come from 2 different assemblies and hence have two SN tags
    ],
)
def test_order_gfa_customgraph(tmp_path, gfa_file):
    input_gfa = f"tests/data/{gfa_file}"
    run_order_gfa(
        gfa_filename=input_gfa,
        outdir=str(tmp_path),
        by_chrom=False,
        without_sequence=True,
    )
    outname = gfa_file.split("/")[-1].split(".")[0] + "-complete.gfa"
    output_gfa = f"{str(tmp_path)}/{outname}"
    graph1 = GFA("tests/data/order_gfa/customgraph-ordered.gfa", low_memory=True)
    graph2 = GFA(output_gfa, low_memory=True)
    compare_bo_no_tags(graph1, graph2)


def test_order_gfa_missing_sn_inside_node_still_orders(tmp_path):
    input_gfa = tmp_path / "missing-sn-inside.gfa"
    with open("tests/data/smallgraph.gfa", "r") as infile:
        lines = infile.readlines()

    replaced = False
    with open(input_gfa, "w") as outfile:
        for line in lines:
            if (
                not replaced
                and line.startswith("S\ts644045\t")
                and "SN:Z:HG01106#2#JAHAMB010000116.1" in line
            ):
                line = line.replace("\tSN:Z:HG01106#2#JAHAMB010000116.1", "", 1)
                replaced = True
            outfile.write(line)

    assert replaced

    run_order_gfa(
        gfa_filename=str(input_gfa),
        outdir=str(tmp_path),
        by_chrom=False,
        without_sequence=True,
    )

    output_gfa = str(tmp_path) + "/missing-sn-inside-complete.gfa"
    graph1 = GFA("tests/data/order_gfa/smallgraph-ordered.gfa", low_memory=True)
    graph2 = GFA(output_gfa, low_memory=True)
    compare_bo_no_tags(graph1, graph2)
    assert "SN" not in graph2.nodes["s644045"].tags


def test_order_gfa_missing_sn_scaffold_node_skips_component(tmp_path, caplog):
    input_gfa = tmp_path / "missing-sn-scaffold.gfa"
    with open("tests/data/smallgraph.gfa", "r") as infile:
        lines = infile.readlines()

    replaced = False
    with open(input_gfa, "w") as outfile:
        for line in lines:
            if not replaced and line.startswith("S\ts2\t") and "SN:Z:chr1" in line:
                line = line.replace("\tSN:Z:chr1", "", 1)
                replaced = True
            outfile.write(line)

    assert replaced

    run_order_gfa(
        gfa_filename=str(input_gfa),
        outdir=str(tmp_path),
        by_chrom=False,
        without_sequence=True,
    )

    assert "scaffold node(s) without an SN tag; cannot order this component" in caplog.text
    assert "Chromosome chr1 was skipped" in caplog.text
    assert not (tmp_path / "missing-sn-scaffold-complete.gfa").exists()


def test_order_gfa_missing_so_scaffold_node_skips_component(tmp_path, caplog):
    input_gfa = tmp_path / "missing-so-scaffold.gfa"
    with open("tests/data/smallgraph.gfa", "r") as infile:
        lines = infile.readlines()

    replaced = False
    with open(input_gfa, "w") as outfile:
        for line in lines:
            if not replaced and line.startswith("S\ts2\t") and "SO:i:3293" in line:
                line = line.replace("\tSO:i:3293", "", 1)
                replaced = True
            outfile.write(line)

    assert replaced

    run_order_gfa(
        gfa_filename=str(input_gfa),
        outdir=str(tmp_path),
        by_chrom=False,
        without_sequence=True,
    )

    assert "scaffold node(s) without an SO tag; cannot order this component" in caplog.text
    assert "Chromosome chr1 was skipped" in caplog.text
    assert not (tmp_path / "missing-so-scaffold-complete.gfa").exists()


def test_order_gfa_empty_inside_ref_warns_and_continues(caplog):
    graph = GFA()
    for node_id, tags in [
        ("a", ["SN:Z:chr1", "SO:i:0", "SR:i:0"]),
        ("b", ["SN:Z:alt", "SO:i:1", "SR:i:1"]),
        ("c", ["SN:Z:alt", "SO:i:2", "SR:i:1"]),
        ("d", ["SN:Z:chr1", "SO:i:10", "SR:i:0"]),
    ]:
        graph.add_node(node_id, seq="A", tags=tags)

    for node1, node2 in [("a", "b"), ("b", "c"), ("c", "a"), ("a", "d")]:
        graph.add_edge(node1, "+", node2, "+", 0, [0])

    scaffold_nodes, inside_nodes, node_order, bo, bubble_count = decompose_and_order(
        graph,
        {"a", "b", "c", "d"},
        "chr1",
        ignore_branching=True,
        scaffold_file=None,
        bo_start=0,
    )

    assert (
        "terminal biconnected component with no internal node matching the reference SN/SO tags"
        in caplog.text
    )
    assert scaffold_nodes == {"a", "d"}
    assert inside_nodes == {"b", "c", "d"}
    assert node_order["a"] == (1, 0)
    assert node_order["d"] == (3, 0)
    assert bo == 4
    assert bubble_count == 2


def test_name_comps_resets_current_tag():
    graph = GFA()
    graph.add_node("n1", tags=["SN:Z:chr1"])
    graph.add_node("n2")

    with pytest.raises(ValueError, match="SN tags could be missing"):
        name_comps(graph, [{"n1"}, {"n2"}])
