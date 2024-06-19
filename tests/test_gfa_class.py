from gaftools.gfa import GFA

"""
The example graph I am using looks like this
            ---->node_2----->
           |                 |
node_1 --->                   ---->node_4
           |                 |
            --->rev(node_3)->

node 1 connects to end of node 3
node 4 connect to beginning of node 3
"""


def test_load_gfa():
    graph = GFA("tests/data/test_GFA_class.gfa")
    assert len(graph) == 4
    # check if the node sequences were loaded correctly
    assert graph["1"].seq == "AGGTCG"
    assert graph["2"].seq == "T"
    assert graph["3"].seq == "G"
    assert graph["4"].seq == "TGGC"
    # check that connectivity is correct
    assert "1" in graph["2"].neighbors()
    assert "2" in graph["1"].neighbors()

    assert "1" in graph["3"].neighbors()
    assert "3" in graph["1"].neighbors()

    assert "2" in graph["4"].neighbors()
    assert "4" in graph["2"].neighbors()

    assert "4" in graph["3"].neighbors()
    assert "3" in graph["4"].neighbors()

    # check directionality
    assert graph["1"].in_direction("2", 1)  # node 2 connects to node 1 from end of node 1
    assert graph["2"].in_direction("1", 0)  # node 1 connects to node 2 from start of node 2
    assert graph["3"].in_direction("1", 1)  # node 1 connects to node 3 from end of node 3
    assert graph["3"].in_direction("4", 0)  # node 4 connects to node 3 from start of node 3
    assert not graph.nodes["1"].in_direction("2", 0)  # checking the negative

    # check children in a certain direction
    assert "2" in graph["1"].children(1)
    assert "4" in graph["3"].children(0)

    # test file not found
    try:
        graph = GFA("test/data/nonexistant_graph.gfa")
    except FileNotFoundError:
        pass

    try:
        graph = GFA("tests/data/test_GFA_class_wrong_graph.gfa")
    except ValueError:
        pass


def test_node_tags():
    graph = GFA("tests/data/test_GFA_class.gfa")
    assert "SG" in graph["1"].tags
    assert graph["1"].tags["SG"] == ("Z", "testing_tags")
    assert len(graph["2"].tags) == 2


def test_delete_node():
    graph = GFA("tests/data/test_GFA_class.gfa")
    del graph["1"]
    assert "1" not in graph.nodes
    # checking that the edges were also removed properly
    assert "1" not in graph["2"].neighbors()
    assert "1" not in graph["3"].neighbors()


def test_add_node():
    graph = GFA("tests/data/test_GFA_class.gfa")
    node_line = "S\t5\tCCCC".split()
    graph.add_node(node_line[1], node_line[2])
    assert "5" in graph
    assert graph["5"].seq == "CCCC"
    assert graph["5"].seq_len == 4
    # adding an edge between end of node 4 and start of node 5 with 0 overlap
    graph.add_edge("4", "+", "5", "+", 0)
    assert "5" in graph["4"].neighbors()
    assert graph["4"].in_direction("5", 1)
    assert graph["5"].in_direction("4", 0)


def test_path_extraction():
    graph = GFA("tests/data/test_GFA_class.gfa")
    assert graph.path_exists(">1>2>4")
    assert graph.extract_path(">1>2>4") == "AGGTCGTTGGC"
    assert graph.extract_path(">1<2>4") == "AGGTCGATGGC"
    assert graph.extract_path("<1") == "CGACCT"

    # testing non-ACTG characters
    graph["1"].seq = graph["1"].seq + "N"
    assert graph.extract_path(">1>2>4") == graph["1"].seq + graph["2"].seq + graph["4"].seq


def test_components():
    graph = GFA("tests/data/test_GFA_class.gfa")
    components = graph.all_components()
    assert len(components) == 1
    assert components[0] == set(graph.nodes.keys())
    graph.add_node("5", "GGCC")
    graph.add_node("6", "CCGG")
    graph.add_edge("5", "+", "6", "+", 2)
    components = graph.all_components()
    assert len(components) == 2
    comp1, comp2 = components
    if len(comp1) > len(comp2):
        assert comp1 == {"1", "2", "3", "4"}
        assert comp2 == {"5", "6"}
    else:
        assert comp2 == {"1", "2", "3", "4"}
        assert comp1 == {"5", "6"}


def test_bicc():
    graph = GFA("tests/data/smallgraph-noseq.gfa")
    biccs, artic_points = graph.biccs()
    assert set(artic_points) == {"s8", "s6", "s4", "s2"}
    assert len(biccs) == 5


def test_contig_length():
    graph = GFA("tests/data/smallgraph.gfa")
    assert 21837 == graph.get_contig_length("chr1")


def test_dfs():
    graph = GFA("tests/data/test_GFA_class.gfa")
    ordered_dfs = graph.dfs("1")
    assert ordered_dfs == ["1", "3", "4", "2"]
