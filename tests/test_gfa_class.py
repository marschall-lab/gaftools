from gaftools.GFA import GFA

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
    assert graph['1'].seq == "AGGTCG"
    assert graph['2'].seq == "T"
    assert graph['3'].seq == "G"
    assert graph['4'].seq == "TGGC"

    # check that connectivity is correct
    assert '1' in graph['2'].neighbors()
    assert '2' in graph['1'].neighbors()

    assert '1' in graph['3'].neighbors()
    assert '3' in graph['1'].neighbors()

    assert '2' in graph['4'].neighbors()
    assert '4' in graph['2'].neighbors()

    assert '4' in graph['3'].neighbors()
    assert '3' in graph['4'].neighbors()

    # check directionality
    assert graph['1'].in_direction('2', 1)  # node 2 connects to node 1 from end of node 1
    assert graph['2'].in_direction('1', 0)  # node 1 connects to node 2 from start of node 2
    assert graph['3'].in_direction('1', 1)  # node 1 connects to node 3 from end of node 3
    assert graph['3'].in_direction('4', 0)  # node 4 connects to node 3 from start of node 3
    assert not graph.nodes['1'].in_direction('2', 0)  # checking the negative

    # check children in a certain direction
    assert '2' in graph['1'].children(1)
    assert '4' in graph['3'].children(0)


def test_node_tags():
    graph = GFA("tests/data/test_GFA_class.gfa")
    assert 'SG:Z' in graph['1'].optional
    assert graph['1'].optional['SG:Z'] == "testing_tags"
    assert len(graph['2'].optional) == 2

def test_delete_node():
    graph = GFA("tests/data/test_GFA_class.gfa")
    del graph['1']
    assert '1' not in graph.nodes
    assert '1' not in graph['2'].neighbors()
    assert '1' not in graph['3'].neighbors()


def test_add_node():
    graph = GFA("tests/data/test_GFA_class.gfa")
    graph.add_node("5", "CCCC")
    assert '5' in graph
    assert graph['5'].seq == "CCCC"
    assert graph['5'].seq_len == 4
    # adding an edge between end of node 4 and start of node 5 with 0 overlap
    graph.add_edge('4', 1, '5', 0, 0)
    assert '5' in graph['4'].neighbors()
    assert graph['4'].in_direction('5', 1)
    assert graph['5'].in_direction('4', 0)


def test_path_extraction():
    graph = GFA("tests/data/test_GFA_class.gfa")
    assert graph.path_exists(">1>2>4")
    assert graph.extract_path(">1>2>4") == "AGGTCGTTGGC"
    assert graph.extract_path(">1<2>4") == "AGGTCGATGGC"
    assert graph.extract_path("<1") == "CGACCT"
