"""
Adds BO and NO tags to GFA
"""

import sys
import logging
from collections import namedtuple, defaultdict
from gaftools.graph import ComponentFinder
import networkx as nx
from argparse import ArgumentParser

default_chromosome_order = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM'
logger = logging.getLogger(__name__)

def run_order_gfa(
    gfa_filename,
    outdir,
    chromosome_order=None,
    with_sequence=False,
):

    chromosome_order = chromosome_order.split(sep=',')

    logger.info('Reading %s', gfa_filename)

    nodes, edges = parse_gfa(gfa_filename, with_sequence)

    logger.info('Nodes: %d', len(nodes))
    logger.info('Edges: %d', len(edges))

    cf = ComponentFinder(nodes.keys())
    for (from_node,to_node),e in edges.items():
        cf.merge(from_node,to_node)

    connected_components = set((cf.find(node) for node in nodes.keys()))
    logger.info('Connected components: %d', len(connected_components))


    name_to_component = dict((name,component) for (component,name) in component_names(cf, nodes, connected_components).items())
    if set(name_to_component.keys()) == set(chromosome_order):
        logger.info('Found one connected component per expected chromosome.')
    else:
        logger.info('Chromsome set mismatch:')
        logger.info('  Expected: %s', ','.join(chromosome_order))
        logger.info('  Found: %s', ','.join(sorted(name_to_component.keys())))
        sys.exit(1)
        
    # running index for the bubble index (BO) already used
    bo = 0
    total_bubbles = 0
    for chromosome in chromosome_order:
        logger.info('Processing %s', chromosome)
        representative_node = name_to_component[chromosome]

        # Initialize files
        f_gfa = open(outdir+'/'+gfa_filename.split("/")[-1][:-4]+'-'+chromosome+'.gfa', 'w')
        f_colors = open(outdir+'/'+gfa_filename.split("/")[-1][:-4]+'-'+chromosome+'.csv', 'w')
        f_colors.write('Name,Color,SN,SO,BO,NO\n')

        component_nodes = set()
        for node_name in sorted(nodes.keys()):
            if cf.find(node_name) == representative_node:
                component_nodes.add(node_name)

        scaffold_nodes, inside_nodes, node_order, bo, bubble_count = decompose_and_order(nodes, edges, component_nodes, bo)
        total_bubbles += bubble_count

        for node_name in sorted(component_nodes):
            node = nodes[node_name]
            bo_tag, no_tag = node_order[node_name]
            node.tags['BO'] = bo_tag
            node.tags['NO'] = no_tag
            f_gfa.write(node.to_line() + '\n')

            if node_name in scaffold_nodes:
                color = 'orange'
            elif node_name in inside_nodes:
                color = 'blue'
            else:
                color = 'gray'
            f_colors.write('{},{},{},{},{},{}\n'.format(node_name,color,node.tags['SN'], node.tags['SO'], bo_tag, no_tag))

        for edge_key in sorted(edges.keys()):
            from_node, to_node = edge_key
            if cf.find(from_node) == representative_node:
                for edge in edges[edge_key]:
                    f_gfa.write(edge.to_line() + '\n')

        f_gfa.close()
        f_colors.close()

    logger.info('Total bubbles: %d', total_bubbles)


def tag_to_str(tag):
    name, value = tag
    if type(value) is int:
        return ':'.join([name,'i',str(value)])
    elif type(value) is str:
        return ':'.join([name,'Z',value])
    else:
        assert False

class Node:
    def __init__(self, name, tags, sequence=None):
        self.name = name
        self.tags = tags
        self.sequence = sequence
    def to_line(self):
        return '\t'.join(['S', self.name, '*' if self.sequence is None else self.sequence] + [tag_to_str(t) for t in self.tags.items()])

class Edge:
    def __init__(self, from_node, from_dir, to_node, to_dir, overlap, tags):
        self.from_node = from_node
        self.from_dir = from_dir
        self.to_node = to_node
        self.to_dir = to_dir
        self.overlap = overlap
        self.tags = tags
    def to_line(self):
        return '\t'.join(['L', self.from_node, self.from_dir, self.to_node, self.to_dir, self.overlap,] + [tag_to_str(t) for t in self.tags.items()])

def parse_tag(s):
    name, type_id, value = s.split(':')
    assert len(name) == 2
    if type_id == 'i':
        return name, int(value)
    elif type_id == 'Z':
        return name, value
    else:
        assert False

def parse_gfa(gfa_filename, with_sequence=False):
    nodes = {}
    edges = defaultdict(list)

    for nr, line in enumerate(open(gfa_filename)):
        fields = line.split('\t')
        if fields[0] == 'S':
            name = fields[1]
            tags = dict(parse_tag(s) for s in fields[3:])
            sequence = None
            if with_sequence and (fields[2] != '*'):
                sequence = fields[2]
            nodes[name] = Node(name,tags,sequence)
        elif fields[0] == 'L':
            from_node = fields[1]
            from_dir = fields[2]
            to_node = fields[3]
            to_dir = fields[4]
            overlap = fields[5]
            tags = dict(parse_tag(s) for s in fields[6:])
            e = Edge(from_node,from_dir,to_node,to_dir,overlap, tags)
            edges[(from_node,to_node)].append(e)

    return nodes, edges


def decompose_and_order(nodes, edges, node_subset, bubble_order_start=0):
    logger.info('  Input graph: %d nodes', len(node_subset))
    if len(node_subset) == 1:
        node = list(node_subset)[0]
        scaffold_nodes = set([node])
        inside_nodes = set()
        node_order = {node: (bubble_order_start,0)}
        bo = bubble_order_start + 1
        return scaffold_nodes, inside_nodes, node_order, bo, 0

    graph = nx.Graph()
    for node_name in node_subset:
        graph.add_node(node_name)
    for (from_node,to_node) in edges.keys():
        if (from_node in node_subset) and (to_node in node_subset):
            graph.add_edge(from_node, to_node)

    # We create a scaffold graph with two types of nodes:
    # s-nodes (scaffold) corresponding to the articulation points and
    # b-nodes corresponding to bubbles (= biconnected components without the articulation nodes)
    # if the graph is "overall linear" then the scaffold graph is a line
    scaffold_graph = nx.Graph()
    # list of sets of nodes representing the bubbles
    bubbles = []
    scaffold_nodes = set(nx.articulation_points(graph))
    inside_nodes = set()
    for bc in nx.biconnected_components(graph):
        bc_inside_nodes = bc.difference(scaffold_nodes)
        bc_end_nodes = bc.intersection(scaffold_nodes)
        inside_nodes.update(bc_inside_nodes)
        if len(bc_inside_nodes) == 0:
            assert len(bc_end_nodes) == 2
            node1, node2 = tuple(bc_end_nodes)
            scaffold_graph.add_edge(('s',node1), ('s',node2))
        else:
            bubble_index = len(bubbles)
            bubbles.append(bc_inside_nodes)
            for end_node in bc_end_nodes:
                scaffold_graph.add_edge(('s',end_node), ('b',bubble_index))
    logger.info('  Bubbles: %d', len(bubbles))
    logger.info('  Scaffold graph: %d nodes', len(scaffold_graph.nodes))

    # Find start/end points of the line by looking for nodes with degree 1
    degree_one = list(node for node in scaffold_graph.nodes if scaffold_graph.degree(node) == 1)
    degree_two = list(node for node in scaffold_graph.nodes if scaffold_graph.degree(node) == 2)
    # Perform DFS traversal
    if len(scaffold_graph.nodes) == 1:
        dfs_tree = nx.dfs_tree(scaffold_graph, source=scaffold_graph.nodes[0])
    else:
        # For a line graph, there should be exactly two nodes with degree one and all other nodes have degree two.
        assert len(degree_one) == 2
        assert len(degree_two) == len(scaffold_graph.nodes) - 2
        dfs_tree = nx.dfs_tree(scaffold_graph, source=degree_one[0])
    traversal = list(dfs_tree.nodes())
    traversal_scaffold_only = [node_name for (node_type,node_name) in traversal if node_type == 's']
    # check that all scaffold nodes carry the same sequence name (SN), i.e. all came for the linear reference
    assert len(set(nodes[n].tags['SN'] for n in  traversal_scaffold_only)) == 1
    coordinates = list(nodes[n].tags['SO'] for n in  traversal_scaffold_only)
    # make sure to that the traversal is in ascending order
    if coordinates[0] > coordinates[-1]:
        traversal.reverse()
        traversal_scaffold_only.reverse()
        coordinates.reverse()
    for i in range(len(coordinates)-1):
        assert coordinates[i] < coordinates[i+1]
    # compute dictionary mapping each node name to the corresponding "bubble order" and "node order" (BO,NO)
    node_order = dict()
    bo = bubble_order_start
    for node in traversal:
        node_type, node_name = node
        if node_type == 's':
            node_order[node_name] = (bo,0)
        elif node_type == 'b':
            for i,n in enumerate(sorted(bubbles[node_name])):
                node_order[n] = (bo,i+1)
        else:
            assert False
        bo += 1

    return scaffold_nodes, inside_nodes, node_order, bo, len(bubbles)


def component_names(cf, nodes, connected_components):
    """
    Returns dictionary mapping names of representative nodes to
    names in SN tags occuring most often in that component.
    """
    d = {}
    for representative_node in connected_components:
        sn_counts = defaultdict(int)
        for node in nodes.values():
            if cf.find(node.name) == representative_node:
                sn_counts[node.tags['SN']]+=1
        component_name = sorted(sn_counts.items(), key=lambda t: t[1], reverse=True)[0][0]
        d[representative_node] = component_name
    return d


def add_arguments(parser):
    arg = parser.add_argument
    arg('--chromosome_order', default=default_chromosome_order,
        help='Order in which to arrange chromosomes in terms of BO sorting. '
        'Expecting comma-separated list. Default: chr1,...,chr22,chrX,chrY,chrM')
    arg('--with-sequence', default=False, action='store_true',
        help='Retain sequences in output (default is to strip sequences)')
    arg('gfa_filename', metavar='GRAPH', help='Input GFA file')
    arg('--outdir', default="./out", help='Output Directory to store all the GFA and CSV files. Default location is a "out" folder from the directory of execution.')


def main(args):
    run_order_gfa(**vars(args))
