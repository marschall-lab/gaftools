"""
Adds BO and NO tags to GFA
"""

import sys
import logging
from collections import namedtuple, defaultdict, Counter
from gaftools.GFA import GFA
from argparse import ArgumentParser
import pdb

default_chromosome_order = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM'
logger = logging.getLogger(__name__)

def run_order_gfa(
    gfa_filename,
    outdir,
    chromosome_order=None,
    with_sequence=False,
):

    ########################### old part
    # chromosome_order = chromosome_order.split(sep=',')

    # logger.info('Reading %s', gfa_filename)

    # to remove later
    nodes, edges = parse_gfa(gfa_filename, with_sequence)

    chromosome_order = chromosome_order.split(sep=",")

    logger.info(f"Reading {gfa_filename}")
    graph = GFA(gfa_filename, low_memory=True)

    # the __str__ functions print the number of nodes and edges
    logger.info("The graph has:")
    logger.info(graph.__str__())

    components = graph.all_components()
    logger.info(f"Connected components: {len(components)}")
    components = name_comps(graph, components)
    if set(components.keys()) == set(chromosome_order):
        logger.info("Found one connected component per expected chromosome.")
    else:
        logger.info("Chromosome set mismatch:")
        logger.info(f"  Expected: {','.join(chromosome_order)}")
        logger.info(f"  Found: {','.join(sorted(components.keys()))}")
        sys.exit(1)


    # running index for the bubble index (BO) already used
    bo = 0
    total_bubbles = 0
    for chromosome in chromosome_order:
        logger.info('Processing %s', chromosome)
        component_nodes = components[chromosome]
        # Initialize files
        f_gfa = open(outdir+'/'+gfa_filename.split("/")[-1][:-4]+'-'+chromosome+'.gfa', 'w')
        f_colors = open(outdir+'/'+gfa_filename.split("/")[-1][:-4]+'-'+chromosome+'.csv', 'w')
        f_colors.write('Name,Color,SN,SO,BO,NO\n')

        scaffold_nodes, inside_nodes, node_order, bo, bubble_count = decompose_and_order(graph, component_nodes, bo)

        # scaffold_nodes2, inside_nodes2, node_order2, bo2, bubble_count2 = decompose_and_order2(graph, component_nodes, bo)
        total_bubbles += bubble_count

        # to adjust later to write according to the GFA I have
        # need to fix that one thing abou the edge tags
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


def decompose_and_order(graph, component, bo_start=0):
    # component_nodes I already have, I can just get form components[chromosome] :D
    # scaffold_nodes are simply the articulation points that bicc will return
    # bubble_count is the number of biccs
    # inside_nodes are the component nodes without the scaffold nodes
    # the order is the tricky one that I need to solve
    logger.info(f" Input graph: {len(component)} nodes")
    if len(component) == 1:
        node = list(component)[0]
        return component, set(), {node: (bo_start, 0)}, bo_start + 1, 0
    logger.info(f" Finding Biconnected Components of the component")
    all_biccs, artic_points = graph.bicc()
    bubbles = []
    scaffold_graph = GFA()
    scaffold_node_types = dict()
    for n in artic_points:
        scaffold_graph.add_node(n)
        scaffold_node_types[n] = 's'
    inside_nodes = set()

    for bc in all_biccs:
        bc_inside_nodes = bc.difference(artic_points)
        bc_end_nodes = bc.intersection(artic_points)
        inside_nodes.update(bc_inside_nodes)

        if len(bc_inside_nodes) == 0:
            assert len(bc_end_nodes) == 2
            node1, node2 = tuple(bc_end_nodes)
            scaffold_graph.add_edge(('s',node1), ('s',node2))
        else:
            bubble_index = len(bubbles)
            bubbles.append(bc_inside_nodes)
            scaffold_graph.add_node(str(bubble_index))
            scaffold_node_types[str(bubble_index)] = 'b'
            for end_node in bc_end_nodes:
                scaffold_graph.add_edge(end_node, "+", str(bubble_index), "+", 0)
    logger.info('  Bubbles: %d', len(bubbles))
    logger.info('  Scaffold graph: %d nodes', len(scaffold_graph))

    # Find start/end points of the line by looking for nodes with degree 1
    degree_one = [x.id for x in scaffold_graph.nodes.values() if len(x.neighbors()) == 1]
    degree_two = [x.id for x in scaffold_graph.nodes.values() if len(x.neighbors()) == 2]

    # Perform DFS traversal
    assert len(degree_one) == 2
    assert len(degree_two) == len(scaffold_graph) - 2
    traversal = scaffold_graph.dfs(degree_one[0])


    traversal_scaffold_only = [node_name for node_name in traversal if scaffold_node_types[node_name] == 's']
    # check that all scaffold nodes carry the same sequence name (SN), i.e. all came for the linear reference
    assert len(set(graph[n].tags['SN'] for n in  traversal_scaffold_only)) == 1
    # I save tags as key:(type, value), so "SO":(i, '123')
    coordinates = list(int(graph[n].tags['SO'][1]) for n in  traversal_scaffold_only)

    # make sure to that the traversal is in ascending order
    if coordinates[0] > coordinates[-1]:
        traversal.reverse()
        traversal_scaffold_only.reverse()
        coordinates.reverse()
    for i in range(len(coordinates)-1):
        assert coordinates[i] < coordinates[i+1]
    # compute dictionary mapping each node name to the corresponding "bubble order" and "node order" (BO,NO)
    node_order = dict()
    bo = bo_start
    for node in traversal:
        node_type = scaffold_node_types[node]
        if node_type == 's':
            node_order[node] = (bo,0)
        elif node_type == 'b':
            for i, n in enumerate(sorted(bubbles[int(node)])):
                node_order[n] = (bo,i+1)
        else:
            assert False
        bo += 1
    return artic_points, inside_nodes, node_order, bo, len(bubbles)


def name_comps(graph, components):
    """
    the graph is GFA object and components is a list of sets of the node ids of each component
    for each component we take a majority vote of the SN tage and name the component accordingly
    """
    named_comps = dict()
    current_tag = ""
    for comp in components:
        counts = Counter([graph[n].tags["SN"][1] for n in comp])
        most_freq = 0
        for tag, count in counts.items():
            if most_freq <= count:
                current_tag, most_freq = tag, count
        if current_tag == "":
            raise ValueError("Were not able to assign a chromosome to component, SN tags could be missing")
        named_comps[current_tag] = comp
    return named_comps


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
