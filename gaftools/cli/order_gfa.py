"""
Adds BO and NO tags to GFA
"""

import sys
import os
import logging
import gzip
import time
from collections import namedtuple, defaultdict
from gaftools.GFA import GFA
from argparse import ArgumentParser

default_chromosome_order = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM'
logger = logging.getLogger(__name__)

def run_order_gfa(
    gfa_filename,
    outdir,
    chromosome_order=None,
    with_sequence=False,
):

    chromosome_order = chromosome_order.split(sep=",")

    if not os.path.isdir(outdir):
        logging.error(f"The directory {outdir} does not exist")
        sys.exit()
    logger.info(f"Reading {gfa_filename}")
    graph = GFA(gfa_filename, low_memory=True)
    # the __str__ functions print the number of nodes and edges
    logger.info("The graph has:")
    logger.info(graph.__str__())

    components = graph.all_components()
    logger.info(f"Connected components: {len(components)}")
    # name_comps checks the most frequent SN tag for the node in the component
    # and returns a dict of chromosome_name: {nodes...}
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
        # f_gfa = open(outdir+'/'+gfa_filename.split("/")[-1][:-4]+'-'+chromosome+'.gfa', 'w')
        f_gfa = outdir+'/'+gfa_filename.split("/")[-1][:-4]+'-'+chromosome+'.gfa'
        f_colors = open(outdir+'/'+gfa_filename.split("/")[-1][:-4]+'-'+chromosome+'.csv', 'w')
        f_colors.write('Name,Color,SN,SO,BO,NO\n')

        scaffold_nodes, inside_nodes, node_order, bo, bubble_count = decompose_and_order(graph, component_nodes, chromosome, bo)

        total_bubbles += bubble_count
        for node_name in sorted(component_nodes):
            node = graph.nodes[node_name]
            bo_tag, no_tag = node_order[node_name]

            node.tags['BO'] = ("i", bo_tag)
            node.tags['NO'] = ("i", no_tag)

            if node_name in scaffold_nodes:
                color = 'orange'
            elif node_name in inside_nodes:
                color = 'blue'
            else:
                color = 'gray'
            f_colors.write('{},{},{},{},{},{}\n'.format(node_name,color,node.tags['SN'][1], node.tags['SO'][1], bo_tag, no_tag))

        graph.write_gfa(set_of_nodes=component_nodes, output_file=f_gfa, append=False, order_bo=True)

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


def parse_tag(s):
    name, type_id, value = s.split(':')
    assert len(name) == 2
    if type_id == 'i':
        return name, int(value)
    elif type_id == 'Z':
        return name, value
    else:
        assert False

def decompose_and_order(graph, component, component_name, bo_start=0):
    """
    This function takes the graph and a component
    detects all biconnected components, order the scaffold nodes and starts ordering
    """
    logger.info(f" Input component: {len(component)} nodes")
    logger.info(f" Finding Biconnected Components of the component {component_name}")
    if len(component) == 1:
        node = list(component)[0]
        return component, set(), {node: (bo_start, 0)}, bo_start + 1, 0

    new_graph = graph.graph_from_comp(component)
    start = time.perf_counter()
    all_biccs, artic_points = new_graph.biccs()
    logger.info(f" It took {time.perf_counter() - start} seconds to find the Biconnected Components")
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
            scaffold_graph.add_edge(node1, "+", node2, "+", 0)
            # scaffold_graph.add_edge(('s',node1), ('s',node2))
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
    assert len(set(new_graph[n].tags['SN'] for n in  traversal_scaffold_only)) == 1
    # I save tags as key:(type, value), so "SO":(i, '123')
    coordinates = list(int(new_graph[n].tags['SO'][1]) for n in  traversal_scaffold_only)

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
            node_order[node] = (bo, 0)
        elif node_type == 'b':
            for i, n in enumerate(sorted(bubbles[int(node)])):
                node_order[n] = (bo, i+1)
        else:
            assert False
        bo += 1
    return artic_points, inside_nodes, node_order, bo, len(bubbles)


def count_sn(graph, comp):
    """
    counts which SN tag is the majority in that chromosome to name the chromosome accordingly
    """
    counts = defaultdict(int)
    for n in comp:
        if not "SN" in graph[n].tags:
            continue
        counts[graph[n].tags["SN"][1]] += 1
    return counts

def name_comps(graph, components):
    """
    the graph is GFA object and components is a list of sets of the node ids of each component
    for each component we take a majority vote of the SN tage and name the component accordingly
    """
    named_comps = dict()
    current_tag = ""
    for comp in components:
        counts = count_sn(graph, comp)
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
