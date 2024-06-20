"""
Ordeing the bubble of the GFA by adding BO and NO tags.

The BO (Bubble Order) tags order the bubbles in the GFA.
The NO (Node Order) tags order the nodes in a bubble (in a lexicographic order).
"""

import sys
import os
import logging
import time
from collections import defaultdict
from gaftools.gfa import GFA

logger = logging.getLogger(__name__)

DEFAULT_CHROMOSOME = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
    "chrM",
]


def run_order_gfa(
    gfa_filename,
    outdir,
    chromosome_order=None,
    with_sequence=False,
):
    if chromosome_order is not None:
        chromosome_order = chromosome_order.split(sep=",")

    if not os.path.isdir(outdir):
        logging.warning(f"The directory {outdir} does not exist, creating one")
        try:
            os.makedirs(outdir)
        except PermissionError:
            logging.error(f"were not able to create directory {outdir}. Permission Denied")
            sys.exit()
        except FileNotFoundError:
            logging.error(f"were not able to create directory {outdir}, File not found")
            sys.exit()
        except OSError:
            logging.error(f"were not able to create directory {outdir}, OSError")
            sys.exit()

    logger.info(f"Reading {gfa_filename}")
    if with_sequence:
        graph = GFA(gfa_filename, low_memory=False)
    else:
        graph = GFA(gfa_filename, low_memory=True)
    # the __str__ functions print the number of nodes and edges
    logger.info("The graph has:")
    logger.info(graph.__str__())

    components = graph.all_components()
    logger.info(f"Connected components: {len(components)}")

    # name_comps checks the most frequent SN tag for the node in the component
    # and returns a dict of chromosome_name: {nodes...}
    components = name_comps(graph, components)
    # if not chromosome_order is None:  # user gave a list
    if chromosome_order != [""]:  # user gave a list
        for c in chromosome_order:
            if c not in set(components.keys()):
                logger.error(
                    f"The chromosome name provided {c} did not match with a component in the graph"
                )
                logger.error(f" What was Found: {','.join(sorted(components.keys()))}")
                sys.exit(1)

    else:  # user did not give a
        try:
            assert set(components.keys()) == set(DEFAULT_CHROMOSOME)
        except AssertionError:
            logger.error(
                f"chromosome order was not provided, so the default was taken, but the default did not match"
                f" what was found in the graph, which is {','.join(sorted(components.keys()))}"
            )
            sys.exit(1)
        chromosome_order = DEFAULT_CHROMOSOME
    # running index for the bubble index (BO) already used
    bo = 0
    total_bubbles = 0
    # todo output final GFA with all the chromosomes ordered
    out_files = []
    for chromosome in chromosome_order:
        logger.info("Processing %s", chromosome)
        component_nodes = components[chromosome]
        # Initialize files
        # f_gfa = open(outdir+'/'+gfa_filename.split("/")[-1][:-4]+'-'+chromosome+'.gfa', 'w')

        scaffold_nodes, inside_nodes, node_order, bo, bubble_count = decompose_and_order(
            graph, component_nodes, chromosome, bo
        )

        # skip a chromosome if something went wrong
        if scaffold_nodes:
            f_gfa = (
                outdir
                + os.sep
                + gfa_filename.split(os.sep)[-1].split(".")[0]
                + "-"
                + chromosome
                + ".gfa"
            )
            out_files.append(f_gfa)
            f_colors = open(
                outdir + os.sep + gfa_filename.split(os.sep)[-1][:-4] + "-" + chromosome + ".csv",
                "w",
            )
            f_colors.write("Name,Color,SN,SO,BO,NO\n")
            total_bubbles += bubble_count
            for node_name in sorted(component_nodes):
                node = graph.nodes[node_name]
                bo_tag, no_tag = node_order[node_name]

                node.tags["BO"] = ("i", bo_tag)
                node.tags["NO"] = ("i", no_tag)

                if node_name in scaffold_nodes:
                    color = "orange"
                elif node_name in inside_nodes:
                    color = "blue"
                else:
                    color = "gray"

                if "SN" in node.tags:
                    sn_tag = node.tags["SN"][1]
                else:
                    sn_tag = "NA"
                if "SO" in node.tags:
                    so_tag = node.tags["SO"][1]
                else:
                    so_tag = "NA"
                f_colors.write(
                    "{},{},{},{},{},{}\n".format(node_name, color, sn_tag, so_tag, bo_tag, no_tag)
                )

            graph.write_gfa(
                set_of_nodes=component_nodes,
                output_file=f_gfa,
                append=False,
                order_bo=True,
            )

            f_colors.close()

        else:
            logger.warning(f"Chromosome {chromosome} was skipped")
    final_gfa = (
        outdir + os.sep + gfa_filename.split(os.sep)[-1].split(".")[0] + "-complete" + ".gfa"
    )
    with open(final_gfa, "w") as outfile:
        # outputting all the S lines first
        for f in out_files:
            with open(f, "r") as infile:
                for l in infile:
                    if l.startswith("S"):
                        outfile.write(l)
        # outputting all the S lines
        for f in out_files:
            with open(f, "r") as infile:
                for l in infile:
                    if l.startswith("L"):
                        outfile.write(l)

    logger.info("Total bubbles: %d", total_bubbles)


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
    logger.info(
        f" It took {time.perf_counter() - start} seconds to find the Biconnected Components"
    )
    bubbles = []
    scaffold_graph = GFA()
    scaffold_node_types = dict()
    for n in artic_points:
        scaffold_graph.add_node(n)
        scaffold_node_types[n] = "s"
    inside_nodes = set()

    for bc in all_biccs:
        # the components still have the articulation nodes in, need to be removed
        bc_inside_nodes = bc.difference(artic_points)
        bc_end_nodes = bc.intersection(artic_points)
        inside_nodes.update(bc_inside_nodes)

        if len(bc_inside_nodes) == 0:
            assert len(bc_end_nodes) == 2
            node1, node2 = tuple(bc_end_nodes)
            scaffold_graph.add_edge(node1, "+", node2, "+", 0)

        else:
            bubble_index = len(bubbles)
            bubbles.append(bc_inside_nodes)
            scaffold_graph.add_node(str(bubble_index))
            scaffold_node_types[str(bubble_index)] = "b"
            for end_node in bc_end_nodes:
                scaffold_graph.add_edge(str(bubble_index), "+", end_node, "+", 0)

    logger.info(f"  Bubbles: {len(bubbles)}")
    logger.info(f"  Scaffold graph: {len(scaffold_graph)} nodes")

    # Find start/end points of the line by looking for nodes with degree 1
    degree_one = [x.id for x in scaffold_graph.nodes.values() if len(x.neighbors()) == 1]
    degree_two = [x.id for x in scaffold_graph.nodes.values() if len(x.neighbors()) == 2]

    try:
        assert len(degree_one) == 2
    except AssertionError:
        logger.warning(
            f"Error: In Chromosome {component_name}, we found more or less than two nodes with degree 1. Skipping this chromosome"
        )
        # hacky but for now maybe ok
        return None, None, None, None, None

    try:
        assert len(degree_two) == len(scaffold_graph) - 2
    except AssertionError:
        logger.warning(
            f"Error: In Chromosome {component_name}, the number of nodes with degree 2 did not mach the expected number"
        )
        return None, None, None, None, None

    # the scaffold graph should be a line graph here
    traversal = scaffold_graph.dfs(degree_one[0])
    traversal_scaffold_only = [
        node_name for node_name in traversal if scaffold_node_types[node_name] == "s"
    ]
    # check that all scaffold nodes carry the same sequence name (SN), i.e. all came for the linear reference
    assert len(set(new_graph[n].tags["SN"] for n in traversal_scaffold_only)) == 1
    # I save tags as key:(type, value), so "SO":(i, '123')
    coordinates = list(int(new_graph[n].tags["SO"][1]) for n in traversal_scaffold_only)

    # make sure that the traversal is in ascending order
    if coordinates[0] > coordinates[-1]:
        traversal.reverse()
        traversal_scaffold_only.reverse()
        coordinates.reverse()
    for i in range(len(coordinates) - 1):
        assert coordinates[i] < coordinates[i + 1]
    # compute dictionary mapping each node name to the corresponding "bubble order" and "node order" (BO,NO)
    node_order = dict()
    bo = bo_start
    for node in traversal:
        node_type = scaffold_node_types[node]
        if node_type == "s":
            node_order[node] = (bo, 0)
        elif node_type == "b":
            for i, n in enumerate(sorted(bubbles[int(node)])):
                node_order[n] = (bo, i + 1)
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
        if "SN" not in graph[n].tags:
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
            raise ValueError(
                "Were not able to assign a chromosome to component, SN tags could be missing"
            )
        named_comps[current_tag] = comp
    return named_comps


def add_arguments(parser):
    arg = parser.add_argument
    arg(
        "--chromosome_order",
        default="",
        help="Order in which to arrange chromosomes in terms of BO sorting. "
        "Expecting comma-separated list. Default: chr1,...,chr22,chrX,chrY,chrM",
    )
    arg(
        "--with-sequence",
        default=False,
        action="store_true",
        help="Retain sequences in output (default is to strip sequences)",
    )
    arg("gfa_filename", metavar="GRAPH", help="Input rGFA file")
    arg(
        "--outdir",
        default="./out",
        help='Output Directory to store all the GFA and CSV files. Default location is a "out" folder from the directory of execution.',
    )


def main(args):
    run_order_gfa(**vars(args))
