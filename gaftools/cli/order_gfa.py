"""
Order bubbles in the GFA by adding BO and NO tags.

The BO (Bubble Order) tags order bubbles in the GFA.
The NO (Node Order) tags order the nodes in a bubble (in a lexicographic order).
"""

import sys
import os
import logging
import time
from collections import defaultdict
from gaftools.gfa import GFA
from gaftools.utils import DEFAULT_CHROMOSOME

logger = logging.getLogger(__name__)


def run_order_gfa(
    gfa_filename,
    outdir,
    by_chrom,
    chromosome_order=None,
    without_sequence=False,
    ignore_branching=False,
    output_scaffold=None,
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
    # if true then sequences are stripped
    graph = GFA(gfa_filename, low_memory=without_sequence)
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
                    f"The chromosome name provided '{c}' did not match with a component in the graph"
                )
                logger.error(
                    f"The components found in the graph: {','.join(sorted(components.keys()))}"
                )
                logger.error("Please check the chromosome order provided")
                sys.exit(1)
    else:  # user did not give chromosome order
        chromosome_order = sorted(components.keys())
        logger.info(
            f"Chromosome order was not provided as input. Searching for the default chromosomes: {', '.join(DEFAULT_CHROMOSOME)}"
        )
        logger.info(f"Found the following component names: {', '.join(chromosome_order)}")
        logger.info(
            "Using the component names as chromosome names to write rGFA. If you wish to subset components or order them, please provide the --chromosome_order argument."
        )
    # running index for the bubble index (BO) already used
    bo = 0
    total_bubbles = 0
    out_gfa = []
    out_csv = []
    for chromosome in chromosome_order:
        logger.info("Processing %s", chromosome)
        component_nodes = components[chromosome]
        if output_scaffold:
            scaffold_file = (
                outdir
                + os.sep
                + gfa_filename.split("/")[-1].split(".")[0]
                + "_"
                + chromosome
                + "_scaffold_graph.gfa"
            )
        else:
            scaffold_file = None

        scaffold_nodes, inside_nodes, node_order, bo, bubble_count = decompose_and_order(
            graph, component_nodes, chromosome, ignore_branching, scaffold_file, bo
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
            out_gfa.append(f_gfa)
            csv_file = (
                outdir + os.sep + gfa_filename.split(os.sep)[-1][:-4] + "-" + chromosome + ".csv"
            )
            out_csv.append(csv_file)
            f_colors = open(csv_file, "w")
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
    if not by_chrom:
        if out_gfa:
            final_gfa = (
                outdir
                + os.sep
                + gfa_filename.split(os.sep)[-1].split(".")[0]
                + "-complete"
                + ".gfa"
            )
            final_csv = (
                outdir
                + os.sep
                + gfa_filename.split(os.sep)[-1].split(".")[0]
                + "-complete"
                + ".csv"
            )
            with open(final_gfa, "w") as outfile:
                # outputting all the S lines first
                for f in out_gfa:
                    with open(f, "r") as infile:
                        for l in infile:
                            if l.startswith("S"):
                                outfile.write(l)
                # outputting all the L lines
                for f in out_gfa:
                    with open(f, "r") as infile:
                        for l in infile:
                            if l.startswith("L"):
                                outfile.write(l)
                for f in out_gfa:
                    os.remove(f)

            with open(final_csv, "w") as outfile:
                for f in out_csv:
                    with open(f, "r") as infile:
                        for l in infile:
                            outfile.write(l)
                for f in out_csv:
                    os.remove(f)

    logger.info("Total bubbles: %d", total_bubbles)


def force_graph_order(
    graph, scaffold_graph, bubbles, artic_points, component_name, inside_nodes, bo_start=0
):
    """
    This function forces the order of a none-linear scaffold graph by only ordering the reference paths
    and giving -1 BO and NO tags for the rest of the paths.
    returns artic_points, inside_nodes, node_order, bo, len(bubbles)
    """

    # nodes to start the traversal from
    candidates = [x.id for x in scaffold_graph.nodes.values() if len(x.neighbors()) != 2]
    new_candidates = set()
    for n in candidates:
        for nn in scaffold_graph[n].neighbors():
            new_candidates.add(nn)
    # starting from the node next to the candidate, otherwise the traversal stops short
    # as it is designed to stop at junction points
    candidates = new_candidates

    traversals = []
    for n in candidates:
        traversals.append(scaffold_graph.dfs_line(n))

    ref_traversals = []
    non_ref_traversals = []
    # now we need to filter the traversals to only keep the ones that have the reference SN tag
    for trav in traversals:
        trav_sn_tags = set(graph[n].tags["SN"][1] for n in trav if not n.startswith("bb_"))
        if len(trav_sn_tags) != 1:
            non_ref_traversals.append(trav)
        else:
            if list(trav_sn_tags)[0] == component_name:
                ref_traversals.append(trav)
            else:
                non_ref_traversals.append(trav)

    if not ref_traversals:
        logger.error(
            f"all traversals in the graph had mixed reference articulation point, cannot order chromosome {component_name}"
        )
        return None, None, None, None, None

    # we deal with the non-reference nodes and give them -1 for BO and NO tag, un-ordered
    node_order = dict()
    for trav in non_ref_traversals:
        for n in trav:
            if not n.startswith("bb_"):
                node_order[n] = (-1, -1)
            else:
                for i, nn in enumerate(sorted(bubbles[int(n[3:])])):
                    node_order[nn] = (-1, -1)

    # sort each ref traversal according to the SO tag
    for trav in ref_traversals:
        # I already know that the traversals in ref_traversals have the reference SN tag
        trav_scaffold = [n for n in trav if not n.startswith("bb_")]
        coordinates = list(int(scaffold_graph[n].tags["SO"][1]) for n in trav_scaffold)
        # make sure that the traversal is in ascending order
        if coordinates[0] > coordinates[-1]:
            trav.reverse()
            coordinates.reverse()
        for i in range(len(coordinates) - 1):
            assert coordinates[i] < coordinates[i + 1]
        trav.append(coordinates[0])

    # we need to sort the traversals between each other using the added coordinates at the end
    ref_traversals.sort(key=lambda x: x[-1])
    # add tags
    bo = bo_start
    for trav in ref_traversals:
        # the last item was an added integer to sort the traversals, is not part of the traversal
        for node in trav[0:-1]:
            if not node.startswith("bb_"):
                node_order[node] = (bo, 0)
            else:
                for i, n in enumerate(sorted(bubbles[int(node[3:])])):
                    node_order[n] = (bo, i + 1)
            bo += 1

    return artic_points, inside_nodes, node_order, bo, len(bubbles)


def decompose_and_order(
    graph, component, component_name, ignore_branching=False, scaffold_file=None, bo_start=0
):
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
        scaffold_graph[n].tags = graph[n].tags
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
            tags = [0]
            scaffold_graph.add_edge(node1, "+", node2, "+", 0, tags)

        else:
            bubble_index = len(bubbles)
            bubbles.append(bc_inside_nodes)
            bubble_size = len(bc_inside_nodes)
            scaffold_graph.add_node("bb_" + str(bubble_index))
            scaffold_graph["bb_" + str(bubble_index)].tags = {"BS": ("i", bubble_size)}
            scaffold_node_types["bb_" + str(bubble_index)] = "b"

            # adding edges
            for end_node in bc_end_nodes:
                tags = [0]
                scaffold_graph.add_edge("bb_" + str(bubble_index), "+", end_node, "+", 0, tags)

    if scaffold_file:
        logger.info(f"Writing scaffold graph to {scaffold_file}")
        scaffold_graph.write_gfa(
            set_of_nodes=None, output_file=scaffold_file, append=False, order_bo=False
        )
    # maybe remove divergent paths where the scaffold nodes do not belong to the same SN tag that is the reference one
    logger.info(f"  Bubbles: {len(bubbles)}")
    logger.info(f"  Scaffold graph: {len(scaffold_graph)} nodes")

    # Find start/end points of the line by looking for nodes with degree 1
    degree_one = [x.id for x in scaffold_graph.nodes.values() if len(x.neighbors()) == 1]
    degree_two = [x.id for x in scaffold_graph.nodes.values() if len(x.neighbors()) == 2]

    if len(degree_one) != 2:
        if ignore_branching:
            logger.info("Scaffold graph is not a line, user chose to ignore the branching")
            return force_graph_order(
                graph, scaffold_graph, bubbles, artic_points, component_name, inside_nodes, bo_start
            )
            # function for forcing the order
        else:
            logger.error(
                f"In Chromosome {component_name}, we expect only two nodes with degree one for a line graph, that was not the case."
                "Ordering can be forced with --ignore-branching"
            )
            sys.exit(1)
            return None, None, None, None, None

    if len(degree_two) != len(scaffold_graph) - 2:
        if ignore_branching:
            logger.info("Scaffold graph is not a line, forcing the order of the graph")
            return force_graph_order(
                graph, scaffold_graph, bubbles, artic_points, component_name, inside_nodes, bo_start
            )

        else:
            logger.error(
                f"Error: In Chromosome {component_name}, the number of nodes with degree 2 did not match the expected number. Skipping this chromosome"
                "Ordering can be forced with --ignore-branching"
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
        # traversal_scaffold_only.reverse()
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
            for i, n in enumerate(sorted(bubbles[int(node[3:])])):
                node_order[n] = (bo, i + 1)
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


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    arg("--chromosome-order", default="",
        help="Order in which to arrange chromosomes in terms of BO sorting. By default, it is arranged in the lexicographic order of identified component names. "
        "Expecting comma-separated list. Example: 'chr1,chr2,chr3'")
    arg("--without-sequence", default=False, action="store_true",
        help="Strip sequences from the output graph (for less memory usage and easier visualization). Default = False")
    arg("gfa_filename", metavar="GRAPH",
        help="Input rGFA file")
    arg("--outdir", default="./out",
        help='Output Directory to store all the GFA and CSV files. Default location is a "out" folder from the directory of execution.')
    arg("--by-chrom", default=False, action="store_true",
        help="Outputs each chromosome as a separate GFA, otherwise, all chromosomes in one GFA file")
    arg("--ignore-branching", action="store_true",
        help="Force the order even when branching paths occur in the scaffold graph. Alternative alleles will not be ordered")
    arg("--output-scaffold", action="store_true",
        help="Output the scaffold graph in GFA format. The scaffold graph is the graph created from collapsing all the biconnected components.")
# fmt: on


def main(args):
    run_order_gfa(**vars(args))
