"""
Order bubbles in the GFA by adding BO and NO tags.

The BO (Bubble Order) tags order bubbles in the GFA.
The NO (Node Order) tags order the nodes in a bubble (in a lexicographic order).
"""

import os
import logging
import time
from collections import defaultdict
from gaftools.gfa import GFA
from gaftools.utils import DEFAULT_CHROMOSOME
from gaftools.errors import ChromosomeNotFoundError, BranchedGfaComponentError

logger = logging.getLogger(__name__)


def _format_node_examples(node_ids, limit=5):
    examples = node_ids[:limit]
    if len(node_ids) > limit:
        examples.append("...")
    return ",".join(examples)


def _has_required_scaffold_sn(graph, scaffold_node_types, component_name):
    # Only scaffold nodes must carry SN for component ordering to proceed.
    scaffold_nodes = [n for n, node_type in scaffold_node_types.items() if node_type == "s"]
    missing_sn_nodes = [n for n in scaffold_nodes if "SN" not in graph[n].tags]
    if missing_sn_nodes:
        logger.error(
            "Chromosome %s has %d scaffold node(s) without an SN tag; cannot order this component. "
            "Example node IDs: %s",
            component_name,
            len(missing_sn_nodes),
            _format_node_examples(missing_sn_nodes),
        )
        return False
    return True


def _has_required_articulation_so(graph, artic_points, component_name):
    missing_so_nodes = [n for n in artic_points if "SO" not in graph[n].tags]
    if missing_so_nodes:
        logger.error(
            "Chromosome %s has %d scaffold node(s) without an SO tag; cannot order this component. "
            "Example node IDs: %s",
            component_name,
            len(missing_so_nodes),
            _format_node_examples(missing_so_nodes),
        )
        return False
    return True


def _order_two_node_component(graph, component, component_name, bo_start):
    if not _has_required_articulation_so(graph, component, component_name):
        return None, None, None, None, None

    ordered_nodes = sorted(component, key=lambda x: int(graph[x].tags["SO"][1]))
    node_order = {
        ordered_nodes[0]: (bo_start, 0),
        ordered_nodes[1]: (bo_start + 1, 0),
    }
    return set(component), set(), node_order, bo_start + 2, 0


def _mark_component_unordered(component, bo_start):
    node_order = {node: (-1, -1) for node in component}
    return set(), set(component), node_order, bo_start, 0


def _mark_nodes_unordered(node_order, unordered_nodes):
    for node in unordered_nodes:
        node_order[node] = (-1, -1)
    return node_order


def _order_traversal_by_so(scaffold_graph, trav):
    # dfs_line gives us the nodes in a linear stretch, but not necessarily in reference order.
    # Order the scaffold nodes by SO, then place each collapsed bubble between the scaffold
    # nodes it connects to in this stretch.
    trav_nodes = set(trav)
    trav_scaffold = sorted(
        [n for n in trav if not n.startswith("bb_")],
        key=lambda n: int(scaffold_graph[n].tags["SO"][1]),
    )
    scaffold_pos = {n: i for i, n in enumerate(trav_scaffold)}
    last_pos = len(trav_scaffold) - 1
    sort_key = {n: float(i) for i, n in enumerate(trav_scaffold)}

    for node in trav:
        if not node.startswith("bb_"):
            continue
        neigh_pos = sorted(
            scaffold_pos[n]
            for n in scaffold_graph[node].neighbors()
            if n in trav_nodes and not n.startswith("bb_")
        )
        if len(neigh_pos) == 2:
            sort_key[node] = (neigh_pos[0] + neigh_pos[1]) / 2
        elif len(neigh_pos) == 1 and neigh_pos[0] == 0:
            sort_key[node] = -0.5
        elif len(neigh_pos) == 1 and neigh_pos[0] == last_pos:
            sort_key[node] = last_pos + 0.5
        else:
            raise ValueError("Could not place a bubble node relative to the scaffold SO order")

    return sorted(trav, key=lambda n: sort_key[n])


def run_order_gfa(
    gfa_filename,
    outdir="./out",
    by_chrom=False,
    chromosome_order="",
    without_sequence=False,
    ignore_branching=False,
    output_scaffold=None,
):
    chromosome_order = chromosome_order.split(sep=",")
    if not os.path.isdir(outdir):
        logging.warning(f"The directory {outdir} does not exist, creating one")
        try:
            os.makedirs(outdir)
        except PermissionError:
            raise PermissionError(f"were not able to create directory {outdir}. Permission Denied")
        except FileNotFoundError:
            raise FileNotFoundError(f"were not able to create directory {outdir}, File not found")
        except OSError:
            raise OSError(f"were not able to create directory {outdir}, OSError")

    logger.info(f"Reading {gfa_filename}")
    # if true then sequences are stripped
    if not os.path.exists(gfa_filename):
        raise FileNotFoundError(f"The file {gfa_filename} does not exist")
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
                raise ChromosomeNotFoundError(
                    f"The chromosome name provided '{c}' did not match with a component in the graph.\n"
                    f"The components found in the graph: {','.join(sorted(components.keys()))}\n"
                    "Please check the chromosome order provided"
                )
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

        # skip a chromosome only if ordering failed completely
        if node_order is not None:
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
        # Each candidate traversal should be collected independently.
        scaffold_graph.set_visited(False)
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
        trav[:] = _order_traversal_by_so(scaffold_graph, trav)
        # I already know that the traversals in ref_traversals have the reference SN tag
        trav_scaffold = [n for n in trav if not n.startswith("bb_")]
        coordinates = list(int(scaffold_graph[n].tags["SO"][1]) for n in trav_scaffold)
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
    if len(component) == 2:
        return _order_two_node_component(graph, component, component_name, bo_start)

    new_graph = graph.graph_from_comp(component)
    start = time.perf_counter()
    all_biccs, artic_points = new_graph.biccs()
    logger.info(
        f" It took {time.perf_counter() - start} seconds to find the Biconnected Components"
    )
    if not artic_points:
        logger.warning(
            "Chromosome %s has no articulation points after biconnected decomposition; "
            "likely a circular or fully biconnected component; emitting BO=-1 and NO=-1 "
            "for all nodes in this component.",
            component_name,
        )
        return _mark_component_unordered(component, bo_start)
    if not _has_required_articulation_so(new_graph, artic_points, component_name):
        return None, None, None, None, None
    highest_artic_point = max(artic_points, key=lambda x: int(new_graph[x].tags["SO"][1]))
    lowest_artic_point = min(artic_points, key=lambda x: int(new_graph[x].tags["SO"][1]))

    bubbles = []
    scaffold_graph = GFA()
    scaffold_node_types = dict()

    for n in artic_points:
        scaffold_graph.add_node(n)
        scaffold_graph[n].tags = graph[n].tags
        scaffold_node_types[n] = "s"

    inside_nodes = set()
    unordered_nodes = set()
    for bc in all_biccs:
        # the components still have the articulation nodes in, need to be removed
        bc_inside_nodes = bc.difference(artic_points)
        bc_end_nodes = bc.intersection(artic_points)
        ref_end_nodes = {
            n
            for n in bc_end_nodes
            if "SN" in graph[n].tags and graph[n].tags["SN"][1] == component_name
        }
        inside_nodes.update(bc_inside_nodes)

        if len(ref_end_nodes) == 0 or (len(ref_end_nodes) == 1 and len(bc_end_nodes) > 1):
            unordered_nodes.update(bc)
            inside_nodes.update(bc_end_nodes)
            continue

        if len(bc_inside_nodes) == 0:
            if len(bc_end_nodes) != 2:
                unordered_nodes.update(bc)
                inside_nodes.update(bc_end_nodes)
                continue
            node1, node2 = tuple(bc_end_nodes)
            tags = [0]
            scaffold_graph.add_edge(node1, "+", node2, "+", 0, tags)

        else:
            if len(bc_end_nodes) == 1:
                # For terminal reference blocks, try to promote the closest reference-tagged
                # internal node into the scaffold so the chromosome ends remain orderable.
                to_adjust = None
                inside_ref = [
                    n
                    for n in bc_inside_nodes
                    if "SN" in graph[n].tags
                    and graph[n].tags["SN"][1] == component_name
                    and "SO" in graph[n].tags
                ]
                if bc_end_nodes == {lowest_artic_point} and inside_ref:
                    to_adjust = min(inside_ref, key=lambda x: int(new_graph[x].tags["SO"][1]))
                elif bc_end_nodes == {highest_artic_point} and inside_ref:
                    to_adjust = max(inside_ref, key=lambda x: int(new_graph[x].tags["SO"][1]))
                elif bc_end_nodes in [{lowest_artic_point}, {highest_artic_point}]:
                    logger.warning(
                        "Chromosome %s has a terminal biconnected component with no internal node "
                        "matching the reference SN/SO tags; leaving it unpromoted in the scaffold",
                        component_name,
                    )

                if to_adjust:
                    artic_points.add(to_adjust)
                    scaffold_graph.add_node(to_adjust)
                    scaffold_graph[to_adjust].tags = graph[to_adjust].tags
                    scaffold_node_types[to_adjust] = "s"
                    bc_inside_nodes.remove(to_adjust)
                    bc_end_nodes.add(to_adjust)

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

    for node in unordered_nodes:
        if node in scaffold_node_types and len(scaffold_graph[node].neighbors()) == 0:
            del scaffold_graph[node]
            del scaffold_node_types[node]

    if len(scaffold_graph) == 0:
        logger.warning(
            "Chromosome %s has no scaffold blocks with at least two reference endpoints; "
            "emitting BO=-1 and NO=-1 for all nodes in this component.",
            component_name,
        )
        return _mark_component_unordered(component, bo_start)

    if not _has_required_scaffold_sn(graph, scaffold_node_types, component_name):
        return None, None, None, None, None

    # import pdb
    # pdb.set_trace()
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
            scaffold_nodes, inside_nodes, node_order, bo, bubble_count = force_graph_order(
                graph, scaffold_graph, bubbles, artic_points, component_name, inside_nodes, bo_start
            )
            if node_order is not None:
                _mark_nodes_unordered(node_order, unordered_nodes)
            return scaffold_nodes, inside_nodes, node_order, bo, bubble_count
            # function for forcing the order
        else:
            raise BranchedGfaComponentError(
                f"In Chromosome {component_name}, we expect only two nodes with degree one for a line graph, that was not the case. "
                "Ordering can be forced with --ignore-branching"
            )

    if len(degree_two) != len(scaffold_graph) - 2:
        if ignore_branching:
            logger.info("Scaffold graph is not a line, forcing the order of the graph")
            scaffold_nodes, inside_nodes, node_order, bo, bubble_count = force_graph_order(
                graph, scaffold_graph, bubbles, artic_points, component_name, inside_nodes, bo_start
            )
            if node_order is not None:
                _mark_nodes_unordered(node_order, unordered_nodes)
            return scaffold_nodes, inside_nodes, node_order, bo, bubble_count

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
    _mark_nodes_unordered(node_order, unordered_nodes)
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
    for comp in components:
        current_tag = ""
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
