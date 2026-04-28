"""
Core scripts for converting coordinates.
"""

import logging
import re

import gaftools.utils as utils
from gaftools.gaf import GAF
from gaftools.errors import CommandLineError, IncorrectGfaFormatError

logger = logging.getLogger(__name__)


# Stable Node representation of the Nodes.
# Used in the coversion from unstable to stable coordinate system.
class StableNode:
    def __init__(self, contig_id, start, end):
        self.contig_id = contig_id
        self.start = start
        self.end = end

    def to_string(self, orient):
        return "%s%s:%d-%d" % (orient, self.contig_id, self.start, self.end)


def merge_nodes(node1, node2, orient1, orient2):
    if (not isinstance(node1, StableNode)) or (not isinstance(node2, StableNode)):
        return False
    if (node1.contig_id != node2.contig_id) or (orient1 != orient2):
        return False
    if (orient1 == ">") and (node1.end != node2.start):
        return False
    if (orient1 == "<") and (node1.start != node2.end):
        return False
    if orient1 == "<":
        node = StableNode(node1.contig_id, node2.start, node1.end)
    else:
        node = StableNode(node1.contig_id, node1.start, node2.end)
    return [node, orient1]


def stable_to_unstable(gaf_path, reference):
    """
    This function converts a GAF file (mappings to a pangenome graph) into unstable coordinate.
    It does not expect sorted input however it strictly assumes that SO, LN and SN tags are available in the rGFA.
    """
    gaf_input = GAF(gaf_path)
    for gaf_line in gaf_input.read_file():
        yield to_unstable(gaf_line, reference)
    gaf_input.close()


def unstable_to_stable(gaf_path, nodes, ref_contig, contig_len):
    """This function converts a gaf file (mappings to a pangenome graph). It does not need sorted
    input however it strictly assumes that SO, LN and SN tags are available...
    """
    gaf_input = GAF(gaf_path)
    for gaf_line in gaf_input.read_file():
        yield to_stable(gaf_line, nodes, ref_contig, contig_len)
    gaf_input.close()


# separate function for converting lines to unstable coordinate
def to_unstable(gaf_line, reference):
    gaf_contigs = list(filter(None, re.split("(>)|(<)", gaf_line.path)))
    assert len(gaf_contigs) >= 1

    unstable_coord = ""
    orient = None

    new_start = -1
    new_total = 0
    query_contig_name = None
    query_start = None
    query_end = None
    split_contig = None
    for nd in gaf_contigs:
        is_unstable_node = False
        if nd == ">" or nd == "<":
            orient = nd
            continue
        if ":" in nd and "-" in nd:
            tmp = nd.rstrip().split(":")
            query_contig_name = tmp[0]
            (query_start, query_end) = tmp[1].rstrip().split("-")
            split_contig = True
        else:
            if len(gaf_contigs) == 1:
                query_start = gaf_line.path_start
                query_end = gaf_line.path_end
                query_contig_name = nd
                split_contig = False
            else:
                # we found a vertex node.
                is_unstable_node = True
        if not orient:
            if gaf_line.strand == "+":
                orient = ">"
            else:
                orient = "<"
        if is_unstable_node:
            unstable_coord += orient + nd
            split_contig = True
            continue
        """Find the matching nodes from the reference genome here"""
        if reference[query_contig_name] == []:
            raise CommandLineError(
                f"Contig name {query_contig_name} in the GAF file was not found in the GFA. "
                f"Check if the GFA appropriate tags like SN, SR and SO."
            )
        node_indices = get_nodes_from_region(
            [query_contig_name, query_start, query_end], reference[query_contig_name]
        )
        nodes_tmp = []
        for i in node_indices:
            node = reference[query_contig_name][i]
            s = int(node.tags["SO"][1])
            e = int(node.tags["SO"][1]) + int(node.seq_len)
            cases = -1
            if s <= int(query_start) < e:
                cases = 1
                if new_start == -1:
                    if split_contig:
                        new_start = int(query_start)
                    else:
                        new_start = int(query_start) - s
            elif s < int(query_end) <= e:
                cases = 2
            elif int(query_start) < s < e < int(query_end):
                cases = 3

            if cases != -1:
                nodes_tmp.append(node.id)
                new_total += e - s

        if orient == "<":
            for i in reversed(nodes_tmp):
                unstable_coord += orient + i
        else:
            for i in nodes_tmp:
                unstable_coord += orient + i

    if gaf_line.strand == "-":
        if split_contig:
            new_total = gaf_line.path_length
            new_start = new_total - gaf_line.path_end
            new_end = new_total - gaf_line.path_start
        else:
            new_end = new_total - new_start
            new_start = new_end - (gaf_line.path_end - gaf_line.path_start)
    else:
        if split_contig:
            new_total = gaf_line.path_length
            new_start = gaf_line.path_start
            new_end = gaf_line.path_end
        else:
            new_end = new_start + (gaf_line.path_end - gaf_line.path_start)

    new_line = "%s\t%s\t%s\t%s\t+\t%s\t%d\t%d\t%d\t%d\t%d\t%d" % (
        gaf_line.query_name,
        gaf_line.query_length,
        gaf_line.query_start,
        gaf_line.query_end,
        unstable_coord,
        new_total,
        new_start,
        new_end,
        gaf_line.residue_matches,
        gaf_line.alignment_block_length,
        gaf_line.mapping_quality,
    )

    # Add cigar in reverse
    if gaf_line.strand == "-":
        new_cigar = utils.reverse_cigar(gaf_line.cigar)
        gaf_line.tags["cg:Z:"] = new_cigar

    for k in gaf_line.tags.keys():
        new_line += "\t%s%s" % (k, gaf_line.tags[k])

    return new_line


def to_stable(gaf_line, tagged_nodes, ref_contig, contig_len):
    reverse_flag = False
    new_total = None
    new_start = None
    if not (">" in gaf_line.path or "<" in gaf_line.path) or ":" in gaf_line.path:
        raise CommandLineError(
            f"Read {gaf_line.query_name} has path {gaf_line.path} which does not follow unstable coordiante format. The file has to be in unstable coordinate."
        )
    gaf_nodes = list(filter(None, re.split("(>)|(<)", gaf_line.path)))
    node_list = []
    stable_coord = ""
    orient = None
    new_line = ""

    for nd in gaf_nodes:
        if nd == ">" or nd == "<":
            orient = nd
            continue
        if not orient:
            orient = ">"
        if nd in tagged_nodes:
            node_list.append([tagged_nodes[nd], orient])
        else:
            node_list.append([nd, orient])
    out_node = [node_list[0]]

    for i in range(len(node_list) - 1):
        n1 = out_node[-1][0]
        o1 = out_node[-1][1]
        n2 = node_list[i + 1][0]
        o2 = node_list[i + 1][1]
        node_merge = merge_nodes(n1, n2, o1, o2)
        if node_merge is False:
            stable_coord += n1.to_string(o1) if isinstance(n1, StableNode) else f"{o1}{n1}"
            out_node.append([n2, o2])
        else:
            out_node[-1] = node_merge

    if (
        len(out_node) == 1
        and isinstance(out_node[0][0], StableNode)
        and out_node[0][0].contig_id in ref_contig
    ):
        if out_node[0][1] == "<":
            reverse_flag = True
            gaf_line.strand = "-"
            new_start = out_node[0][0].start + gaf_line.path_length - gaf_line.path_end
        else:
            new_start = out_node[0][0].start + gaf_line.path_start

        stable_coord = out_node[0][0].contig_id
        new_total = contig_len[stable_coord]
    else:
        node = out_node[-1]
        stable_coord += (
            node[0].to_string(node[1]) if isinstance(node[0], StableNode) else f"{node[1]}{node[0]}"
        )
        new_start = gaf_line.path_start
        new_total = gaf_line.path_length

    new_line += "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" % (
        gaf_line.query_name,
        gaf_line.query_length,
        gaf_line.query_start,
        gaf_line.query_end,
        gaf_line.strand,
        stable_coord,
        new_total,
        new_start,
        new_start + gaf_line.path_end - gaf_line.path_start,
        gaf_line.residue_matches,
        gaf_line.alignment_block_length,
        gaf_line.mapping_quality,
    )

    # Add cigar in reverse
    if reverse_flag:
        new_cigar = utils.reverse_cigar(gaf_line.cigar)
        gaf_line.tags["cg:Z:"] = new_cigar

    # adding tags in the sequence it was found
    for k in gaf_line.tags.keys():
        new_line += "\t%s%s" % (k, gaf_line.tags[k])

    return new_line


def get_nodes_from_region(region, tagged_contig_nodes):
    """
    Finds the unstable node indices in tagged_contig_nodes from the region and the nodes belonging to the contig.

    region: [contig, start, end] encodes the region information
    tagged_contig_nodes: Iterable[gfa.Node]
    """

    node_list = []
    for node in tagged_contig_nodes:
        if len(node_list) > 0:
            assert node_list[-1][2] < int(node.tags["SO"][1])
        node_list.append(
            [
                node.id,
                node.tags["SN"][1],
                int(node.tags["SO"][1]),
                int(node.tags["SO"][1]) + node.seq_len,
            ]
        )

    unstable_node_indices = search(region, node_list)

    # checking if the unstable nodes cover the whole region or there are gaps.
    # cannot deal with gaps.
    for i in range(1, len(unstable_node_indices)):
        prev_unstable_node = node_list[unstable_node_indices[i] - 1]
        curr_unstable_node = node_list[unstable_node_indices[i]]
        if prev_unstable_node[3] != curr_unstable_node[2]:
            raise IncorrectGfaFormatError(
                f"Tried to convert region {region[0]}:{region[1]}-{region[2]}. "
                f"Found discontinuous nodes {', '.join([node_list[x][0] for x in unstable_node_indices])}. Region is not fully annotated in the GFA."
            )

    return unstable_node_indices


def search(region, node_list):
    """
    Find the unstable node indices in node_list from region information

    region: [contig, start, end] encodes the region information
    node_list: sorted list of tagged nodes of format Iterable[node_name, contig, start, end] (sort index is start)

    returns the nodes that belong to the region
    """

    if region[1] is None:
        # Only contig is given
        assert region[2] is None
        return node_list
    s = 0
    pos = 0
    e = len(node_list) - 1
    q_s = int(region[1])
    q_e = int(region[2])
    while s != e:
        m = int((s + e) / 2)
        if (q_s >= node_list[m][2]) and (q_s < node_list[m][3]):
            pos = m
            break
        elif q_s >= node_list[m][3]:
            s = m + 1
        else:
            e = m - 1
        pos = s

    # This if-elif statement will cover some edge cases where the node_list coming from the index has gaps.
    # This can happen when the alignment does not cover a node.
    # In the examples below ------ will denote nodes that are not covered and ++++++ are the nodes that are covered.
    # In the binary search, if the query state does not belong to a +++++ node, then it will report either the node
    #   to the left or right of the absent node.
    if q_s < node_list[pos][2] and q_s < node_list[pos][3]:
        #  Case 1:  ------------+++++++++++++
        #               ^              ^
        #           query start    node reported
        if q_e <= node_list[pos][2]:
            # query region ends before the reported node starts
            return []
        result = [pos]
        pos += 1
        while pos < len(node_list):
            if q_e <= node_list[pos][2]:
                break
            result.append(pos)
            pos += 1
        return result
    elif q_s > node_list[pos][2] and q_s >= node_list[pos][3]:
        #  Case 2:  ++++++++++++---------------
        #               ^              ^
        #           node reported  query start
        if pos == len(node_list) - 1:
            # reported node is the last node
            return []
        pos += 1
        if q_e <= node_list[pos][2]:
            # query region ends before the current node starts
            return []
        result = [pos]
        pos += 1
        while pos < len(node_list):
            if q_e <= node_list[pos][2]:
                break
            result.append(pos)
            pos += 1
        return result

    # now we are back to the normal case
    # -----------+++++++++++++++++---------------
    #                   ^
    #  query starts here and node is reported here
    assert q_s >= node_list[pos][2] and q_s < node_list[pos][3]
    result = [pos]
    pos += 1
    while pos < len(node_list):
        if q_e <= node_list[pos][2]:
            break
        result.append(pos)
        pos += 1

    return result
