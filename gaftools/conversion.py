"""
Core scripts for converting coordinates.
"""

import logging
import re

import gaftools.utils as utils
from gaftools.gaf import GAF


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
    for nd in gaf_contigs:
        if nd == ">" or nd == "<":
            orient = nd
            continue
        if ":" in nd and "-" in nd:
            tmp = nd.rstrip().split(":")
            query_contig_name = tmp[0]
            (query_start, query_end) = tmp[1].rstrip().split("-")
            split_contig = True
        else:
            query_start = gaf_line.path_start
            query_end = gaf_line.path_end
            query_contig_name = nd
            split_contig = False
        if not orient:
            if gaf_line.strand == "+":
                orient = ">"
            else:
                orient = "<"

        """Find the matching nodes from the reference genome here"""
        start, end = utils.search_intervals(
            reference[query_contig_name],
            int(query_start),
            int(query_end),
            0,
            len(reference[query_contig_name]),
        )
        nodes_tmp = []
        for i in reference[query_contig_name][start : end + 1]:
            s = int(i.tags["SO"][1])
            e = int(i.tags["SO"][1]) + int(i.tags["LN"][1])
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
                nodes_tmp.append(i.id)
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


def to_stable(gaf_line, nodes, ref_contig, contig_len):
    reverse_flag = False
    new_total = None
    new_start = None
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
        node_list.append([nodes[nd], orient])
    out_node = [node_list[0]]

    for i in range(len(node_list) - 1):
        n1 = out_node[-1][0]
        o1 = out_node[-1][1]
        n2 = node_list[i + 1][0]
        o2 = node_list[i + 1][1]
        node_merge = merge_nodes(n1, n2, o1, o2)
        if node_merge is False:
            stable_coord += n1.to_string(o1)
            out_node.append([n2, o2])
        else:
            out_node[-1] = node_merge

    if len(out_node) == 1 and out_node[0][0].contig_id in ref_contig:
        if out_node[0][1] == "<":
            reverse_flag = True
            gaf_line.strand = "-"
            new_start = out_node[0][0].start + gaf_line.path_length - gaf_line.path_end
        else:
            new_start = out_node[0][0].start + gaf_line.path_start

        stable_coord = out_node[0][0].contig_id
        new_total = contig_len[stable_coord]
    else:
        stable_coord += out_node[-1][0].to_string(out_node[-1][1])
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
