"""
Index a GAF file for the view functionality.

This script creates an inverse look-up table where:
    - key: the node information
    - value: the offsets in the GAF file where the node is present
"""

import re
import pickle
from pysam import libcbgzf
from collections import defaultdict
import logging

import gaftools.utils as utils
from gaftools.gfa import GFA
from gaftools.timer import StageTimer
from gaftools.cli import log_memory_usage
from gaftools.conversion import get_nodes_from_region
from gaftools.errors import IncorrectGfaFormatError

logger = logging.getLogger(__name__)


def run(gaf_path, gfa_path, output=None):
    timers = StageTimer()
    if output is None:
        output = gaf_path + ".gvi"

    nodes = {}
    reference = defaultdict(lambda: [])
    ref_contig = []

    gfa_file = GFA(graph_file=gfa_path, low_memory=True)
    contigs = list(gfa_file.contigs.keys())
    ref_contig = [contig for contig in gfa_file.contigs if gfa_file.contigs[contig] == 0]
    for contig in contigs:
        path = gfa_file.get_path(contig, throw_warning=False)
        for node in path:
            reference[contig].append(gfa_file[node])
    nodes = gfa_file.nodes
    del gfa_file

    if utils.is_file_gzipped(gaf_path):
        gaf_file = libcbgzf.BGZFile(gaf_path, "rb")
    else:
        gaf_file = open(gaf_path, "rt")

    out_dict = {}
    offset = 0
    while True:
        offset = gaf_file.tell()
        mapping = gaf_file.readline()
        if not mapping:
            break
        try:
            gaf_line = mapping.rstrip().split("\t")
        except TypeError:
            gaf_line = mapping.decode("utf-8").rstrip().split("\t")

        with timers("convert_coord"):
            unstable_nodes = convert_coord(gaf_line, reference)
        for nd in unstable_nodes:
            try:
                out_dict[
                    (
                        nodes[nd].id,
                        nodes[nd].tags["SN"][1],
                        int(nodes[nd].tags["SO"][1]),
                        int(nodes[nd].tags["SO"][1]) + int(nodes[nd].seq_len),
                    )
                ].append(offset)
            except KeyError:
                out_dict[
                    (
                        nodes[nd].id,
                        nodes[nd].tags["SN"][1],
                        int(nodes[nd].tags["SO"][1]),
                        int(nodes[nd].tags["SO"][1]) + int(nodes[nd].seq_len),
                    )
                ] = [offset]
    out_dict["ref_contig"] = ref_contig

    gaf_file.close()

    with open(output, "wb") as handle:
        with timers("write_file"):
            pickle.dump(out_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info(
        "Time to convert coordinates:                 %9.2f s", timers.elapsed("convert_coord")
    )
    logger.info(
        "Time to write the file:                      %9.2f s", timers.elapsed("write_file")
    )
    logger.info("Time spent on rest:                          %9.2f s", total_time - timers.sum())
    logger.info("Total time:                                  %9.2f s", total_time)


def convert_coord(line, ref):
    has_stable = False
    if (">" not in line[5] and "<" not in line[5]) or (":" in line[5]):
        has_stable = True
    if not has_stable:
        return list(re.split(">|<", line[5]))[1:]

    unstable_coord = []
    path = list(filter(None, re.split("(>)|(<)", line[5])))
    for nd in path:
        is_unstable_node = False
        if nd == ">" or nd == "<":
            continue
        if ":" in nd and "-" in nd:
            tmp = nd.rstrip().split(":")
            query_contig_name = tmp[0]
            (query_start, query_end) = tmp[1].rstrip().split("-")
        else:
            if len(path) == 1:
                query_start = line[7]
                query_end = line[8]
                query_contig_name = nd
            else:
                # we found a vertex node.
                is_unstable_node = True

        if is_unstable_node:
            unstable_coord.append(nd)
            continue
        """Find the matching nodes from the reference genome here"""
        if ref[query_contig_name] == []:
            raise IncorrectGfaFormatError(
                f"Found stable cooridnates for contig {query_contig_name} in the GAF file "
                "but appropriate tags not found in GFA. "
                "Check if the GFA provided is the same as the one used for alignment."
            )
        node_indices = get_nodes_from_region(
            [query_contig_name, query_start, query_end], ref[query_contig_name]
        )
        for i in node_indices:
            node = ref[query_contig_name][i]
            s = int(node.tags["SO"][1])
            e = int(node.tags["SO"][1]) + int(node.seq_len)
            cases = -1
            if s <= int(query_start) < e:
                cases = 1
            elif s < int(query_end) <= e:
                cases = 2
            elif int(query_start) < s < e < int(query_end):
                cases = 3

            if cases != -1:
                unstable_coord.append(node.id)

    return unstable_coord


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_path', metavar='GAF',
        help='GAF file (can be bgzip-compressed)')
    arg('gfa_path', metavar='rGFA',
        help='rGFA file (can be bgzip-compressed)')
    arg('-o', '--output', default=None,
        help='Output GAF View Index (GVI) file. Default: <GAF file>.gvi')
# fmt: on


def main(args):
    run(**vars(args))
