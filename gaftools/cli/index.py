"""
Indexing the GAF file for the view functionality.

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
from gaftools.gaf import GAF
from gaftools.gfa import GFA
from gaftools.timer import StageTimer
from gaftools.cli import log_memory_usage

logger = logging.getLogger(__name__)


def run(gaf_path, gfa_path, output=None):
    timers = StageTimer()
    if output is None:
        output = gaf_path + ".gvi"

    # Detecting if GAF has stable or unstable coordinate
    gaf = GAF(gaf_path)
    stable = None
    # checking format in the first 10 lines.
    for i, gaf_line in enumerate(gaf.read_file()):
        if i == 10:
            break
        if i == 0:
            stable = gaf_line.detect_path_format()
        assert stable == gaf_line.detect_path_format()
    gaf.close()

    nodes = {}
    reference = defaultdict(lambda: [])
    ref_contig = []

    gfa_file = GFA(graph_file=gfa_path, low_memory=True)
    contigs = list(gfa_file.contigs.keys())
    ref_contig = [contig for contig in gfa_file.contigs if gfa_file.contigs[contig] == 0]
    nodes = gfa_file.nodes
    if stable:
        for contig in contigs:
            path = gfa_file.get_path(contig)
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
            val = mapping.rstrip().split("\t")
        except TypeError:
            val = mapping.decode("utf-8").rstrip().split("\t")

        if stable:
            with timers("convert_coord"):
                alignment = convert_coord(val, reference)
        else:
            alignment = list(re.split(">|<", val[5]))[1:]
        for a in alignment:
            try:
                out_dict[
                    (
                        nodes[a].id,
                        nodes[a].tags["SN"][1],
                        int(nodes[a].tags["SO"][1]),
                        int(nodes[a].tags["SO"][1]) + int(nodes[a].tags["LN"][1]),
                    )
                ].append(offset)
            except KeyError:
                out_dict[
                    (
                        nodes[a].id,
                        nodes[a].tags["SN"][1],
                        int(nodes[a].tags["SO"][1]),
                        int(nodes[a].tags["SO"][1]) + int(nodes[a].tags["LN"][1]),
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
    logger.info("Time to sort gfa:                            %9.2f s", timers.elapsed("sort_gfa"))
    logger.info(
        "Time to store contig info:                   %9.2f s", timers.elapsed("store_contig_info")
    )
    logger.info("Total time:                                  %9.2f s", total_time)


def convert_coord(line, ref):
    unstable_coord = []
    gaf_contigs = list(filter(None, re.split("(>)|(<)", line[5])))
    for nd in gaf_contigs:
        if nd == ">" or nd == "<":
            continue
        if ":" in nd and "-" in nd:
            tmp = nd.rstrip().split(":")
            query_contig_name = tmp[0]
            (query_start, query_end) = tmp[1].rstrip().split("-")
        else:
            query_start = line[7]
            query_end = line[8]
            query_contig_name = nd

        """Find the matching nodes from the reference genome here"""
        start, end = utils.search_intervals(
            ref[query_contig_name], int(query_start), int(query_end), 0, len(ref[query_contig_name])
        )

        for node in ref[query_contig_name][start : end + 1]:
            cases = -1
            if (
                int(node.tags["SO"][1])
                <= int(query_start)
                < int(node.tags["SO"][1]) + int(node.tags["LN"][1])
            ):
                cases = 1
            elif (
                int(node.tags["SO"][1])
                < int(query_end)
                <= int(node.tags["SO"][1]) + int(node.tags["LN"][1])
            ):
                cases = 2
            elif (
                int(query_start)
                < int(node.tags["SO"][1])
                < int(node.tags["SO"][1]) + int(node.tags["LN"][1])
                < int(query_end)
            ):
                cases = 3

            if cases != -1:
                unstable_coord.append(node.id)

    return unstable_coord


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_path', metavar='GAF',
        help='Input GAF file (can be bgzip-compressed)')
    arg('gfa_path', metavar='rGFA',
        help='Reference rGFA file')
    arg('-o', '--output', default=None,
        help='Path to the output Indexed GAF file. If omitted, use <GAF File>.gvi')


# fmt: on
def validate(args, parser):
    return True


def main(args):
    run(**vars(args))
