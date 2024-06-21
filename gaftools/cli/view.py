"""
View the GAF file (subsetting and converting GAF format available). GAF file should first be index using gaftools index.

The view command allows subsetting the GAF file based on node IDs or regions available.
"""

import logging
import pickle
import os
import sys
from collections import defaultdict

from gaftools.cli import log_memory_usage, CommandLineError
from gaftools.timer import StageTimer
from gaftools.gaf import GAF
from gaftools.gfa import GFA
from gaftools.conversion import (
    StableNode,
    stable_to_unstable,
    unstable_to_stable,
    to_stable,
    to_unstable,
)

logger = logging.getLogger(__name__)


def run(gaf_path, gfa=None, output=None, index=None, nodes=[], regions=[], format=None):
    timers = StageTimer()

    if output is None:
        writer = sys.stdout
    else:
        writer = open(output, "w")
    # Need to detect the format of the input gaf
    gaf = GAF(gaf_path)
    gaf_format = None
    # checking format in the first 10 lines.
    for i, gaf_line in enumerate(gaf.read_file()):
        if i == 10:
            break
        if i == 0:
            gaf_format = gaf_line.detect_path_format()
        assert gaf_format == gaf_line.detect_path_format()
    gaf.close()

    ref_contig = []
    gfa_nodes = None
    contig_len = {}
    # if format is given, prepare some objects for use later
    if format:
        if format == "stable":
            if gaf_format is True:
                raise CommandLineError(
                    "Input GAF already has stable coordinates. Please remove the --format stable option"
                )
            gfa_file = GFA(graph_file=gfa, low_memory=True)
            gfa_nodes = {
                id: StableNode(
                    contig_id=gfa_file[id].tags["SN"][1],
                    start=int(gfa_file[id].tags["SO"][1]),
                    end=int(gfa_file[id].tags["SO"][1]) + int(gfa_file[id].tags["LN"][1]),
                )
                for id in gfa_file.nodes
            }
            ref_contig = [contig for contig in gfa_file.contigs if gfa_file.contigs[contig] == 0]
            for contig in gfa_file.contigs:
                contig_len[contig] = gfa_file.get_contig_length(contig)
            print(contig_len)
            del gfa_file
        else:
            assert format == "unstable"
            if gaf_format is False:
                raise CommandLineError(
                    "Input GAF already has unstable coordinates. Please remove the --format unstable option"
                )
            reference = defaultdict(lambda: [])
            # Assuming that the gfa is sorted using the order_gfa function.
            gfa_file = GFA(graph_file=gfa, low_memory=True)
            contigs = list(gfa_file.contigs.keys())
            for contig in contigs:
                path = gfa_file.get_path(contig)
                for node in path:
                    reference[contig].append(gfa_file[node])
            del gfa_file

    # now find out what lines to view and how to view
    if len(nodes) != 0 or len(regions) != 0:
        if index is None:
            index = gaf_path + ".gvi"
            if not os.path.exists(index):
                raise CommandLineError(
                    "No index found. Please provide the path to the index or create one with gaftools index."
                )

        ind = None
        with open(index, "rb") as tmp:
            ind = pickle.load(tmp)

        ind_key = sorted(list(ind.keys()), key=lambda x: (x[1], x[2]))
        ind_dict = {}
        for i in ind_key:
            ind_dict[i[0]] = i

        if regions:
            assert nodes == []
            nodes = get_unstable(regions, ind)
        offsets = ind[ind_dict[nodes[0]]]
        for nd in nodes[1:]:
            # extracting all the lines that touches at least one of the nodes
            offsets = list(set(offsets) | set(ind[ind_dict[nd]]))
        offsets.sort()
        if len(offsets) == 0:
            raise CommandLineError("No alignments found for the given nodes/regions")
        gaf = GAF(gaf_path)
        # if format specified, have to make the changes.
        if format:
            if format == "stable":
                for ofs in offsets:
                    line = gaf.read_line(ofs)
                    print(to_stable(line, gfa_nodes, ref_contig, contig_len), file=writer)
            else:
                assert format == "unstable"
                for ofs in offsets:
                    line = gaf.read_line(ofs)
                    print(to_unstable(line, reference), file=writer)
        # if no format given, then just print the selected lines
        else:
            for ofs in offsets:
                line = gaf.read_line(ofs)
                print(line, file=writer)
        gaf.close()
    else:
        # No nodes or regions indicates the entire file will be viewed
        # converting the format of the entire file
        if format:
            if format == "stable":
                for line in unstable_to_stable(gaf_path, gfa_nodes, ref_contig, contig_len):
                    print(line, file=writer)
            else:
                for line in stable_to_unstable(gaf_path, reference):
                    print(line, file=writer)
        else:
            # No format also given. So just need to print the file.
            gaf = GAF(gaf_path)
            for line in gaf.file:
                if gaf.gz_flag:
                    print(line.decode("utf-8").rstrip(), file=writer)
                else:
                    print(line.rstrip(), file=writer)

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


def get_unstable(regions, index):
    """Takes the regions and returns the node IDs"""

    contig = [x.split(":")[0] for x in regions]
    node_dict = {}
    start = [x.split(":")[1].split("-")[0] for x in regions]
    end = [x.split(":")[1].split("-")[-1] for x in regions]

    result = []
    for n, c in enumerate(contig):
        try:
            node_list = node_dict[c]
        except KeyError:
            node_list = list(filter(lambda x: (x[1] == c), list(index.keys())))
            node_list.sort(key=lambda x: x[2])
            node_dict[c] = node_list

        node = search([contig[n], start[n], end[n]], node_list)
        if len(node) > 1:
            logger.info("INFO: Region %s spans multiple nodes.\nThe nodes are:" % (node[n]))
            for n in node:
                logger.info("INFO: %s\t%s\t%d\t%d" % (n[0], n[1], n[2], n[3]))

        result.append(node[0][0])

    return result


def search(node, node_list):
    """Find the unstable node id from the region"""

    s = 0
    pos = 0
    e = len(node_list) - 1
    q_s = int(node[1])
    q_e = int(node[2])
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
    # if there is only one node for the entire contig (case for non-reference nodes)
    # then the above loop is not executed and we extract the only node with pos=0
    result = [node_list[pos]]
    while True:
        if q_e < node_list[pos][3]:
            break
        pos += 1
        result.append(node_list[pos])

    return result


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_path', metavar='GAF',
        help='Input GAF file (can be bgzip-compressed)')
    arg('-g', '--gfa', dest='gfa', metavar='GFA', default=None,
        help='Input GFA file (can be gzip-compressed). Required when converting from one coordinate system to another.')
    arg('-o', '--output', dest='output', metavar='OUTPUT', default=None,
        help='Output file. Default is stdout.')
    arg('-i', '--index', default=None,
        help='Path to GAF Index file. This index is created using gaftools index. '
        'If path is not provided, it is assumed to be in the same directory as GAF file with the same name and .gvi extension (default location of the index script)')
    arg('-n', '--node', dest='nodes', metavar='NODE', default=[], action='append',
        help='Nodes to search. '
        'Multiple can be provided (Eg. gaftools view .... -n s1 -n s2 -n s3 .....).')
    arg('-r', '--region', dest='regions', metavar='REGION', default=[], action='append',
        help='Regions to search. '
        'Multiple can be provided (Eg. gaftools view .... -r chr1:10-20 -r chr1:50-60 .....).')
    arg('-f', '--format', dest='format', metavar='FORMAT',
        help='format of output path (unstable | stable)')

# fmt: on


def validate(args, parser):
    if args.format and (args.format not in ["unstable", "stable"]):
        parser.error("--format only accepts unstable or stable as input.")
    if args.nodes and args.regions:
        parser.error("provide either of the --regions and --nodes options and not both.")
    if args.format and not args.gfa:
        parser.error("GFA file has to be provided along with --format.")


def main(args):
    run(**vars(args))
