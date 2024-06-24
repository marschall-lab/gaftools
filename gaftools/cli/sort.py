"""
Sorting GAF alignments using BO and NO tags of the corresponding graph

The script uses the BO and NO tags defined by the order_gfa command. Using the bubble ordering done, the alignments are sorted.

The index is dictionary created using pickle library which contains the reference contig names as keys and the offset the alignments begin and end..

It adds some tags into the sorted gaf file. The tags are:
    1. bo:i: - This is the BO tags from the GFA carried forward. The sorting uses a BO tag for each alignment, the BO tag for the start node (or end node based on overall orientation) of the alignment path.
    2. sn:Z: - The name of the reference contig it mapped to.
    3. iv:i: - 1 if the alignment path has an inversion. 0 otherwise.
"""

import sys
import logging
import functools
import re
import resource
import pickle as pkl
from pysam import libcbgzf
from collections import defaultdict, namedtuple

import gaftools.utils as utils
from gaftools.timer import StageTimer
from gaftools.gfa import GFA

logger = logging.getLogger(__name__)
timers = StageTimer()


def run_sort(gfa, gaf, outgaf=None, outind=None, bgzip=False):
    if outgaf is None:
        writer = sys.stdout
        index_file = None
    else:
        if bgzip:
            writer = libcbgzf.BGZFile(outgaf, "wb")
        else:
            writer = open(outgaf, "w")
        if outind:
            index_file = outind
        else:
            index_file = outgaf + ".gsi"
    index_dict = defaultdict(lambda: [None, None])
    with timers("read_gfa"):
        gfa_file = GFA(graph_file=gfa, low_memory=True)
    with timers("total_sort"):
        sort(gaf, gfa_file.nodes, writer, index_dict, index_file)
    writer.close()

    memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    logger.info("\nMemory Information")
    logger.info("  Maximum memory usage:              %.3f GB", memory_kb / 1e6)
    logger.info("\nTime Summary:")
    logger.info("  Time to parse GFA file:            %.3f seconds" % (timers.elapsed("read_gfa")))
    logger.info(
        "  Total time to sort GAF file:       %.3f seconds" % (timers.elapsed("total_sort"))
    )
    logger.info("    Time to parse GAF file:          %.3f seconds" % (timers.elapsed("read_gaf")))
    logger.info("    Time to sort GAF file:           %.3f seconds" % (timers.elapsed("sort_gaf")))
    logger.info("    Time to write GAF file:          %.3f seconds" % (timers.elapsed("write_gaf")))


def sort(gaf, nodes, writer, index_dict, index_file):
    logger.info("Parsing GAF file and sorting it")
    if utils.is_file_gzipped(gaf):
        reader = libcbgzf.BGZFile(gaf, "rb")
    else:
        reader = open(gaf, "r")

    Alignment = namedtuple("Alignment", ["offset", "BO", "NO", "start", "inv", "sn"])
    gaf_alignments = []
    count_inverse = 0
    # First pass: Store all the alignment lines as minimally. Just storing line offset and alignment string.
    with timers("read_gaf"):
        while True:
            offset = reader.tell()
            line = reader.readline()
            if not line:
                break
            try:
                line = line.rstrip().split("\t")
            except TypeError:
                line = line.decode("utf8").rstrip().split("\t")
            bo, no, start, inv, sn = process_alignment(line, nodes, offset)
            if inv == 1:
                count_inverse += 1
            gaf_alignments.append(
                Alignment(offset=offset, BO=bo, NO=no, start=start, inv=inv, sn=sn)
            )

    logger.info("\tNumber of alignments with inversions: %d" % (count_inverse))
    # Sorting the alignments based on BO and NO tag
    with timers("sort_gaf"):
        logger.info("\tSorting the alignments...")
        gaf_alignments.sort(key=functools.cmp_to_key(compare_gaf))

    # Writing the sorted file
    with timers("write_gaf"):
        logger.debug("Writing Output File...")
        for alignment in gaf_alignments:
            off = alignment.offset
            reader.seek(off)
            line = reader.readline()
            if isinstance(line, bytes):
                line = line.decode("utf-8").rstrip()
            elif isinstance(line, str):
                line = line.rstrip()
            else:
                raise RuntimeError("GAF alignments not in string or byte format.")
            line += "\tbo:i:%d\tsn:Z:%s\tiv:i:%d\n" % (alignment.BO, alignment.sn, alignment.inv)
            if index_file is not None:
                out_off = writer.tell()
                if index_dict[alignment.sn][0] is None:
                    index_dict[alignment.sn][0] = out_off
                    index_dict[alignment.sn][1] = out_off
                else:
                    index_dict[alignment.sn][1] = out_off
            write_to_file(line, writer)
    if index_file is not None:
        index_dict.pop("unknown")
        index_dict = dict(index_dict)
        with open(index_file, "wb") as ind:
            pkl.dump(index_dict, ind)
    reader.close()


def process_alignment(line, nodes, offset):
    path = list(filter(None, re.split("(>)|(<)", line[5])))
    orient = None
    # If there is no scaffold node present in the alignment, then assuming that it is in the correct orientation.
    # TODO: Need to find a way to deal with such alignments.
    bo = None
    no = None
    start = None
    orient_list = []
    inv = 0  # Whether the alignment has an inversion
    sn = None
    for n in path:
        if n in [">", "<"]:
            orient = n
            continue

        sn_tag = nodes[n].tags["SN"][1]
        bo_tag = int(nodes[n].tags["BO"][1])
        no_tag = int(nodes[n].tags["NO"][1])
        sr_tag = int(nodes[n].tags["SR"][1])

        # Finding the chromosome where the alignment is
        if sn is None and sr_tag == 0:
            sn = sn_tag
        elif sr_tag == 0:
            assert sn == sn_tag

        if bo_tag == -1 or no_tag == -1:
            logger.debug("[ERR]\tOF:i:%d\tND:Z:%s" % (offset, n))
            continue
        # Skipping the non-scaffold nodes
        if no_tag != 0:
            continue
        # Only keeping the orientation of the scaffold nodes in the path matching
        orient_list.append(orient)
    # If there are scaffold nodes in both direction, this indicates an inversion
    if orient_list.count(">") != 0 and orient_list.count("<") != 0:
        inv = 1
    # TODO: What to do when the start node that we have is an untagged node?
    # One argument is to just sort them to the end and completely ignore them since they dont belong to any vcf bubble.
    # That works for my tools but not ideal for general purpose.
    if orient_list.count(">") < orient_list.count("<"):
        l = int(line[6])
        e = int(line[8])
        start = l - e
        n = path[-1]
        bo = int(nodes[n].tags["BO"][1])
        no = int(nodes[n].tags["NO"][1])
    else:
        start = int(line[7])
        n = path[1]
        bo = int(nodes[n].tags["BO"][1])
        no = int(nodes[n].tags["NO"][1])
    if sn is None:
        sn = "unknown"
    return bo, no, start, inv, sn


def compare_gaf(al1, al2):
    # Since we are sorting in ascending order, al1 is above al2 if it comes before al2.
    # Have to consider the case when the start node is untagged and has BO and NO has -1. These should be sorted to the end and not the start
    if al1.BO == -1 and al2.BO == -1:
        if al1.offset < al2.offset:
            return -1
        else:
            return 1
    elif al1.BO == -1:
        return 1
    elif al2.BO == -1:
        return 1

    # Comparing BO tags
    if al1.BO < al2.BO:
        return -1
    if al1.BO > al2.BO:
        return 1

    # Comparing NO tags
    if al1.NO < al2.NO:
        return -1
    if al1.BO > al2.BO:
        return 1

    # Comparing start position in node
    if al1.start < al2.start:
        return -1
    if al1.start > al2.start:
        return 1

    # This will be execulted only when two reads start at the exact same position at the same node. Then we use offset values. So the read which comes first in the GAF file will come higher up.
    if al1.offset < al2.offset:
        return -1
    if al1.offset > al2.offset:
        return 1


def write_to_file(line, writer):
    try:
        writer.write(line)
    except TypeError:
        writer.write(str.encode(line))


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg("gaf", metavar='GAF',
        help="Input GAF File (can be bgzip-compressed)")
    arg("gfa", metavar='GFA',
        help="GFA file with the sort keys (BO and NO tagged). This is done with gaftools order_gfa")
    arg("--outgaf", default=None,
        help="Output GAF File path (Default: sys.stdout)")
    arg("--outind", default=None,
        help="Output Index File path for the GAF file. "
        "When --outgaf is not given, no index is created. "
        "If it is given and --outind is not specified, it will have same file name with .gsi extension)")
    arg("--bgzip", action='store_true',
        help="Flag to bgzip the output. Can only be given with --outgaf.")


# fmt: on
def validate(args, parser):
    if args.bgzip and not args.outgaf:
        parser.error(
            "--bgzip flag has been specified but not output path has been defined. Please define the output path."
        )
    if args.outind and not args.outgaf:
        parser.error("index path specified but no output gaf path. Please provide an output path.")


def main(args):
    run_sort(**vars(args))
