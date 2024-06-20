"""
This uses the BO and NO tags defined by Samarendra Pani through processing the bubbles in the VCF produced by Glenn Hickey for the MC graph. The script to add these BO and NO tags are in add-bubble-info.py in this directory.

Here since the bubbles are already defined, there is no need to process the entire alignment to find orientation of alignment since inversions are already taken into consideration.
"""

import sys
import argparse
import logging
from collections import defaultdict, namedtuple
import functools
import gzip
import re
from pysam import libcbgzf
import resource
from whatshap.timer import StageTimer

logger = logging.getLogger(__name__)
timers = StageTimer()


def setup_logging(debug):
    handler = logging.StreamHandler()
    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if debug else logging.INFO)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    parser.add_argument(
        "--gfa", required=True, help="GFA file with the sort keys (BO and NO tagged)"
    )
    parser.add_argument("--gaf", required=True, help="Input GAF File")
    parser.add_argument("--output", default=None, help="Output GAF File path (Default: sys.stdout)")
    parser.add_argument(
        "--bgzip",
        action="store_true",
        help="Flag to bgzip the output. Can only be given with --output.",
    )

    options = parser.parse_args()
    setup_logging(options.debug)

    validate_arguments(options)

    bubble_sort(options)

    memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    logger.info("\nMemory Information")
    logger.info("  Maximum memory usage: %.3f GB", memory_kb / 1e6)
    logger.info("\nTime Summary:")
    logger.info("  Time to parse GFA file: %.3f" % (timers.elapsed("read_gfa")))
    logger.info("  Total time to sort GAF file: %.3f" % (timers.elapsed("total_sort")))
    logger.info("    Time to parse GAF file: %.3f" % (timers.elapsed("read_gaf")))
    logger.info("    Time to sort GAF file: %.3f" % (timers.elapsed("sort_gaf")))
    logger.info("    Time to write GAF file: %.3f" % (timers.elapsed("write_gaf")))


def bubble_sort(options):
    if options.output is None:
        writer = sys.stdout
    else:
        if options.bgzip:
            writer = libcbgzf.BGZFile(options.output, "wb")
        else:
            writer = open(options.output, "w")
    nodes = defaultdict(lambda: [-1, -1])
    with timers("read_gfa"):
        read_gfa(options.gfa, nodes)
    with timers("total_sort"):
        sort(options.gaf, nodes, writer)
    writer.close()


def sort(gaf, nodes, writer):
    logger.info("\n##### Parsing GAF file and sorting it #####")
    if is_file_gzipped(gaf):
        reader = gzip.open(gaf, "rt")
    else:
        reader = open(gaf, "r")

    Alignment = namedtuple("Alignment", ["offset", "BO", "NO", "start"])
    gaf_alignments = []

    # First pass: Store all the alignment lines as minimally. Just storing line offset and alignment string.
    with timers("read_gaf"):
        logger.debug("Processing and Storing GAF Alignments...")
        while True:
            offset = reader.tell()
            line = reader.readline()
            if not line:
                break
            line = line.split("\t")
            bo, no, start = process_alignment(line, nodes, offset)
            gaf_alignments.append(Alignment(offset=offset, BO=bo, NO=no, start=start))

    # Sorting the alignments based on BO and NO tag
    with timers("sort_gaf"):
        logger.debug("Sorting the alignments...")
        gaf_alignments.sort(key=functools.cmp_to_key(compare_gaf))

    # Writing the sorted file
    with timers("write_gaf"):
        logger.debug("Writing Output File...")
        for alignment in gaf_alignments:
            off = alignment.offset
            reader.seek(off)
            line = reader.readline()
            write_to_file(line, writer)

    reader.close()


def read_gfa(gfa, node):
    logger.info("\n##### Parsing GFA file and reading sort key information #####")
    if is_file_gzipped(gfa):
        reader = gzip.open(gfa, "rt")
    else:
        reader = open(gfa, "r")

    total_nodes = 0
    tagged_nodes = 0
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] != "S":
            continue
        total_nodes += 1
        fields = line.split("\t")
        for f in fields:
            if f.startswith("BO:i:"):
                node[fields[1]][0] = int(f[5:])
                tagged_nodes += 1
            if f.startswith("NO:i:"):
                node[fields[1]][1] = int(f[5:])

    logger.info("Total Nodes Processed: %d" % (total_nodes))
    logger.info("Nodes with tags: %d" % (tagged_nodes))
    reader.close()


def process_alignment(line, nodes, offset):
    path = list(filter(None, re.split("(>)|(<)", line[5])))
    orient = None
    # If there is no scaffold node present in the alignment, then assuming that it is in the correct orientation.
    # TODO: Need to find a way to deal with such alignments.
    rv = False
    bo = None
    no = None
    start = None
    for n in path:
        if n in [">", "<"]:
            orient = n
            continue
        if nodes[n][0] == -1 or nodes[n][1] == -1:
            logger.debug("[ERR]\tOF:i:%d\tND:Z:%s" % (offset, n))
            # raise RuntimeError("Found a node which was not in VCF.")
            # These nodes exist.
            continue
        if nodes[n][0] % 2 == 1 and nodes[n][1] == 0:
            # Here we assume that the orientation of the first scaffold node found is the orientation for all the scaffold nodes. (Here we will have to keep in mind that some scaffold nodes are inside the bubbles also and their orientation can be anything there.)
            # TODO: Have to check for these scaffold nodes that can be present inside bubbles (Cannot trust them)
            # TODO: Have given the NO tag of 1 to these scaffold nodes to distinguish them from normal scaffold nodes which are tagged with 0.
            if orient == ">":
                break
            else:
                rv = True
                break
    # TODO: What to do when the start node that we have is an untagged node?
    # One argument is to just sort them to the end and completely ignore them since they dont belong to any vcf bubble.
    # That works for my tools but not ideal for general purpose.
    if rv:
        l = int(line[6])
        e = int(line[8])
        start = l - e
        n = path[-1]
        bo = nodes[n][0]
        no = nodes[n][1]
    else:
        start = int(line[7])
        n = path[1]
        bo = nodes[n][0]
        no = nodes[n][1]
    # print(bo, no, start, sep='\t')
    return bo, no, start


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


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b"\x1f\x8b"


def validate_arguments(options):
    if options.bgzip and not options.output:
        raise RuntimeError(
            "--bgzip flag has been specified but not output path has been defined. Please define the output path."
        )


if __name__ == "__main__":
    main()
