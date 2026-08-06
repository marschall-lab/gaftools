"""
Sort the GAF alignments using BO and NO tags of the corresponding graph.
"""

import logging
import re
import pickle as pkl
from pysam import libcbgzf
from collections import defaultdict
import tempfile
import heapq
import os
import bisect

from gaftools.utils import is_file_gzipped, FileWriter
from gaftools.timer import StageTimer
from gaftools.gfa import GFA
from gaftools.cli import log_memory_usage
from gaftools.errors import IncorrectGafFormatError

logger = logging.getLogger(__name__)
timers = StageTimer()
PATH_SPLIT_REGEX = re.compile(r"(>)|(<)")

"""
The script uses the BO and NO tags defined by the order_gfa command. Using the bubble ordering done, the alignments are sorted.

The index is dictionary created using pickle library which contains the reference contig names as keys and the offset the alignments begin and end.

It adds some tags into the sorted GAF file. The tags are:
    1. bo:i: - This is the BO tags from the GFA carried forward. The sorting uses a BO tag for each alignment, the BO tag for the start node (or end node based on overall orientation) of the alignment path.
    2. sn:Z: - The name of the reference contig it mapped to.
    3. iv:i: - 1 if the alignment path has an inversion. 0 otherwise.
"""


def run_sort(gfa, gaf, outgaf=None, outind=None):
    index_file = None
    if outgaf:
        if outind:
            index_file = outind
        else:
            index_file = outgaf + ".gsi"
    index_dict = defaultdict(lambda: [None, None])
    with timers("read_gfa"):
        gfa_file = GFA(graph_file=gfa, low_memory=True)
    sort(gaf, gfa_file.nodes, outgaf, index_dict, index_file)

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    # fmt: off
    logger.info("Time to parse GFA file:                      %9.2f s" % (timers.elapsed("read_gfa")))
    logger.info("Time to parse GAF file:                      %9.2f s" % (timers.elapsed("read_gaf")))
    logger.info("Time to sort GAF file:                       %9.2f s" % (timers.elapsed("sort_gaf")))
    logger.info("Time to write GAF file:                      %9.2f s" % (timers.elapsed("write_gaf")))
    logger.info("Time spent on rest:                          %9.2f s", total_time - timers.sum())
    logger.info("Total time:                                  %9.2f s", total_time)
    # fmt: on


"""Sorts and flushes all current in-memory buckets to disk."""


def flush_buffer(buf, counts, temp_directory):
    for sn, alignments in buf.items():
        if not alignments:
            continue
        alignments.sort(key=lambda x: x[0])
        chunk_id = counts[sn]
        filepath = os.path.join(temp_directory.name, f"{sn}_{chunk_id}.tmp")
        logger.info(f"Writing sorted temp file {filepath}...")
        with open(filepath, "w") as out:
            for key_tuple, final_line in alignments:
                out.write(
                    f"{key_tuple[0]}\t{key_tuple[1]}\t{key_tuple[2]}\t{key_tuple[3]}\t{key_tuple[4]}\t{final_line}"
                )
        counts[sn] += 1
    buf.clear()


def sort(gaf, nodes, outgaf, index_dict, index_file):
    # Fast binary search to find which chromosome a BO belongs to
    def get_bucket_sn(bo_val):
        if bo_val == -1 or not bo_keys:
            return "unknown"
        # Find the insertion point
        idx = bisect.bisect_right(bo_keys, bo_val) - 1
        if idx < 0:
            return "unknown"  # BO is smaller than the smallest known chromosome BO
        return bo_boundaries[idx][1]

    # stream a chunk file
    def stream_chunk(filepath):
        with open(filepath, "r") as f:
            for line in f:
                k0, k1, k2, k3, k4, final_string = line.split("\t", 5)
                yield ((int(k0), int(k1), int(k2), int(k3), int(k4)), final_string)

    logger.info("Parsing GAF file and sorting it")

    CHUNK_LIMIT = 5_000_000
    total_buffered = 0
    # buffer[sn] = list of tuples
    buffer = defaultdict(list)
    # tracks how many chunks have been written per chromosome
    chunk_counts = defaultdict(int)
    temp_dir = tempfile.TemporaryDirectory()

    if is_file_gzipped(gaf):
        reader = libcbgzf.BGZFile(gaf, "rb")
    else:
        reader = open(gaf, "r")

    try:
        # These intervals will be used to bucket the lines based on chromosome
        logger.info("Extracting BO intervals for chromosomes from GFA...")
        chr_min_bo = {}
        for _, node in nodes.items():
            n_sn = node.tags.get("SN", [None, None])[1]
            n_sr = int(node.tags.get("SR", [None, None])[1])
            n_bo = int(node.tags.get("BO", [None, -1])[1])
            if n_sr == 0 and n_sn and n_bo != -1:
                if n_sn not in chr_min_bo or n_bo < chr_min_bo[n_sn]:
                    chr_min_bo[n_sn] = n_bo
        # Create sorted list of boundaries: [(0, 'chr1'), (1500, 'chr2'), ...]
        bo_boundaries = sorted([(bo, sn) for sn, bo in chr_min_bo.items()])
        bo_keys = [b[0] for b in bo_boundaries]
        ordered_chromosomes = [sn for _, sn in bo_boundaries]

        count_inverse = 0
        # First pass: Using a overall chunk size of 5 million, create separate tmp files with sorted chromosomal alignments.
        with timers("read_gaf"):
            while True:
                offset = reader.tell()
                raw_line = reader.readline()
                if not raw_line:
                    break
                if type(raw_line) is bytes:  # noqa: E721
                    raw_line = raw_line.decode("utf8")
                raw_line = raw_line.rstrip()
                line_arr = raw_line.split("\t", 9)
                # checking for stable coordinates
                if not ((">" in line_arr[5]) or ("<" in line_arr[5])):
                    raise IncorrectGafFormatError(
                        "Detected stable coordinates in GAF. Sorting requires unstable coordinates. "
                        "Please convert the GAF with gaftools view."
                    )
                bo, no, start, inv, sn = process_alignment(line_arr, nodes, offset)
                if inv == 1:
                    count_inverse += 1
                sort_key = (1, 0, 0, 0, offset) if bo == -1 else (0, bo, no, start, offset)
                final_line = f"{raw_line}\tbo:i:{bo}\tsn:Z:{sn}\tiv:i:{inv}\n"

                # ROUTE based on topological BO interval
                bucket_sn = get_bucket_sn(bo)
                buffer[bucket_sn].append((sort_key, final_line))
                total_buffered += 1
                if total_buffered >= CHUNK_LIMIT:
                    flush_buffer(buffer, chunk_counts, temp_dir)
                    total_buffered = 0

            if total_buffered > 0:
                flush_buffer(buffer, chunk_counts, temp_dir)

        logger.info("Number of alignments with inversions: %d" % (count_inverse))

        # Sort chromosomes by their starting BO tag
        ordered_chromosomes = sorted(chr_min_bo.keys(), key=lambda k: chr_min_bo[k])

        # Handle the "unknown" SN variants by throwing them at the end
        if "unknown" in chunk_counts and "unknown" not in ordered_chromosomes:
            ordered_chromosomes.append("unknown")

        # Writing the sorted file from sorted external mem chunks
        writer = FileWriter(outgaf)
        logger.info("Sorting temp files...")
        with timers("write_gaf"):
            for sn in ordered_chromosomes:
                if sn not in chunk_counts:
                    continue

                num_chunks = chunk_counts[sn]
                generators = []
                for i in range(num_chunks):
                    filepath = os.path.join(temp_dir.name, f"{sn}_{i}.tmp")
                    generators.append(stream_chunk(filepath))

                # Flag to track the first line of the chromosome
                is_first_line = True

                for _, final_string in heapq.merge(*generators):
                    # Capture the start offset before writing the first line of the chromosomal block
                    if index_file is not None and is_first_line:
                        index_dict[sn][0] = writer.tell()
                        is_first_line = False

                    writer.write(final_string)

                    # Capture the end offset AFTER writing the last line of the chromosomal block
                    if index_file is not None:
                        index_dict[sn][1] = writer.tell()

        writer.close()
        if index_file is not None:
            if "unknown" in index_dict:
                index_dict.pop("unknown")
            index_dict = dict(index_dict)
            with open(index_file, "wb") as ind:
                pkl.dump(index_dict, ind)

    finally:
        if reader:
            reader.close()
        temp_dir.cleanup()


def process_alignment(line, nodes, offset):
    path = [n for n in PATH_SPLIT_REGEX.split(line[5]) if n]
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
        if n in (">", "<"):
            orient = n
            continue

        node_tags = nodes[n].tags
        sn_tag = node_tags["SN"][1]
        bo_tag = int(node_tags["BO"][1])
        no_tag = int(node_tags["NO"][1])
        sr_tag = int(node_tags["SR"][1])

        # Finding the chromosome where the alignment is
        if sn is None and sr_tag == 0:
            sn = sn_tag
        elif sr_tag == 0:
            assert sn == sn_tag

        if bo_tag == -1 or no_tag == -1:
            logger.debug("[ERR]\tOF:i:%d\tND:Z:%s", offset, n)
            continue
        # Skipping the non-scaffold nodes
        if no_tag != 0:
            continue
        # Only keeping the orientation of the scaffold nodes in the path matching
        orient_list.append(orient)
    # If there are scaffold nodes in both direction, this indicates an inversion
    fwd_count = orient_list.count(">")
    rev_count = orient_list.count("<")
    if fwd_count > 0 and rev_count > 0:
        inv = 1
    # TODO: What to do when the start node that we have is an untagged node?
    # One argument is to just sort them to the end and completely ignore them since they dont belong to any vcf bubble.
    # That works for my tools but not ideal for general purpose.
    if fwd_count < rev_count:
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


'''
# Deprecated.
# Now the sort key is generated directly and not using this function
def gaf_sort_key(al):
    """
    Key function for sorting GAF alignments.
    Returns a tuple evaluated natively at C-speed by Python's Timsort.
    """
    if al[0] == -1:
        return (1, 0, 0, 0, al[3])

    return (0, al[0], al[1], al[2], al[3])
'''

"""
# Deprecated
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
    if al1.NO > al2.NO:
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
"""


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg("gaf", metavar='GAF',
        help="GAF File (can be bgzip-compressed)")
    arg("gfa", metavar='GFA',
        help="GFA file with the sort keys (BO and NO tagged). This is done with gaftools order_gfa")
    arg("--outgaf", default=None,
        help="Output GAF (bgzipped if the file ends with .gz). If omitted, use standard output.")
    arg("--outind", default=None,
        help="Output GAF Sorting Index file. "
        "When --outgaf is not given, no index is created. "
        "If it is given and --outind is not specified, it will have same file name with .gsi extension)")


# fmt: on
def validate(args, parser):
    if args.outind and not args.outgaf:
        parser.error("index path specified but no output gaf path. Please provide an output path.")


def main(args):
    run_sort(**vars(args))
