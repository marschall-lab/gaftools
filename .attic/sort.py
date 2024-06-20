"""
Sorting GAF and GFA files
"""

import logging
import sys

from gaftools.cli import log_memory_usage
from gaftools.timer import StageTimer

logger = logging.getLogger(__name__)


def run(gaf_file, gfa_file, output=sys.stdout):
    timers = StageTimer()
    if gaf_file:
        gaf_sort(gaf_file, output)
    else:
        gfa_sort(gfa_file, output)

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


def gfa_sort(gfa_path, out_path=None, return_list=True):
    """This function sorts the given gfa file based on the contig name and start position within the
    contig. Note that it only sorts S lines and leaves the others.
    This can be called from the command line or from another funtion by providing "True" to the
    return_list argument.
    """

    import glob
    from heapq import merge
    import functools
    import gzip
    import os

    gfa_lines = []
    path = "part*.gfa"
    chunk_size = 250000
    chunk_id = 1

    if is_file_gzipped(gfa_path):
        open_gfa = gzip.open
    else:
        open_gfa = open

    with open_gfa(gfa_path, "rt") as gfa_file:
        f_out = open("part_{}.gfa".format(chunk_id), "w")

        for line_num, mapping in enumerate(gfa_file, 1):
            val = mapping.rstrip().split("\t")
            gfa_lines.append(val)

            if not line_num % chunk_size:
                gfa_lines.sort(key=functools.cmp_to_key(compare_gfa))

                for line_count, line in enumerate(gfa_lines):
                    f_out.write("\t".join(line) + "\n")

                logger.info("INFO: Splitting %d" % chunk_id)
                f_out.close()
                gfa_lines = []
                chunk_id += 1
                f_out = open("part_{}.gfa".format(chunk_id), "w")

        if gfa_lines:
            logger.info("INFO: Splitting %d" % chunk_id)
            gfa_lines.sort(key=functools.cmp_to_key(compare_gfa))
            for line_count, line in enumerate(gfa_lines):
                f_out.write("\t".join(line) + "\n")
            f_out.close()
            gfa_lines = []

    chunks = []
    for filename in glob.glob(path):
        chunks += [open(filename, "r")]

    if return_list:
        gfa_s = []
        tmp = merge(*chunks, key=functools.cmp_to_key(compare_gfa2))
        for i in tmp:
            # print(i)
            if i[0] == "S":
                # print(i.rstrip())
                gfa_s.append(i.rstrip().split("\t"))

        for part_file in glob.glob(path):
            if os.path.isfile(part_file):
                os.remove(part_file)

        return gfa_s

    with open(out_path, "w") as f_out:
        f_out.writelines(merge(*chunks, key=functools.cmp_to_key(compare_gfa2)))

    for part_file in glob.glob(path):
        if os.path.isfile(part_file):
            os.remove(part_file)


def gaf_sort(gaf_path, out_path=None):
    """This function sorts a gaf file (mappings to a pangenome graph) in sstable coordinate system based on 1)Contig name 2)Start
    position of the contig's mapping loci.
    """

    import glob
    import functools
    import gzip
    from heapq import merge
    import os

    gaf_lines = []
    path = "part*.gaf"
    chunk_size = 500000
    chunk_id = 1

    if is_file_gzipped(gaf_path):
        open_gaf = gzip.open
    else:
        open_gaf = open

    with open_gaf(gaf_path, "rt") as gaf_file:
        f_out = open("part_{}.gaf".format(chunk_id), "w")

        for line_num, mapping in enumerate(gaf_file, 1):
            val = mapping.rstrip().split("\t")
            gaf_lines.append(val)

            if not line_num % chunk_size:
                gaf_lines.sort(key=functools.cmp_to_key(compare_gaf))

                for line_count, line in enumerate(gaf_lines):
                    f_out.write("\t".join(line) + "\n")

                logger.info("INFO: Splitting", chunk_id)
                f_out.close()
                gaf_lines = []
                chunk_id += 1
                f_out = open("part_{}.gaf".format(chunk_id), "w")

        if gaf_lines:
            logger.info("INFO: Splitting", chunk_id)
            gaf_lines.sort(key=functools.cmp_to_key(compare_gaf))
            for line_count, line in enumerate(gaf_lines):
                f_out.write("\t".join(line) + "\n")
            f_out.close()
            gaf_lines = []

    chunks = []
    for filename in glob.glob(path):
        chunks += [open(filename, "r")]

    with open(out_path, "w") as f_out:
        f_out.writelines(merge(*chunks, key=functools.cmp_to_key(compare_gaf2)))

    for part_file in glob.glob(path):
        if os.path.isfile(part_file):
            os.remove(part_file)


def compare_gfa(ln1, ln2):
    if not ln1[0] == "S" and not ln2[0] == "S":
        # If both are not S lines, then leave it
        return -1
    elif ln1[0] == "S" and not ln2[0] == "S":
        return -1
    elif not ln1[0] == "S" and ln2[0] == "S":
        return 1

    chr1 = [k for k in ln1 if k.startswith("SN:Z:")][0][5:]
    chr2 = [k for k in ln2 if k.startswith("SN:Z:")][0][5:]
    start1 = int([k for k in ln1 if k.startswith("SO:i:")][0][5:])
    start2 = int([k for k in ln2 if k.startswith("SO:i:")][0][5:])

    if chr1 == chr2:
        if start1 < start2:
            return -1
        else:
            return 1
    if chr1 < chr2:
        return -1
    else:
        return 1


def compare_gfa2(ln1, ln2):
    ln1 = ln1.rstrip().split("\t")
    ln2 = ln2.rstrip().split("\t")
    # print(ln1, ln2)
    if not ln1[0] == "S" and not ln2[0] == "S":
        # If both are not S lines, then leave it
        return -1
    elif ln1[0] == "S" and not ln2[0] == "S":
        return -1
    elif not ln1[0] == "S" and ln2[0] == "S":
        return 1

    chr1 = [k for k in ln1 if k.startswith("SN:Z:")][0][5:]
    chr2 = [k for k in ln2 if k.startswith("SN:Z:")][0][5:]
    start1 = int([k for k in ln1 if k.startswith("SO:i:")][0][5:])
    start2 = int([k for k in ln2 if k.startswith("SO:i:")][0][5:])

    if chr1 == chr2:
        if start1 < start2:
            return -1
        else:
            return 1
    if chr1 < chr2:
        return -1
    else:
        return 1


def compare_gaf2(ln1, ln2):
    ln1 = ln1.rstrip().split("\t")
    ln2 = ln2.rstrip().split("\t")

    if ln1[5][0] == ">" or ln1[5][0] == "<":
        chr1 = ln1[5][1:].lower()
    else:
        chr1 = ln1[5].lower()

    if ln2[5][0] == ">" or ln2[5][0] == "<":
        chr2 = ln2[5][1:].lower()
    else:
        chr2 = ln2[5].lower()

    start1 = ln1[7]
    start2 = ln2[7]

    if chr1 == chr2:
        if start1 < start2:
            return -1
        else:
            return 1
    if chr1 < chr2:
        return -1
    else:
        return 1


def compare_gaf(ln1, ln2):
    if ln1[5][0] == ">" or ln1[5][0] == "<":
        chr1 = ln1[5][1:].lower()
    else:
        chr1 = ln1[5].lower()

    if ln2[5][0] == ">" or ln2[5][0] == "<":
        chr2 = ln2[5][1:].lower()
    else:
        chr2 = ln2[5].lower()

    start1 = ln1[7]
    start2 = ln2[7]

    if chr1 == chr2:
        if start1 < start2:
            return -1
        else:
            return 1
    if chr1 < chr2:
        return -1
    else:
        return 1


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b"\x1f\x8b"


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('--gaf', dest = 'gaf_file', metavar='GAF', default = None, required = False, help='GAF File whose coordinates have to be changed')
    arg('--gfa', dest = 'gfa_file', metavar='GFA', default = None, required = False, help='Input GFA file to convert the coordinates')
    arg('-o', '--output', default=sys.stdout, help='Output GAF file. If omitted, use standard output.')
# fmt: on


def validate(args, parser):
    if args.gaf_file and args.gfa_file:
        parser.error("Please input either a GAF or GFA file.")
    if not args.gaf_file and not args.gfa_file:
        parser.error("Please input a GAF or a GFA file to be sorted")


def main(args):
    run(**vars(args))
