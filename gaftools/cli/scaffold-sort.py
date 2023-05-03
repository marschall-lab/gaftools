'''
Scaffold Sort:

The basis of the script is the BO and NO tags defined by Tobias Marschall in his order_gfa.py script. The script defines scaffold nodes as nodes whose removal increases total number of connected components.
Using this definition of scaffold nodes, the sorting has been done.

The main difference between this sorting and the bubble-sort.py script is that here the the entire alignment is processed and used to determine whether the alignment is reveresed.
Here inversions are not taken into account.

Index file created:
The index is dictionary created using pickle library which contains the reference contig names as keys and the offset the alignments begin and end as a list of size 2.

Other things the code does:
It adds some tags into the sorted gaf file. The tags are:
    1. bo:i: - This is the BO tags from the GFA carried forward. The sorting uses a BO tag for each alignment, the BO tag for the start node (or end node based on overall orientation) of the alignment path.
    2. sn:Z: - The name of the reference contig it mapped to.
    3. iv:i: - 1 if the alignment path has an inversion. 0 otherwise.
'''

import sys
import argparse
import logging
from collections import defaultdict, namedtuple
import functools
import gzip
import re
from pysam import libcbgzf
import resource
import pickle as pkl
from gaftools.timer import StageTimer

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
    parser.add_argument("--gfa", required=True, help="GFA file with the sort keys (BO and NO tagged)")
    parser.add_argument("--gaf", required=True, help="Input GAF File")
    parser.add_argument("--outgaf", default=None, help="Output GAF File path (Default: sys.stdout)")
    parser.add_argument("--outind", default=None, help="Output Index File path for the GAF file. (When --outgaf is not given, no index is created. If it is given and --outind is not specified, it will have same file name with .gai extension)")
    parser.add_argument("--bgzip", action='store_true', help="Flag to bgzip the output. Can only be given with --output.")
    
    options = parser.parse_args()
    setup_logging(options.debug)
    
    validate_arguments(options)
    
    bubble_sort(options)

    memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    logger.info("\nMemory Information")
    logger.info("  Maximum memory usage: %.3f GB", memory_kb / 1e6)
    logger.info("\nTime Summary:")
    logger.info("  Time to parse GFA file: %.3f"%(timers.elapsed("read_gfa")))
    logger.info("  Total time to sort GAF file: %.3f"%(timers.elapsed("total_sort")))
    logger.info("    Time to parse GAF file: %.3f"%(timers.elapsed("read_gaf")))
    logger.info("    Time to sort GAF file: %.3f"%(timers.elapsed("sort_gaf")))
    logger.info("    Time to write GAF file: %.3f"%(timers.elapsed("write_gaf")))
    

def bubble_sort(options):
    
    if options.outgaf == None:
        writer = sys.stdout
        index_file = None
    else:
        if options.bgzip:
            writer = libcbgzf.BGZFile(options.outgaf, 'wb')
        else:
            writer = open(options.outgaf, "w")
        if options.outind:
            index_file = options.outind
        else:
            index_file = options.outgaf+".gai"
    nodes = defaultdict(lambda: [-1,-1,-1,-1])
    index_dict = defaultdict(lambda: [None, None])
    with timers("read_gfa"):
        read_gfa(options.gfa, nodes)
    with timers("total_sort"):
        sort(options.gaf, nodes, writer, index_dict)
    writer.close()
    index_dict = dict(index_dict)
    with open(index_file, "wb") as ind:
        pkl.dump(index_dict, ind)
    
    

def sort(gaf, nodes, writer, index_dict):
    
    logger.info("Parsing GAF file and sorting it")
    if is_file_gzipped(gaf):
        reader = gzip.open(gaf, 'rt')
    else:
        reader = open(gaf, 'r')
    
    Alignment = namedtuple('Alignment', ['offset', 'BO', 'NO', 'start', 'inv', 'sn'])
    gaf_alignments = []
    count_inverse = 0
    # First pass: Store all the alignment lines as minimally. Just storing line offset and alignment string.
    with timers("read_gaf"):
        while True:
            offset = reader.tell()
            line = reader.readline()
            if not line:
                break
            line = line.split('\t')
            bo, no, start, inv, sn = process_alignment(line, nodes, offset)
            if inv == 1:
                count_inverse += 1
            gaf_alignments.append(Alignment(offset=offset, BO=bo, NO=no, start=start, inv=inv, sn=sn))
    
    logger.info("\tNumber of alignments with invertions: %d"%(count_inverse))
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
            line = reader.readline().rstrip()
            line += "\tbo:i:%d\tsn:Z:%s\tiv:i:%d\n"%(alignment.BO, alignment.sn, alignment.inv)
            out_off = writer.tell()
            if index_dict[alignment.sn][0] == None:
                index_dict[alignment.sn][0] = out_off
                index_dict[alignment.sn][1] = out_off
            else:
                index_dict[alignment.sn][1] = out_off
            write_to_file(line, writer)
    index_dict.pop('unknown')
    reader.close()

def read_gfa(gfa, node):
    
    logger.info("Parsing GFA file and reading sort key information")
    if is_file_gzipped(gfa):
        reader = gzip.open(gfa, 'rt')
    else:
        reader = open(gfa, 'r')
    
    total_nodes = 0
    tagged_nodes = 0
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] != 'S':
            continue
        total_nodes += 1
        fields = line.split("\t")
        for f in fields:
            if f.startswith("BO:i:"):
                node[fields[1]][0] = int(f[5:])
                tagged_nodes += 1
            if f.startswith("NO:i:"):
                node[fields[1]][1] = int(f[5:])
            if f.startswith("SR:i:"):
                node[fields[1]][2] = int(f[5:])
            if f.startswith("SN:Z:"):
                node[fields[1]][3] = f[5:]
            
    logger.info("\tTotal Nodes Processed: %d"%(total_nodes))
    logger.info("\tNodes with tags: %d"%(tagged_nodes))
    reader.close()


def process_alignment(line, nodes, offset):
    path = list(filter(None, re.split('(>)|(<)', line[5])))
    orient = None
    # If there is no scaffold node present in the alignment, then assuming that it is in the correct orientation.
    # TODO: Need to find a way to deal with such alignments.
    rv = False
    bo = None
    no = None
    start = None
    orient_list = []
    inv = 0     # Whether the alignment has an inversion
    sn = None
    for n in path:
        if n in [">", "<"]:
            orient = n
            continue
        
        # Finding the chromosome where the alignment is
        if sn == None and nodes[n][2] == 0:
            sn = nodes[n][3]
        elif nodes[n][2] == 0:
            assert sn == nodes[n][3]
        
        if nodes[n][0] == -1 or nodes[n][1] == -1:
            logger.debug("[ERR]\tOF:i:%d\tND:Z:%s"%(offset,n))
            continue
        # Skipping the non-scaffold nodes
        if nodes[n][1] != 0:
            continue
        # Only keeping the orientation of the scaffold nodes in the path matching
        orient_list.append(orient)
    # If there are scaffold nodes in both direction, this indicates an inversion
    if orient_list.count('>') !=0 and orient_list.count('<') != 0:
        inv = 1
    # TODO: What to do when the start node that we have is an untagged node?
    # One argument is to just sort them to the end and completely ignore them since they dont belong to any vcf bubble.
    # That works for my tools but not ideal for general purpose.
    if orient_list.count('>') < orient_list.count('<'):
        l = int(line[6])
        e = int(line[8])
        start = l-e
        n = path[-1]
        bo = nodes[n][0]
        no = nodes[n][1]
    else:
        start = int(line[7])
        n = path[1]
        bo = nodes[n][0]
        no = nodes[n][1]
    if sn == None:
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


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'


def validate_arguments(options):
    if options.bgzip and not options.outgaf:
        raise RuntimeError("--bgzip flag has been specified but not output path has been defined. Please define the output path.")
    if options.outind and not options.outgaf:
        raise RuntimeError("index path specified but no output gaf path. Please provide an output path.")


if __name__ == "__main__":
    main()
