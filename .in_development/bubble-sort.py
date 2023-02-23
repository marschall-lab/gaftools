import sys
import argparse
import logging
from collections import defaultdict, namedtuple
import functools
import gzip
import re
from pysam import libcbgzf
import resource

logger = logging.getLogger(__name__)

def setup_logging(debug):
    handler = logging.StreamHandler()
    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if debug else logging.INFO)

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    parser.add_argument("--gfa", required=True, help="GFA file with the sort keys (SK tagged)")
    parser.add_argument("--gaf", required=True, help="Input VCF file (output of the GFA-to-VCF conversion pipeline)")
    parser.add_argument("--output", default=None, help="Output GFA File path (Default: sys.stdout)")
    parser.add_argument("--bgzip", action='store_true', help="Flag to bgzip the output. Can only be given with --output.")
    
    options = parser.parse_args()
    setup_logging(options.debug)
    
    validate_arguments(options)
    
    bubble_sort(options)

    memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    logger.info("\n##### Memory Information #####")
    logger.info("Maximum memory usage: %.3f GB", memory_kb / 1e6)
    

def bubble_sort(options):
    """
    Adding the bubble tags to the nodes in the GFA file.
    Sort Keys will be added to the GFA S-lines. Based on the level of the bubble the variant is a part of, the value and the number of values will differ.
    Odd sort keys will represent the backbone while the bubbles will have even sort keys. This is true for all levels of variant bubbles.
    If the sort keys for a node is 62, 3 then we can say that it is in Bubble 31 on the main backbone and there is a level 1 bubble before it and it is not part of the level 1 bubble.
    """
    
    if options.output == None:
        writer = sys.stdout
    else:
        if options.bgzip:
            writer = libcbgzf.BGZFile(options.output, 'wb')
        else:
            writer = open(options.output, "w")
    nodes = defaultdict(lambda: [-1,-1])
    read_gfa(options.gfa, nodes)
    sort(options.gaf, nodes, writer)
    writer.close()
    

def sort(gaf, nodes, writer):
    
    logger.info("\n##### Parsing GFA file and adding sort keys #####")
    if is_file_gzipped(gaf):
        reader = gzip.open(gaf, 'rt')
    else:
        reader = open(gaf, 'r')
    
    Alignment = namedtuple('Alignment', ['offset', 'BO', 'NO', 'start'])
    gaf_alignments = []

    # First pass: Store all the alignment lines as minimally. Just storing line offset and alignment string.
    while True:
        offset = reader.tell()
        line = reader.readline()
        if not line:
            break
        line = line.split('\t')
        bo, no, start = process_alignment(line, nodes)
        gaf_alignments.append(Alignment(offset=offset, BO=bo, NO=no, start=start))
    
    # Sorting the alignments based on BO and NO tag
    gaf_alignments.sort(key=functools.cmp_to_key(compare))
    
    # Writing the sorted file
    for alignment in gaf_alignments:
        off = alignment.offset
        reader.seek(off)
        line = reader.readline()
        write_to_file(line, writer)

def read_gfa(gfa, node):
    
    logger.info("\n##### Parsing GFA file and reading sort key information #####")
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
        if total_nodes == 100000:
            break
        fields = line.split("\t")
        for f in fields:
            if f.startswith("BO:i:"):
                node[fields[1]][0] = int(f[5:])
                tagged_nodes += 1
            if f.startswith("NO:i:"):
                node[fields[1]][1] = int(f[5:])
            
    logger.info("Total Nodes Processed: %d"%(total_nodes))
    logger.info("Nodes with SK tag: %d"%(tagged_nodes))
    reader.close()


def process_alignment(line, nodes):
    path = list(filter(None, re.split('(>)|(<)', line[7])))
    orient = None
    rv = False
    bo = None
    no = None
    start = None
    for n in path:
        if n in [">", "<"]:
            orient = n
            continue
        if nodes[n][0] == -1 or nodes[n][1] == -1:
            raise RuntimeError("Found a node which was not in VCF.")
        if nodes[n][0]%2 == 1 and nodes[n][1] == 0:
            # Here we assume that the orientation of the first scaffold node found is the orientation for all the scaffold nodes. (Here we will have to keep in mind that some scaffold nodes are inside the bubbles also and their orientation can be anything there.)
            # TODO: Have to check for these scaffold nodes that can be present inside bubbles (Cannot trust them)
            # TODO: Have given the NO tag of 1 to these scaffold nodes to distinguish them from normal scaffold nodes which are tagged with 0.
            if orient == ">":
                break
            else:
                rv = True
                break
    if rv:
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

    return bo, no, start


def compare(al1, al2):
    # Since we are sorting in ascending order, al1 is above al2 if it comes before al2.
    # Comparing BO tags
    if al1.BO < al2.BO:
        return 1
    if al1.BO > al2.BO:
        return -1
    
    # Comparing NO tags
    if al1.NO < al2.NO:
        return 1
    if al1.BO > al2.BO:
        return -1
    
    # Comparing start position in node
    if al1.start < al2.start:
        return 1
    if al1.start > al2.start:
        return -1
    
    # This will be execulted only when two reads start at the exact same position at the same node. Then we use offset values. So the read which comes first in the GAF file will come higher up.
    if al1.offset < al2.offset:
        return 1
    if al1.offset > al2.offset:
        return -1


def write_to_file(line, writer):
    try:
        writer.write(line)
    except TypeError:
        writer.write(str.encode(line))


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'


def validate_arguments(options):
    if options.bgzip and not options.output:
        raise RuntimeError("--bgzip flag has been specified but not output path has been defined. Please define the output path.")



if __name__ == "__main__":
    main()