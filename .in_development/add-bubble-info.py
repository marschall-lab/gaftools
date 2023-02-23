import sys
import argparse
import logging
from collections import defaultdict
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
    parser.add_argument("--gfa", required=True, help="GFA file which was used to create the VCF in the GFA-to-VCF pipeline")
    parser.add_argument("--vcf", required=True, help="Input VCF file (output of the GFA-to-VCF conversion pipeline)")
    parser.add_argument("--output", default=None, help="Output GFA File path (Default: sys.stdout)")
    parser.add_argument("--bgzip", action='store_true', help="Flag to bgzip the output. Can only be given with --output.")
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    options = parser.parse_args()
    setup_logging(options.debug)
    
    # Checks for argument parsing
    if options.bgzip and not options.output:
        raise RuntimeError("--bgzip flag has been specified but not output path has been defined. Please define the output path.")

    add_bubble_tags(options)

    memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    logger.info("\n##### Memory Information #####")
    logger.info("Maximum memory usage: %.3f GB", memory_kb / 1e6)
    

def add_bubble_tags(options):
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
    node = defaultdict(lambda: defaultdict(lambda: -1))
    read_vcf(options.vcf, node)
    modify_gfa(options.gfa, node, writer)
    writer.close()
    
def write_to_file(line, writer):
    try:
        writer.write(line)
    except TypeError:
        writer.write(str.encode(line))
    

def modify_gfa(gfa, node, writer):
    
    logger.info("\n##### Parsing GFA file and adding sort keys #####")
    if is_file_gzipped(gfa):
        reader = gzip.open(gfa, 'rt')
    else:
        reader = open(gfa, 'r')
    
    S_lines_not_in_vcf = 0
    line_count = defaultdict(lambda: 0)
    while True:
        line = reader.readline()
        if not line:
            break
        line_count[line[0]] += 1
        if line[0] != "S":
            write_to_file(line, writer)
            continue
        fields = line.split("\t")
        n = fields[1]
        if len(node[n].keys()) == 0:
            S_lines_not_in_vcf += 1
            write_to_file(line, writer)
            continue
        sk_tag = "SK:B:I"
        for lv in node[n].keys():
            assert node[n][lv] != -1, "Node %s has a defined sort key at level %d but the key is -1. If the level has been defined, it should have proper key."%(n, lv)
            sk_tag += ","+str(node[n][lv])
        new_line = line.rstrip()+"\t"+sk_tag+"\n"
        write_to_file(new_line, writer)
    
    for line_type, value in line_count.items():
        logger.info("Number of %s lines: %d"%(line_type, value))
    logger.info("Number of S lines not present in VCF: %d"%(S_lines_not_in_vcf))
    
# TODO: Work on inconsitent nodes (nodes present in multiple bubbles) at level 1 and above.
def read_vcf(vcf, node):
    
    logger.info("\n##### Parsing VCF file and reading bubble information #####")
    if is_file_gzipped(vcf):
        reader = gzip.open(vcf, 'rt')
    else:
        reader = open(vcf, 'r')
    
    total_variants = 0
    nested_bubble = 0
    inconsistent_bubble_starts = 0
    top_bubble = 0
    inconsistent_nodes = defaultdict(lambda: 0)
    sk = 1
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] == '#':
            continue
        total_variants += 1
        fields = line.split("\t")
        flanks = list(filter(None, re.split('(>)|(<)', fields[2])))
        assert flanks[0] == flanks[2], "Bad orientation for the bubble."
        if flanks[0] == "<":
            flanks = [">", flanks[-1], ">", flanks[1]]
        info = fields[7].split(";")
        lv = None
        lv = int([x[3:] for x in info if x[0:2] == "LV"][0])
        assert lv != None, "Variant at chromosome %s position %s does not have a LV field. Need LV info to process"%(fields[0], fields[1])
        at = [x[3:].split(",") for x in info if x[0:2] == "AT"][0]
        assert len(flanks) == 4
        if lv == 0:
            top_bubble += 1
            if sk == 1:
                node[flanks[1]][0] = sk    #Since there is no bubble, this value has to be assigned
            else:
                if node[flanks[1]][0] != sk:    #By definition, the end node of bubble 1 should be start node of bubble 2. Asserting that
                    inconsistent_bubble_starts += 1
                    node[flanks[1]][0] = sk
            
            for traversal in at:
                tn = list(filter(None, re.split('(>)|(<)', traversal)))
                assert (tn[0] == flanks[0] and tn[1] == flanks[1]) and (tn[-2] == flanks[2] and tn[-1] == flanks[3])
                tn = tn[2:-2]
                for n in tn:
                    if n == ">" or n == "<":
                        continue
                    if node[n][lv] == sk:
                        continue
                    if not (node[n][0] == -1 or node[n][0] == sk+1):
                        logger.debug("[0] Node %s present in Bubble %s is also present in Bubble %s at level 0. Anomaly found in variant at chromosome %s:%s."%(n, str((sk+1)/2), str(node[n][0]/2), fields[0], fields[1]))        #Asserting that the node n has not been called before. This will violate the definition of top level bubble.
                        inconsistent_nodes[0] += 1
                        node[n][0] = 0
                        continue
                    node[n][0] = sk+1
            node[flanks[3]][0] = sk+2  #Assigning the sort key to the end node of the buble
            sk += 2
        else:
            nested_bubble += 1
            if  node[flanks[1]][lv] == -1:
                new_key = 1
                node[flanks[1]][lv] = new_key
            else:
                # If the node already exists and was first initiallizeed inside a bubble, then make new key as 1. More importance to flanks than inside bubbles.
                if (node[flanks[1]][lv]%2 != 1):
                    new_key = 1
                    node[flanks[1]][lv] = new_key
                # If the node already exists and was defined in the flank region, then keep is like that and update new_key.
                else:
                    new_key = node[flanks[1]][lv]
                assert new_key%2 == 1
            for traversal in at:
                tn = list(filter(None, re.split('(>)|(<)', traversal)))[2:-2]
                for n in tn:
                    if n == ">" or n == "<":
                        continue
                    if node[n][lv] == new_key:
                        continue
                    if not(node[n][lv] == -1 or node[n][lv] == new_key+1):
                        inconsistent_nodes[lv] += 1
                        if node[n][lv]%2 == 0:
                            node[n][lv] = 0
                        logger.debug("[%d] Node %s present in Bubble %s is also present in Bubble %s at level %d. Anomaly found in variant at chromosome %s:%s."%(lv, n, str((new_key+1)/2), str(node[n][lv]/2), lv, fields[0], fields[1]))
                        continue
                    node[n][lv] = new_key+1
            node[flanks[3]][lv] = new_key+2
    logger.info("Total Variant Bubbles Processed: %d"%(total_variants))
    logger.info("Top Level Bubbles: %d"%(top_bubble))
    logger.info("Nested Bubbles: %d"%(nested_bubble))
    logger.info("Inconsistent Bubble Starts: %d"%(inconsistent_bubble_starts))
    for key, value in inconsistent_nodes.items():
        logger.info("Number of nodes found in 2 separate bubbles (at level %d): %d"%(key, value))
    reader.close()


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'

if __name__ == "__main__":
    main()