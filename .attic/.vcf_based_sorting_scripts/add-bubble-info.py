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
    parser.add_argument(
        "--gfa",
        required=True,
        help="GFA file which was used to create the VCF in the GFA-to-VCF pipeline",
    )
    parser.add_argument(
        "--vcf", required=True, help="Input VCF file (output of the GFA-to-VCF conversion pipeline)"
    )
    parser.add_argument("--output", default=None, help="Output GFA File path (Default: sys.stdout)")
    parser.add_argument(
        "--bgzip",
        action="store_true",
        help="Flag to bgzip the output. Can only be given with --output.",
    )
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    options = parser.parse_args()
    setup_logging(options.debug)

    # Checks for argument parsing
    if options.bgzip and not options.output:
        raise RuntimeError(
            "--bgzip flag has been specified but not output path has been defined. Please define the output path."
        )

    add_bubble_tags(options)

    memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    logger.info("\n##### Memory Information #####")
    logger.info("Maximum memory usage: %.3f GB", memory_kb / 1e6)


def add_bubble_tags(options):
    """
    Adding the bubble tags to the nodes in the GFA file.
    Sort Keys will be added to the GFA S-lines. For this bubble tag, only the top level bubbles (LV=0) are considered.

    Two tags will be added:
    1. BO:i:<integer> - This tag gives the bubble order. Odd tag means that the node is a scaffold node (at the start or end of the bubble). Even tag means that the node is inside a bubble.
    2. NO:i:<integer> - This tag gives the node order. For nodes with the same BO tag, they are distinguished based on this NO tag to order each node uniquely.

    Properties of BO tag:
    1. BO tag of the form 2n indicates that node is present in the nth bubble defined in the VCF. Consequently, BO tags 2n-1 and 2n+1 are the scaffold nodes for the nth bubble of the VCF.
    2. This tag has a proper ordering based on the genome coordinate system (Node with BO tag n+1 can be confidently put after Node with BO tag n). There might be some cases where the bubble definition is violated and this order cannot be established. Nodes can be part of two bubbles. Such inconsistent nodes are given a 0 tag.

    Properties of NO tag:
    1. This does not give ordering based on the genome coordinate system.
    2. The ordering is based on the appearance of the node in the Allele Traversal INFO in the VCF file. Sequentially, the non scaffold nodes will be provided a tag starting from 1.
    3. Though this is not ordered based on genome coordinate, this will be a deterministic ordering of some sort.
    4. NO tag will be 0 for scaffold nodes and inconsistent nodes (nodes present in two top level variants)

    General Notes:
    1. So inconsistent nodes will have 0 for both BO and NO tag.
    """

    if options.output is None:
        writer = sys.stdout
    else:
        if options.bgzip:
            writer = libcbgzf.BGZFile(options.output, "wb")
        else:
            writer = open(options.output, "w")
    node = defaultdict(lambda: [-1, -1])
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
        reader = gzip.open(gfa, "rt")
    else:
        reader = open(gfa, "r")

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
        if node[n] == [-1, -1]:
            S_lines_not_in_vcf += 1
            write_to_file(line, writer)
            continue
        bo_tag = "BO:i:"
        no_tag = "NO:i:"
        assert node[n][0] != -1 and node[n][1] != -1, "Node %s has -1,-1 tag" % (n)
        bo_tag += str(node[n][0])
        no_tag += str(node[n][1])
        new_line = line.rstrip() + "\t" + bo_tag + "\t" + no_tag + "\n"
        write_to_file(new_line, writer)

    for line_type, value in line_count.items():
        logger.info("Number of %s lines: %d" % (line_type, value))
    logger.info("Number of S lines not present in VCF: %d" % (S_lines_not_in_vcf))


def read_vcf(vcf, node):
    logger.info("\n##### Parsing VCF file and reading bubble information #####")
    if is_file_gzipped(vcf):
        reader = gzip.open(vcf, "rt")
    else:
        reader = open(vcf, "r")

    total_variants = 0
    inconsistent_bubble_starts = 0
    inconsistent_nodes = defaultdict(lambda: 0)
    bo = 1
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] == "#":
            continue
        fields = line.split("\t")
        info = fields[7].split(";")
        lv = None
        lv = int([x[3:] for x in info if x[0:2] == "LV"][0])
        assert lv is not None, (
            "Variant at chromosome %s position %s does not have a LV field. Need LV info to process"
            % (fields[0], fields[1])
        )

        # Not considering any bubble at LV>0
        if lv != 0:
            continue
        total_variants += 1

        at = None
        at = [x[3:].split(",") for x in info if x[0:2] == "AT"][0]
        assert at is not None, (
            "Variant at chromosome %s position %s does not have a AT field. Need AT info to process"
            % (fields[0], fields[1])
        )

        flanks = list(filter(None, re.split("(>)|(<)", fields[2])))
        assert flanks[0] == flanks[2], "Bad orientation for the bubble."
        if flanks[0] == "<":
            flanks = [">", flanks[-1], ">", flanks[1]]
        assert len(flanks) == 4

        if bo == 1:
            node[flanks[1]][0] = (
                bo  # Since there is no previous bubble, this value has to be assigned
            )
            node[flanks[1]][1] = 0
        else:
            if (
                node[flanks[1]][0] != bo
            ):  # By definition, the end node of bubble 1 should be start node of bubble 2. Checking that and counting such inconsistent starts.
                inconsistent_bubble_starts += 1
                node[flanks[1]][0] = bo
                node[flanks[1]][1] = 0

        no = 1
        for traversal in at:
            tn = list(filter(None, re.split("(>)|(<)", traversal)))
            assert (tn[0] == flanks[0] and tn[1] == flanks[1]) and (
                tn[-2] == flanks[2] and tn[-1] == flanks[3]
            )
            tn = tn[2:-2]
            for n in tn:
                if n == ">" or n == "<":
                    continue
                # This condition checks for the case where the scaffold node is also present inside the bubble.
                if node[n][0] == bo:
                    node[n][1] = (
                        1  # Assigning an NO tag of 1. Usually scaffold nodes are tagged with 0. Here it is tagged with 1 to show that this is also inside a bubble. Hence its orientation cannot be trusted.
                    )
                    continue
                if not (node[n][0] == -1 or node[n][0] == bo + 1):
                    logger.debug(
                        "[0] Node %s present in Bubble %s is also present in Bubble %s at level 0. Anomaly found in variant at chromosome %s:%s."
                        % (n, str((bo + 1) / 2), str(node[n][0] / 2), fields[0], fields[1])
                    )  # Asserting that the node n has not been called before. This will violate the definition of top level bubble.
                    inconsistent_nodes[0] += 1
                    node[n][0] = 0
                    node[n][1] = 0
                    continue
                node[n][0] = bo + 1
                if node[n][1] != -1:
                    continue
                node[n][1] = no
                no += 1
        node[flanks[3]][0] = bo + 2  # Assigning the sort key to the end node of the buble
        node[flanks[3]][1] = 0
        bo += 2
    logger.info("Total Variant Bubbles Processed: %d" % (total_variants))
    logger.info("Inconsistent Bubble Starts: %d" % (inconsistent_bubble_starts))
    for key, value in inconsistent_nodes.items():
        logger.info("Number of nodes found in 2 separate bubbles (at level %d): %d" % (key, value))
    reader.close()


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b"\x1f\x8b"


if __name__ == "__main__":
    main()
