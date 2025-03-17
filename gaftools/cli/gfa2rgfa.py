"""
Converting a GFA file to rGFA format using the W-lines
"""

import os
import sys
import logging
import re
from pysam import libcbgzf
from gaftools.timer import StageTimer
from gaftools.cli import log_memory_usage

logger = logging.getLogger(__name__)
timers = StageTimer()


class Node:
    def __init__(self, line):
        self.ln: int = len(line[2])


class OutputNode:
    def __init__(self, ln, so, sn, sr):
        self.ln: int = ln
        self.so: int = so
        self.sn: str = sn
        self.sr: int = sr

    # converting tag dictionary to string format.
    def tags_to_string(self):
        return f"\tLN:i:{self.ln}\tSO:i:{self.so}\tSN:Z:{self.sn}\tSR:i:{self.sr}"


def run(gfa=None, reference_name="CHM13", reference_tagged=False, seqfile=None, output=None):
    # if reference node info is missing, then need to tag the reference nodes based on the W-line given by user.
    with timers("read_gfa"):
        logger.info("Reading GFA file.")
        node_dict, walks = process_gfa(gfa)
    if not reference_tagged:
        with timers("tag_reference"):
            create_ref_tags(node_dict, walks, reference_name, gfa)
    # if reference node info is already tagged (as is the case with the minigraph-cactus gfa), then just need to fill the rest of the information.
    if seqfile:
        sample_order = get_sample_order(seqfile)
        with timers("tag_assemblies"):
            for index, sample in enumerate(sample_order):
                if sample == reference_name:
                    continue
                create_assembly_tags(node_dict, walks, sample, index, gfa)
    with timers("write_rGFA"):
        stats_counter = write_rGFA(gfa, node_dict, output)
    for key, value in stats_counter.items():
        logger.info(f"Number of {key}: {value}")
    if stats_counter["rank-0 nodes"] == 0:
        if reference_tagged:
            logger.warning(
                "No rank-0 nodes found in the GFA file. Either --reference-tagged was given but no reference nodes were tagged or no. Output rGFA will not have any rank-0 nodes."
            )
        else:
            logger.warning(
                "No rank-0 nodes found in the GFA file. Please check the W-lines for the reference genome. Output rGFA will not have any rank-0 nodes."
            )

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Time to read GFA file:                       %9.2f s", timers.elapsed("read_gfa"))
    logger.info(
        "Time to tag reference nodes:                 %9.2f s", timers.elapsed("tag_reference")
    )
    logger.info(
        "Time to tag assembly nodes:                  %9.2f s", timers.elapsed("tag_assemblies")
    )
    logger.info(
        "Time spend counting arrows:                  %9.2f s", timers.elapsed("counting_arrows")
    )
    logger.info(
        "Time to write rGFA file:                     %9.2f s", timers.elapsed("write_rGFA")
    )
    logger.info("Total time:                                  %9.2f s", total_time)


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gfa', metavar='GFA',
        help='Input GFA file (can be bgzip-compressed). This GFA should have a W-line corresponding to the reference genome or the reference nodes have to be tagged already.')
    arg("--reference-name", metavar='REFERENCE NAME', default='CHM13',
        help="The name of the reference genome given in the W-line. Default: CHM13")
    arg("--reference-tagged", default=False, action="store_true",
        help="Flag to denote reference nodes are already tagged in the GFA.")
    arg("--seqfile", metavar='SEQFILE',
        help='File containing the sequence in which assemblies were given. The first line should be the reference genome. There should be W lines for each assembly in the GFA.')
    arg("--output", metavar='GFA',
        help='Output rGFA.')
# fmt: on


def validate(args, parser):
    if not os.path.exists(args.gfa):
        parser.error(f"GFA file {args.gfa} does not exist.")
    if args.seqfile and not os.path.exists(args.seqfile):
        parser.error(f"Seqfile {args.seqfile} does not exist.")
    if not args.reference_name and not args.reference_tagged:
        parser.error(
            "Either provide the reference name with --reference-name or flag that reference nodes are already tagged with --reference-tagged."
        )


# getting the order of the samples from the seqfile.
def get_sample_order(seqfile):
    order = []
    with open(seqfile, "r") as reader:
        for line in reader:
            sample, _ = line.strip().split("\t")
            order.append(sample)
    return order


# This function processes the GFA file and returns the node dictionary and the walk dictionary.
def process_gfa(gfa):
    nodes = {}
    walks = {}
    gzipped = False
    opened_file = None
    if gfa.endswith(".gz"):
        opened_file = libcbgzf.BGZFile(gfa, "rb")
        gzipped = True
    else:
        opened_file = open(gfa, "r")
    count = 0
    while True:
        offset = opened_file.tell()
        line = opened_file.readline().decode("utf-8") if gzipped else opened_file.readline()
        if not line:
            break
        if line.startswith("S"):
            line = line.strip().split("\t")
            node_id = line[1]
            if len(line) > 3:
                ln = len(line[2])
                for field in line[3:]:
                    tag, _, value = field.split(":")
                    if tag == "LN":
                        ln = int(value)
                    elif tag == "SN":
                        sn = value
                    elif tag == "SO":
                        so = int(value)
                    elif tag == "SR":
                        sr = int(value)
                nodes[int(node_id)] = OutputNode(ln=ln, so=so, sn=sn, sr=sr)
            else:
                nodes[node_id] = Node(line)
        elif line.startswith("W"):
            _, name, hap, _, _, _, _ = line.strip().split(
                "\t",
            )
            try:
                walks[(name, hap)].append(offset)
            except KeyError:
                walks[(name, hap)] = [offset]
        count += 1
        if count % 100000 == 0:
            logger.info(f"Processed {count} lines.")
    opened_file.close()
    return nodes, walks


def yield_walks(file, offsets, gzipped):
    for offset in offsets:
        file.seek(offset)
        line = file.readline().decode("utf-8") if gzipped else file.readline()
        _, _, _, assm_name, start, _, path = line.strip().split(
            "\t",
        )
        yield assm_name, path, int(start)


# creating the rGFA relevant tags for the nodes coming from the reference genome.
# the walk corresponding to reference genome should be (by design) in the 5' to 3' orientation and without cycles.
def create_ref_tags(nodes, walks, reference_name, file):
    gzipped = False
    if file.endswith(".gz"):
        opened_file = libcbgzf.BGZFile(file, "rb")
        gzipped = True
    else:
        opened_file = open(file, "r")
    try:
        walks_list = yield_walks(opened_file, walks[(reference_name, "0")], gzipped)
    except KeyError:
        logger.error(f"No walks found for reference genome {reference_name}.")
        sys.exit(1)
    for assm_name, walk_path, walk_start in walks_list:
        sn = f"{reference_name}#{assm_name}"
        # check how slow this part is.
        walk = list(filter(None, re.split("(>)|(<)", walk_path)))
        so = walk_start  # initializing SO tag
        for i in range(len(walk)):
            if walk[i] in [">", "<"]:
                continue
            id = walk[i]
            assert type(nodes[id]) == Node, f"Reference node {id} is not a rank-0 node."
            nodes[id] = OutputNode(ln=nodes[id].ln, so=so, sn=sn, sr=0)
            so = so + nodes[id].ln
    opened_file.close()


# creating the rGFA relevant tags for the nodes coming from the assemblies.
# for assemblies, the assembly sequence can be in the 3' to 5' orientation. Then the node orientation is affected
# but the SO tags are exclusively based on how the reference sequence looks.
# in hprc-v1.0-minigraph-chm13.gfa.gz, the HG00438#1#JAHBCB010000006.1 contig in the file HG00438.paternal.f1_assembly_v2_genbank.fa is in the reverse orientation.
# you can see how the nodes of that contig (s483177, s483178, s483179) are oriented in the bubble >s13080>s13086
def create_assembly_tags(nodes, walks, sample, index, file):
    if "." in sample:
        sample_name, hap = sample.split(".")
    else:
        sample_name, hap = sample, "0"
    gzipped = False
    if file.endswith(".gz"):
        opened_file = libcbgzf.BGZFile(file, "rb")
        gzipped = True
    else:
        opened_file = open(file, "r")
    try:
        walks_list = yield_walks(opened_file, walks[(sample_name, hap)], gzipped)
    except KeyError:
        logger.error(f"No walks found for sample {sample_name}#{hap}.")
    for assm_name, walk_path, walk_start in walks_list:
        sn = f"{sample_name}#{hap}#{assm_name}"
        walk = list(filter(None, re.split("(>)|(<)", walk_path)))
        so = walk_start  # initializing SO tag
        for i in range(len(walk)):
            if walk[i] in [">", "<"]:
                continue
            id = walk[i]
            if type(nodes[id]) == Node:
                nodes[id] = OutputNode(ln=nodes[id].ln, so=so, sn=sn, sr=index)
            so = so + nodes[id].ln


# writing the rGFA file.
def write_rGFA(gfa, nodes, output):
    stats_counter = {"rank-0 nodes": 0}
    writer_gzipped = False
    writer = None
    if output is None:
        writer = sys.stdout
    elif output.endswith(".gz"):
        writer = libcbgzf.BGZFile(output, "wb")
        writer_gzipped = True
    else:
        writer = open(output, "w")
    reader_gzipped = False
    reader = None
    if gfa.endswith(".gz"):
        reader = libcbgzf.BGZFile(gfa, "rb")
        reader_gzipped = True
    else:
        reader = open(gfa, "r")
    while True:
        line = reader.readline().decode("utf-8") if reader_gzipped else reader.readline()
        if not line:
            break
        if line.startswith("S"):
            # only need to change S lines.
            id = line.split("\t")[1]
            node = nodes[id]
            if type(node) == Node:
                # no tags have been assigned to this node
                node = OutputNode(ln=node.ln, sn="unknown", so=0, sr=-1)
            stats_counter["rank-0 nodes"] += 1 if node.sr == 0 else 0
            new_line = line.strip()
            new_line += node.tags_to_string()
            new_line += "\n"
            writer.write(str.encode(new_line) if writer_gzipped else new_line)
        else:
            writer.write(str.encode(line) if writer_gzipped else line)

    return stats_counter


def main(args):
    run(**vars(args))
