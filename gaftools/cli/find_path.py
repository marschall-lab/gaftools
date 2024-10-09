"""
Find the genomic sequence of a given GFA path.
"""

import logging
from gaftools.cli import log_memory_usage
from gaftools.timer import StageTimer
from gaftools.gfa import GFA
import sys

logger = logging.getLogger(__name__)


def run(gfa_path, input_path, output=None, fasta=False):
    timers = StageTimer()

    graph = GFA(gfa_path)
    if input_path[0] in [">", "<"]:
        # detected node path
        nodes = [input_path]
        path_seqs = [graph.extract_path(input_path)]
    else:
        # detected file
        reader = open(input_path, "r")
        nodes = []
        path_seqs = []
        for line in reader:
            nodes.append(line.strip())
            path_seqs.append(graph.extract_path(nodes[-1]))
        reader.close()

    if output is None:
        writer = sys.stdout
    else:
        writer = open(output, "w")
    if fasta:
        for node, path_seq in zip(nodes, path_seqs):
            print(f">seq_{node}", file=writer)
            print(path_seq, file=writer)
    else:
        for node, path_seq in zip(nodes, path_seqs):
            print(path_seq, file=writer)

    if output is not None:
        writer.close()

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg("gfa_path", metavar="GFA", help="Input GFA file (can be bgzip-compressed)")
    arg(
        "input_path",
        metavar="path",
        help='GFA node path to retrieve the sequence (e.g., ">s82312<s82313") OR a filepath containing node paths in different lines',
    )
    arg("-o", "--output", default=None, help="Output file. If omitted, use standard output.")
    arg(
        "-f",
        "--fasta",
        action="store_true",
        help="Flag to output the sequence as a FASTA file with the seqeunce named seq_<node path>",
    )


def validate(args, parser):
    return True


def main(args):
    run(**vars(args))
