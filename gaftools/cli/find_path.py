"""
Find the genomic sequence of a given connected GFA path.
"""

import logging
import sys
import os
from gaftools.cli import log_memory_usage
from gaftools.timer import StageTimer
from gaftools.gfa import GFA
from gaftools.utils import FileWriter

logger = logging.getLogger(__name__)


def run(gfa_path, path=None, paths_file=None, keep_going=True, output=None, fasta=False):
    timers = StageTimer()

    if not path and not paths_file:
        logger.error("Either an input path (--p, --path) or (--paths_file) must be specified")
        sys.exit(1)

    logger.info(f"Reading GFA file {gfa_path}")
    graph = GFA(gfa_path)
    if path:
        if path[0] in [">", "<"]:
            # detected node path
            nodes = [path]
            seq = graph.extract_path(path)
            if not seq:
                if keep_going:
                    logger.warning(
                        f"The path {path} does not exist in the GFA, will keep going to next path if given"
                    )
                else:
                    logger.error(f"The path {path} does not exist in the GFA")
                    sys.exit(1)

            path_seqs = [seq]
        else:
            logger.error(f"The input path {path} is not a valid node path")
            sys.exit(1)

    else:
        if os.path.exists(paths_file):
            # detected file
            reader = open(paths_file, "r")
            nodes = []
            path_seqs = []
            for line in reader:
                nodes.append(line.strip())
                seq = graph.extract_path(nodes[-1])
                if not seq:
                    if keep_going:
                        logger.warning(
                            f"The path {line} does not exist in the GFA, will keep going to next path if given"
                        )
                    else:
                        logger.error(f"The path {line} does not exist in the GFA")
                        sys.exit(1)
                path_seqs.append(seq)
            reader.close()
        else:
            logger.error(f"The input file {paths_file} does not exist")
            sys.exit(1)

    writer = FileWriter(output)
    if fasta:
        for node, path_seq in zip(nodes, path_seqs):
            if path_seq:
                writer.write(f">seq_{node}\n")
                writer.write(f"{path_seq}\n")
    else:
        for node, path_seq in zip(nodes, path_seqs):
            if path_seq:
                writer.write(f"{path_seq}\n")

    writer.close()

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


# fmt:off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg("gfa_path", metavar="GFA",
        help="GFA file (can be bgzip-compressed)")
    arg("-p", "--path", default=None,
        help='GFA node path to retrieve the sequence (e.g., ">s82312<s82313") with the quotes')
    arg("--paths_file", default=None,
        help="File containing the paths to retrieve the sequences for, each path on new line")
    arg("-k", "--keep-going", action="store_true",
         help="Keep going after instead of stopping when a path does not exist")
    arg("-o", "--output", default=None,
        help="Output file. If omitted, use standard output.")
    arg("-f", "--fasta", action="store_true",
        help="Flag to output the sequence as a FASTA file with the seqeunce named seq_<node path>")

# fmt:on


def main(args):
    run(**vars(args))
