"""
Statistics of a GAF File
"""

import sys
import logging
import itertools
from gaftools.gaf import GAF, Read
from gaftools.cli import log_memory_usage
from gaftools.timer import StageTimer
from gaftools.gfa import GFA

logger = logging.getLogger(__name__)


def run_find_path(gfa_path, input_path, output=None):
	timers = StageTimer()
	
	graph = GFA(gfa_path)

	path_seq = graph.extract_path(input_path)
	print(path_seq)	
	
	logger.info("\n== SUMMARY ==")
	total_time = timers.total()
	log_memory_usage()
	logger.info("Total time:                                  %9.2f s", total_time)


def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg("gfa_path", metavar="GFA", help="Input GFA file (can be bgzip-compressed)")
    arg("input_path", metavar="path", help="GFA path to retrieve the sequence (e.g., \">s82312<s82313\").")
    arg("-o", "--output", default=None, help="Output file. If omitted, use standard output.")

def validate(args, parser):
    return True


def main(args):
    run_find_path(**vars(args))
