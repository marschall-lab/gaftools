"""
Convert Coordinate Systems between the stable system and unstable system
"""
import logging
import sys
import platform

from gaftools import __version__
from gaftools.cli import log_memory_usage
from gaftools.cli import CommandLineError


logger = logging.getLogger(__name__)

def run(
    gaf_file,
    gfa_file,
    output=sys.stdout,
    unstable=False,
    stable=False
):
    assert (unstable != stable)
    if (unstable):
        stable_to_unstable(gaf_file, gfa_file)
    else:
        unstable_to_stable(gaf_file, gfa_file)


def stable_to_unstable(gaf_file, gfa_file):
    
    exit

def unstable_to_stable(gaf_file, gfa_file):
    exit



# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_file', metavar='GAF', help='GAF File whose coordinates have to be changed')
    arg('gfa_file', metavar='GFA', help='Input GFA file to conver the coordinates')
    arg('-o', '--output', default=sys.stdout,
        help='Output GAF file.'
        'If omitted, use standard output.')
    arg('--unstable', dest='unstable', default=False, action='store_true',
        help='Convert to Unstable Coordinates')
    arg('--stable', dest='stable', default=False, action='store_true',
        help='Convert to Stable Coordinates')
    
# fmt: on
def validate(args, parser):
    if not args.unstable and not args.stable:
        parser.error("Either one of --unstable or --stable has to be provided")
    if args.unstable and args.stable:
        parser.error("Specify one of --unstable or --stable. Both cannot be provided")

def main(args):
    run(**vars(args))
