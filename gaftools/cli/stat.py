"""
Statistics of a GAF File
"""

import sys
import logging
import gzip
import itertools
from gaftools.gaf import parse_gaf
from gaftools.cli import log_memory_usage
from gaftools.timer import StageTimer

logger = logging.getLogger(__name__)


def run_stat(
    gaf_path,
    cigar_stat=False,
    output=None,
):
    '''This function outputs some statistics of a GAF file. If you run with "--cigar" option, then
    cigar related statistics are also output. This increases the running time more than 10 folds
    because I need to iterate through the cigar of each alignment making it run in O(n^2)
    '''

    timers = StageTimer()

    if output is None:
        output = sys.stdout
    else:
        output = open(output, 'w')

    read_names = set()
    total_aligned_bases = 0
    total_mapq = 0
    total_primary = 0
    total_secondary = 0

    if cigar_stat:
        total_del = 0
        total_ins = 0
        total_x = 0
        total_del_large = 0
        total_ins_large = 0
        total_x_large = 0
        total_match = 0
        total_match_large = 0
        total_perfect = 0
        
    for alignment_count, mapping in enumerate(parse_gaf(gaf_path), 1):
        hashed_readname = hash(mapping.query_name)
        read_names.add(hashed_readname)
        total_aligned_bases += mapping.residue_matches
        total_mapq += mapping.mapping_quality
        if mapping.is_primary:
            total_primary += 1
        else:
            total_secondary += 1

        if cigar_stat:
            '''Cigar string analysis'''
            cigar = mapping.cigar
            all_cigars = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
            if len(all_cigars) == 2:
                total_perfect += 1
            for cnt in range(0, len(all_cigars)-1, 2):
                if all_cigars[cnt+1] == "D":
                    total_del += 1
                    if int(all_cigars[cnt]) >= 50:
                        total_del_large += 1
                elif all_cigars[cnt+1] == "I":
                    total_ins += 1
                    if int(all_cigars[cnt]) >= 50:
                        total_ins_large += 1
                elif all_cigars[cnt+1] == "X":
                    total_x += 1
                    if int(all_cigars[cnt]) >= 50:
                        total_x_large += 1
                elif all_cigars[cnt+1] == "=":
                    total_match += 1
                    if int(all_cigars[cnt]) >= 50:
                        total_match_large += 1
    

    print("Total alignments:", alignment_count, file=output)
    print("\tPrimary:", total_primary, file=output)
    print("\tSecondary:", total_secondary, file=output)
    print("Reads with at least one alignment:", len(read_names), file=output)
    print("Total aligned bases:", str(total_aligned_bases), file=output)
    print("Average mapping quality:", round((total_mapq/alignment_count), 1), file=output)
    
    if cigar_stat:
        print("Cigar string statistics:\n\tTotal deletion regions: %d (%d >50bps)\n\tTotal insertion regions: %d (%d >50bps)\n\tTotal substitution regions: %d (%d >50bps)\n\tTotal match regions: %d (%d >50bps)" % 
          (total_del, total_del_large, total_ins, total_ins_large, total_x, total_x_large,
           total_match, total_match_large), file=output)
        
        print("Total perfect alignments (exact match):", total_perfect, file=output)

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_path', metavar='GAF', help='Input GAF file')
    arg('-o', '--output', default=None, help='Output file. If omitted, use standard output.')
    arg('--cigar', dest='cigar_stat', default=False, action='store_true', help='Outputs cigar related statistics (takes a lot long)')


def validate(args, parser):
    return True


def main(args):
    run_stat(**vars(args))


