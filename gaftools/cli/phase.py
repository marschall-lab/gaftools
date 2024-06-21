"""
Add phasing information to the GAF file from a haplotag TSV file.

The script uses the TSV file containing the haplotag information generated from WhatsHap's haplotag command.
The H1 and H2 labels for each read are then added to the reads in the GAF file.
"""

import logging
import sys

from gaftools.cli import log_memory_usage
from gaftools.timer import StageTimer
from gaftools.gaf import GAF


class Node:
    def __init__(self, chr_name, haplotype, phase_set):
        self.chr_name = chr_name
        self.phase_set = phase_set
        self.haplotype = haplotype


logger = logging.getLogger(__name__)


def run(gaf_file, tsv_file, output=sys.stdout):
    timers = StageTimer()
    add_phase_info(gaf_file, tsv_file, output)
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


def add_phase_info(gaf_path, tsv_path, out_path):
    """This function adds phasing information to the gaf file using .tsv file of the WhatsHap
    Haplotag... The information is added as two tags ("ps:Z" and "ht:Z") before the cigar string
    """

    logger.info("INFO: Adding phasing information...")

    tsv_file = open(tsv_path, "r")

    phase = {}
    for line in tsv_file:
        line_elements = line.rstrip().split("\t")

        if line_elements[0] not in phase:
            tmp = Node(line_elements[3], line_elements[1], line_elements[2])
            phase[line_elements[0]] = tmp

    gaf_out = open(out_path, "w")

    line_count = 0
    missing_in_tsv = 0
    phased = 0
    gaf_file = GAF(gaf_path)
    for gaf_line in gaf_file.read_file():
        if line_count != 0:
            gaf_out.write("\n")

        line_count += 1

        gaf_out.write(
            "%s\t%s\t%s\t%s\t+\t%s\t%d\t%d\t%d\t%d\t%d\t%d"
            % (
                gaf_line.query_name,
                gaf_line.query_length,
                gaf_line.query_start,
                gaf_line.query_end,
                gaf_line.path,
                gaf_line.path_length,
                gaf_line.path_start,
                gaf_line.path_end,
                gaf_line.residue_matches,
                gaf_line.alignment_block_length,
                gaf_line.mapping_quality,
            )
        )

        in_tsv = True
        if gaf_line.query_name not in phase:
            missing_in_tsv += 1
            in_tsv = False

        if in_tsv and phase[gaf_line.query_name].haplotype != "none":
            gaf_out.write(
                "\tps:Z:%s-%s\tht:Z:%s\t"
                % (
                    phase[gaf_line.query_name].chr_name,
                    phase[gaf_line.query_name].phase_set,
                    phase[gaf_line.query_name].haplotype,
                )
            )
            phased += 1
        else:
            gaf_out.write("\tps:Z:none\tht:Z:none")

        for k in gaf_line.tags.keys():
            gaf_out.write("\t%s:%s" % (k, gaf_line.tags[k]))

        gaf_out.write("\t%s" % gaf_line.cigar)

    logger.info(
        "INFO: Added phasing info (ps:Z and ht:Z) for %d reads out of %d GAF lines - (%d reads are missing in .tsv)"
        % (phased, line_count, missing_in_tsv)
    )
    gaf_file.close()
    tsv_file.close()
    gaf_out.close()


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument

    # Positional arguments
    arg('gaf_file', metavar='GAF',
        help='Input GAF file (can be bgzip-compressed)')
    arg('tsv_file', metavar='TSV',
        help='WhatsHap haplotag TSV file. Refer to https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag')
    arg('-o', '--output', default=sys.stdout,
        help='Output GAF file. If omitted, output is directed to standard output.')


# fmt: on
def validate(args, parser):
    return True


def main(args):
    run(**vars(args))
