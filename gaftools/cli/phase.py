"""
Add phasing information to the GAF file from a haplotag TSV file.

The script uses the TSV file containing the haplotag information generated from WhatsHap's haplotag command.
The H1 and H2 labels for each read are then added to the reads in the GAF file.
"""

import logging
import re
from dataclasses import dataclass

from gaftools.utils import FileWriter
from gaftools.cli import log_memory_usage
from gaftools.timer import StageTimer
from gaftools.gaf import GAF


class IncorrectHaplotypeError(Exception):
    pass


class IncorrectPhaseSetError(Exception):
    pass


class DuplicateHaplotagError(Exception):
    pass


class IncompatibleHaplotagError(Exception):
    pass


HAPLOTYPE_PATTERN = re.compile(r"^H\d+$")


@dataclass
class HaplotagInformation:
    """Class to store the information from haplotag TSV file"""

    chr_name: str
    phase_set: int
    haplotype: str

    def __init__(self, chr_name, haplotype, phase_set):
        self.chr_name = chr_name  # no restriction on this since contig names can be weird
        try:
            if phase_set == "none":
                self.phase_set = -1
            else:
                self.phase_set = int(phase_set)
        except ValueError:
            raise IncorrectPhaseSetError(
                f"Found phaseset {phase_set} in haplotag tsv. Phaseset can only be 'none' or a number."
            )
        if HAPLOTYPE_PATTERN.match(haplotype) or haplotype == "none":
            self.haplotype = haplotype
        else:
            raise IncorrectHaplotypeError(
                f"Found haplotype {haplotype} in haplotag tsv. Haplotag can only be 'none' or of form `H<number>`"
            )

        if self.phase_set == -1 and self.haplotype != "none":
            raise IncompatibleHaplotagError(
                f"Phaseset is 'none' but haplotype is {self.haplotype}. Haplotype should also be 'none'."
            )
        if self.phase_set != -1 and self.haplotype == "none":
            raise IncompatibleHaplotagError(
                f"Haplotype is 'none' but phaseset is {self.phase_set}. Phaseset should also be 'none'."
            )

    def ps_tag(self):
        if self.phase_set == -1:
            return f"ps:Z:{self.chr_name}-none"
        return f"ps:Z:{self.chr_name}-{self.phase_set}"

    def ht_tag(self):
        return f"ht:Z:{self.haplotype}"


logger = logging.getLogger(__name__)
timers = StageTimer()


def run(gaf_file, tsv_file, output=None):
    add_phase_info(gaf_file, tsv_file, output)
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Time to read TSV                             %9.2f s", timers.elapsed("read_tsv"))
    logger.info(
        "Time to read and write GAF                   %9.2f s", timers.elapsed("read_write_gaf")
    )
    logger.info("Time spent on rest:                          %9.2f s", total_time - timers.sum())
    logger.info("Total time:                                  %9.2f s", total_time)


def add_phase_info(gaf_path, tsv_path, out_path):
    """This function adds phasing information to the gaf file using .tsv file of the WhatsHap
    Haplotag... The information is added as two tags ("ps:Z" and "ht:Z") before the cigar string
    """

    logger.info(
        f"Adding phasing information from haplotag TSV {tsv_path} to GAF {gaf_path}. Output written to {out_path if out_path is not None else 'stdout'}."
    )
    tsv_file = open(tsv_path, "r")
    phase = {}
    with timers("read_tsv"):
        for line in tsv_file:
            if line[0] == "#":
                continue
            line_elements = line.rstrip().split("\t")
            if line_elements[0] not in phase:
                phase[line_elements[0]] = HaplotagInformation(
                    line_elements[3], line_elements[1], line_elements[2]
                )
            else:
                raise DuplicateHaplotagError(
                    f"Multiple entries for {line_elements[0]} found. Entries should be unique."
                )
    gaf_out = FileWriter(out_path)
    line_count = 0
    missing_in_tsv = 0
    phased = 0
    gaf_file = GAF(gaf_path)
    with timers("read_write_gaf"):
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
                read_info = phase[gaf_line.query_name]
                gaf_out.write(f"\t{read_info.ps_tag()}\t{read_info.ht_tag()}")
                phased += 1
            else:
                gaf_out.write("\tps:Z:none\tht:Z:none")
            for k in gaf_line.tags.keys():
                gaf_out.write("\t%s%s" % (k, gaf_line.tags[k]))

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
        help='GAF file (can be bgzip-compressed)')
    arg('tsv_file', metavar='TSV',
        help='WhatsHap haplotag TSV file. Refer to https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag')
    arg('-o', '--output', default=None,
        help='Output GAF file (bgzipped if the file ends with .gz). If omitted, output is directed to stdout.')
# fmt: on


def main(args):
    run(**vars(args))
