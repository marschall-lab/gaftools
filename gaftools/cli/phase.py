"""
Add phasing information to the bam file
"""

import logging
import sys

from gaftools import __version__
from gaftools.cli import log_memory_usage
from gaftools.cli import CommandLineError
from gaftools.timer import StageTimer

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
    '''This function adds phasing information to the gaf file using .tsv file of the WhatsHap
    Haplotag... The information is added as two tags ("ps:Z" and "ht:Z") before the cigar string
    '''

    import gzip
    import itertools
    
    logger.info("INFO: Adding phasing information...")
    
    tsv_file = open(tsv_path,"r")
 
    phase = {}
    for line in tsv_file:
        line_elements = line.rstrip().split('\t')
        
        if line_elements[0] not in phase:
            tmp = Node(line_elements[3], line_elements[1], line_elements[2])
            phase[line_elements[0]] = tmp
    

    gz_flag = gaf_path[-2:] == "gz"
    if gz_flag:
        gaf_file = gzip.open(gaf_path,"r")
    else:
        gaf_file = open(gaf_path,"r")
    
    gaf_out = open(out_path, "w")
    
    line_count = 0
    missing_in_tsv = 0
    phased = 0
    for gaf_line in gaf_file:
        if gz_flag:
            gaf_line_elements = gaf_line.decode("utf-8").rstrip().split('\t')
        else:
            gaf_line_elements = gaf_line.rstrip().split('\t')
        
        cigar_pos = -1
        for cnt, k in enumerate(gaf_line_elements):
            if k.startswith("cg:Z:"):
                cigar_pos = cnt

        
        if line_count != 0:
            gaf_out.write("\n")
        
        line_count += 1
        for i in gaf_line_elements[:cigar_pos]:
            gaf_out.write("%s\t"%i)
        
        in_tsv = True
        if gaf_line_elements[0] not in phase:
            #print("%s not in .tsv file" %gaf_line_elements[0])
            missing_in_tsv += 1
            in_tsv = False
        
        if in_tsv and phase[gaf_line_elements[0]].haplotype != "none":
            gaf_out.write("ps:Z:%s-%s\tht:Z:%s\t" % (phase[gaf_line_elements[0]].chr_name,
                                               phase[gaf_line_elements[0]].phase_set, phase[gaf_line_elements[0]].haplotype))
            phased +=1
        else:
            gaf_out.write("ps:Z:none\tht:Z:none")

        for i in gaf_line_elements[cigar_pos:]:
            gaf_out.write("\t%s"%i)
    
    logger.info("INFO: Added phasing info (ps:Z and ht:Z) for %d reads out of %d GAF lines - (%d reads are missing in .tsv)" %(phased, line_count, missing_in_tsv))

    gaf_file.close()
    gaf_out.close()


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    
    # Positional arguments
    arg('gaf_file', metavar='GAF', help='GAF File')
    arg('tsv_file', metavar='phase', help='WhatsHap Haplotag file (.tsv)')
    arg('-o', '--output', default=sys.stdout,
        help='Output GAF file. If omitted, use standard output.')

# fmt: on
def validate(args, parser):
    return True

def main(args):
    run(**vars(args))
