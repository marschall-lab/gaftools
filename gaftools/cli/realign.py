"""
Realign GAF file using wavefront alignment algorithm (WFA)
"""

import logging
import sys
import pysam
import gaftools.gaf
from gaftools import __version__
from gaftools.cli import log_memory_usage
from gaftools.cli import CommandLineError
from gaftools.timer import StageTimer
from gaftools.GFA import GFA
from pywfa.align import (WavefrontAligner, cigartuples_to_str)

logger = logging.getLogger(__name__)


def run_realign(
    gaf,
    graph,
    fasta,
    output=None,
    ext=False
):
    timers = StageTimer()
    
    if output is None:
        output = sys.stdout
    else:
        output = open(output, 'w')
    
    realign_gaf(gaf, graph, fasta, output, ext)
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


def overlap_ratio(x_start,x_end, y_start, y_end):

    overlap = max(0, min(x_end, y_end) - max(x_start, y_start))
    total_length = x_end - x_start + y_end - y_start
    x_length = x_end - x_start
    y_length = y_end - y_start

    return max(2 * (overlap / total_length) , (overlap / x_length), (overlap/y_length))


def filter_duplicates(aln):
    import functools

    for k in aln.keys():
        aln[k].sort(key=functools.cmp_to_key(gaftools.gaf.compare_aln))
     
    for read_name, mappings in aln.items():
        if len(mappings) == 1:
            continue
        for cnt, line in enumerate(mappings):
            if line.duplicate == True:
                continue
            for cnt2, line2 in enumerate(mappings):
                if cnt == cnt2 or line2.duplicate == True:
                    continue
                
                sim = overlap_ratio(line.query_start, line.query_end, line2.query_start, line2.query_end)
                
                if sim > 0.75:
                    if(line.score > line2.score):
                        mappings[cnt2].duplicate = True
                        #print("......Duplicate")
                    else:
                        mappings[cnt].duplicate = True
                        #print("......Duplicate, breaking")
                        break


def write_alignments(aln, output):
    
    #Write the alignment back to the GAF
    for read_name, mappings in aln.items():
        for gaf_line in mappings:
            if not gaf_line.duplicate:


                #Write the alignment back to the GAF
                output.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" %(gaf_line.query_name,
                        gaf_line.query_length, gaf_line.query_start, gaf_line.query_end, gaf_line.strand,
                        gaf_line.path, gaf_line.path_length, gaf_line.path_start, gaf_line.path_end,
                        match, cigar_len, gaf_line.mapping_quality))

                for k in gaf_line.tags.keys():
                    output.write("\t%s:%s"%(k,gaf_line.tags[k]))
 
                output.write("\tcg:Z:%s\n" %gaf_line.cigar)


def wfa_alignment(aln, gaf_line, ref, query, path_start, output, extended):

    aligner = WavefrontAligner(ref)
    if extended:
        res = aligner(query, clip_cigar = True, min_aligned_bases_left = 30, min_aligned_bases_right = 30)
    else:
        res = aligner(query, clip_cigar = False)
    
   
    match, mismatch, cigar_len, ins, deletion, soft_clip = 0, 0, 0, 0, 0, 0
    cigar = ""
    
    for op_type, op_len in res.cigartuples:
        if op_type == 0:
            match += op_len
            cigar += str(op_len) + "="
        elif op_type == 1:
            ins += op_len
            cigar += str(op_len) + "I"
        elif op_type == 2:
            deletion += op_len
            cigar += str(op_len) + "D"
        elif op_type == 4:
            soft_clip += op_len
        elif op_type == 8:
            mismatch += op_len
            cigar += str(op_len) + "X"
        else:
            assert False
        cigar_len += op_len
 

    if extended:
        if match < 30:
            return
        if gaf_line.query_name in aln:
            aln[gaf_line.query_name].append(Alignment("", gaf_line.query_length, res.text_start,
                                               res.text_end, gaf_line.strand, gaf_line.path,
                                               gaf_line.path_length, path_start + res.pattern_start,
                                               path_start + res.pattern_end, match, 0, 0, False, cigar,
                                               cigar_len, res.score))
        else:
            aln[gaf_line.query_name] = [Alignment("", gaf_line.query_length, res.text_start,
                                               res.text_end, gaf_line.strand, gaf_line.path,
                                               gaf_line.path_length, path_start + res.pattern_start,
                                               path_start + res.pattern_end, match, 0, 0, False,
                                               cigar, cigar_len, res.score)]
    else:
        cigar = aligner.cigarstring.replace("M", "=")

        #Write the alignment back to the GAF
        output.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" %(gaf_line.query_name,
                        gaf_line.query_length, gaf_line.query_start, gaf_line.query_end, gaf_line.strand,
                        gaf_line.path, gaf_line.path_length, gaf_line.path_start, gaf_line.path_end,
                        match, cigar_len, gaf_line.mapping_quality))

        for k in gaf_line.tags.keys():
            output.write("\t%s:%s"%(k,gaf_line.tags[k]))
 
        output.write("\tcg:Z:%s\n" %cigar)


def realign_gaf(gaf, graph, fasta, output, extended):
    """
    Uses pyWFA (https://github.com/kcleal/pywfa)
    """

    fastafile = pysam.FastaFile(fasta)
    graph_obj = GFA(graph)

    aln = {}
    for cnt, line in enumerate(gaftools.gaf.parse_gaf(gaf)):
        path_sequence = graph_obj.extract_path(line.path)
        
        if extended:
            extension_start = line.query_start
            extension_end = line.query_length - line.query_end

            # Also add 10% of the extension to each side
            if extension_start > 0:
                extension_start += extension_start // 10
            if extension_end > 0:
                extension_end += extension_end // 10

            path_start = line.path_start - extension_start
            path_end = line.path_end + extension_end
            
            if path_start < 0:
                path_start = 0
            if path_end > line.path_length:
                path_end = line.path_length

            ref = path_sequence[path_start:path_end]
            query = fastafile.fetch(line.query_name)
            wfa_alignment(aln, line, ref, query, path_start, output, extended)

        else:
            ref = path_sequence[line.path_start:line.path_end]
            query = fastafile.fetch(line.query_name, line.query_start, line.query_end)
            wfa_alignment([], line, ref, query, 0, output, False)

    fastafile.close()

    if extended:
        filter_duplicates(aln)
        write_alignments(aln, output)


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf', metavar='GAF', help='GAF File')
    arg('graph', metavar='GFA', help='Input GFA file')
    arg('fasta', metavar='FASTA', help='Input FASTA file')
    arg('-o', '--output', default=None, help='Output file. If omitted, use standard output.')
    arg("--ext", action='store_true', help="Extend the aligned path (use for Minigraph alignments).")
 
def main(args):
    run_realign(**vars(args))
