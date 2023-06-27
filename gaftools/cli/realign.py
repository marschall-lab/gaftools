"""
Realign GAF file using wavefront alignment algorithm (WFA)
"""

import logging
import sys
import re
import pysam

from gaftools import __version__
from gaftools.cli import log_memory_usage
from gaftools.cli import CommandLineError
from gaftools.timer import StageTimer
from pywfa.align import (WavefrontAligner, cigartuples_to_str)
from collections import namedtuple, defaultdict

logger = logging.getLogger(__name__)


def run(gaf, graph, fasta, extended):
    timers = StageTimer()
    
    realign_gaf(gaf, graph, fasta, extended)
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


complement = str.maketrans('ACGT', 'TGCA')

"""class Node:
    def __init__(self, name, tags, sequence=None):
        self.name = name
        self.tags = tags
        self.sequence = sequence

class Edge:
    def __init__(self, from_node, from_dir, to_node, to_dir, overlap, tags):
        self.from_node = from_node
        self.from_dir = from_dir
        self.to_node = to_node
        self.to_dir = to_dir
        self.overlap = overlap
        self.tags = tags
"""

class Alignment:
    def __init__(self, query_length, query_start, query_end, strand, path, path_length, path_start,
                 path_end, residue_matches, cigar_length, score, cigar):
        self.query_length = query_length
        self.query_start = query_start
        self.query_end = query_end
        self.strand = strand
        self.path = path
        self.path_length = path_length
        self.path_start = path_start
        self.path_end = path_end
        self.residue_matches = residue_matches
        self.cigar_length = cigar_length
        self.score = score
        self.cigar = cigar
        self.duplicate = False


def parse_tag(s):
    name, type_id, value = s.split(':')
    assert len(name) == 2
    if type_id == 'i':
        return name, int(value)
    elif type_id == 'Z':
        return name, value
    else:
        assert False


def compare_aln(ln1, ln2):

    if ln1.query_start < ln2.query_start:
        return -1
    elif ln1.query_start == ln2.query_start:
        if ln1.query_end < ln2.query_end:
            return 1
        else:
            return -1
    else:
        return 1


def parse_gfa(gfa_filename, with_sequence=False):
    nodes = {}

    for nr, line in enumerate(open(gfa_filename)):
        fields = line.split('\t')
        if fields[0] == 'S':
            name = fields[1]
            #tags = dict(parse_tag(s) for s in fields[3:])
            sequence = None
            if with_sequence and (fields[2] != '*'):
                sequence = fields[2]
            #nodes[name] = Node(name,tags,sequence)
            nodes[name] = sequence
    return nodes


GafLine = namedtuple("GafLine", "query_name query_length query_start query_end strand path path_length path_start path_end residue_matches alignment_block_length mapping_quality is_primary")

def parse_gaf(filename):
    for line in open(filename):
        fields = line.split('\t')
 
        yield GafLine(
            #If the query name has spaces (e.g., GraphAligner), we get rid of the segment after the space
            query_name = fields[0].split(' ')[0],
            query_length = int(fields[1]),
            query_start = int(fields[2]),
            query_end = int(fields[3]),
            strand = fields[4],
            path = fields[5],
            path_length = int(fields[6]),
            path_start = int(fields[7]),
            path_end = int(fields[8]),
            residue_matches = int(fields[9]),
            alignment_block_length = int(fields[10]),
            mapping_quality = int(fields[11]),
            is_primary = [k for k in fields if k.startswith("tp:A:") or None],
            #cigar = [k for k in fields if k.startswith("cg:Z:")][0][5:].strip()
        )
    return


def get_path(nodes, path):
    l = []
    for s in re.findall('[><][^><]+', path):
        node_seq = nodes[s[1:]]
        #print(node.name, len(node.sequence))
        if s[0] == '>':
            l.append(node_seq)
        elif s[0] == '<':
            l.append(node_seq[::-1].translate(complement))
        else:
            assert False
    return ''.join(l)




def overlap_ratio(x_start,x_end, y_start, y_end):

    overlap = max(0, min(x_end, y_end) - max(x_start, y_start))
    total_length = x_end - x_start + y_end - y_start
    x_length = x_end - x_start
    y_length = y_end - y_start

    return max(2 * (overlap / total_length) , (overlap / x_length), (overlap/y_length))


def filter_duplicates(aln):
    import functools

    #print("Filtering...")
    for k in aln.keys():
        aln[k].sort(key=functools.cmp_to_key(compare_aln))
     
    for read_name, mappings in aln.items():
        if len(mappings) == 1:
            continue
        #print(read_name)
        for cnt, line in enumerate(mappings):
            if line.duplicate == True:
                continue
            #print("Query", line.query_start, line.query_end, line.score)
            for cnt2, line2 in enumerate(mappings):
                if cnt == cnt2 or line2.duplicate == True:
                    continue
                
                sim = overlap_ratio(line.query_start, line.query_end, line2.query_start, line2.query_end)
                
                #print("-", line2.query_start, line2.query_end, line2.score, "sim =", sim) 
                if sim > 0.8:
                    if(line.score > line2.score):
                        mappings[cnt2].duplicate = True
                        #print("......Duplicate")
                    else:
                        mappings[cnt].duplicate = True
                        #print("......Duplicate, breaking")
                        break
        """print()
        for line in mappings:
            print(line.query_start, line.query_end, line.duplicate)
                    
        print()"""


def write_alignments(aln):
    
    #Write the alignment back to the GAF
    for read_name, mappings in aln.items():
        for gaf_line in mappings:
            if not gaf_line.duplicate:
                sys.stdout.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" %(read_name, gaf_line.query_length, 
                     gaf_line.query_start, gaf_line.query_end, gaf_line.strand, gaf_line.path, gaf_line.path_length, 
                     gaf_line.path_start, gaf_line.path_end, gaf_line.residue_matches, gaf_line.cigar_length, 60))
        
                #if gaf_line.is_primary:
                #    sys.stdout.write("\t%s" %gaf_line.is_primary[0])
                sys.stdout.write("\tcg:Z:%s\n" %gaf_line.cigar)


def wfa_alignment(aln, gaf_line, ref, query, path_start, extended):

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
            print("ERRRORRR")
            print(op_type, op_len)
        cigar_len += op_len
    
    """print("Match = ", match, " ins = ", ins, " del = ", deletion, " mismatch ", mismatch, " soft-clip = ", soft_clip)
    print("Query check: ", match + ins + soft_clip + mismatch)
    #print(cigartuples_to_str(res.cigartuples))
    #print(res.cigartuples)
    
    print("Query start-end = ", res.text_start, res.text_end)
    print("Reference start-end = ", res.pattern_start, res.pattern_end)
    print("Score = ", res.score)
    print("Cigar len = ", cigar_len, " Match = ", match)
    #print(cigar)
    print()
    """
    

    if extended:
        if gaf_line.query_name in aln:
            aln[gaf_line.query_name].append(Alignment(gaf_line.query_length, res.text_start,
                                               res.text_end, gaf_line.strand, gaf_line.path,
                                               gaf_line.path_length, path_start + res.pattern_start,
                                               path_start + res.pattern_end, match, cigar_len,
                                               res.score, cigar))
        else:
            aln[gaf_line.query_name] = [Alignment(gaf_line.query_length, res.text_start,
                                               res.text_end, gaf_line.strand, gaf_line.path,
                                               gaf_line.path_length, path_start + res.pattern_start,
                                               path_start + res.pattern_end, match, cigar_len,
                                               res.score, cigar)]
    else:
        cigar = aligner.cigarstring.replace("M", "=")

        #Write the alignment back to the GAF
        sys.stdout.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" %(gaf_line.query_name, gaf_line.query_length, 
                     gaf_line.query_start, gaf_line.query_end, gaf_line.strand, gaf_line.path, gaf_line.path_length, 
                     gaf_line.path_start, gaf_line.path_end, match, cigar_len, gaf_line.mapping_quality))

        if gaf_line.is_primary:
            sys.stdout.write("\t%s" %gaf_line.is_primary[0])
        sys.stdout.write("\tcg:Z:%s\n" %cigar)


def realign_gaf(gaf, graph, fasta, extended):
    """
        Uses pyWFA (https://github.com/kcleal/pywfa)
    """
    

    fastafile = pysam.FastaFile(fasta)
    nodes = parse_gfa(graph, with_sequence=True)
    
    #ref = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
    #query = "ACCCTCGCAAGCGTTGGAGAATG"

    aln = {}
    for cnt, line in enumerate(parse_gaf(gaf)):
        if not line.is_primary:
            continue

        path_sequence = get_path(nodes, line.path)

        if extended: 
            extension_start = line.query_start
            extension_end = line.query_length - line.query_end

            path_start = line.path_start - extension_start
            path_end = line.path_end + extension_end

            if path_start < 0:
                path_start = 0
            if path_end > line.path_length:
                path_end = line.path_length

            ref = path_sequence[path_start:path_end]
            query = fastafile.fetch(line.query_name)
            wfa_alignment(aln, line, ref, query, path_start, extended)

        else:
            ref = path_sequence[line.path_start:line.path_end]
            query = fastafile.fetch(line.query_name, line.query_start, line.query_end)
            wfa_alignment([], line, ref, query, 0, False)

        #if line.query_name == "bf9a7d6b-5595-4bf5-a48a-d93f28e195e1":
        #if line.query_name == "beda39f4-2976-4c0b-9cb2-9bf4ebdf0129":


        if cnt == 2000:
            break

    fastafile.close()

    if extended:
        filter_duplicates(aln)
        write_alignments(aln)


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf', metavar='GAF', help='GAF File')
    arg('graph', metavar='GFA', help='Input GFA file')
    arg('fasta', metavar='FASTA', help='Input FASTA file')
    arg("--extended", action='store_true', help="Realign the whole read instead of the alignment only\'o.")
    #arg('-o', '--output', default=sys.stdout,
    #    help='Output GFA file. If omitted, use standard output.')
 
def main(args):
    run(**vars(args))
