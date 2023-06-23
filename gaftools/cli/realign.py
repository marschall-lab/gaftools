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
from pywfa.align import WavefrontAligner
from collections import namedtuple, defaultdict

logger = logging.getLogger(__name__)


def run(gaf, graph, fasta):
    timers = StageTimer()
    
    realign_gaf(gaf, graph, fasta)
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

def parse_tag(s):
    name, type_id, value = s.split(':')
    assert len(name) == 2
    if type_id == 'i':
        return name, int(value)
    elif type_id == 'Z':
        return name, value
    else:
        assert False


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


def wfa_alignment(gaf_line, ref, query):

    match = 0
    cigar_len = 0
    
    aligner = WavefrontAligner(ref)
    res = aligner(query, clip_cigar=False)
    
    #0 = match
    for op_type, op_len in res.cigartuples:
        if op_type == 0:
            match += op_len
        cigar_len += op_len
        
    cigar = aligner.cigarstring.replace("M", "=")
        
    #Write the alignment back to the GAF
    sys.stdout.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" %(gaf_line.query_name, gaf_line.query_length, 
                     gaf_line.query_start, gaf_line.query_end, gaf_line.strand, gaf_line.path, gaf_line.path_length, 
                     gaf_line.path_start, gaf_line.path_end, match, cigar_len, gaf_line.mapping_quality))
        
    if gaf_line.is_primary:
        sys.stdout.write("\t%s" %gaf_line.is_primary[0])
    sys.stdout.write("\tcg:Z:%s\n" %cigar)
    



def realign_gaf(gaf, graph, fasta):
    """
        Uses pyWFA (https://github.com/kcleal/pywfa)
    """

    fastafile = pysam.FastaFile(fasta)
    nodes = parse_gfa(graph, with_sequence=True)


    #gaf_lines = []
    for cnt, line in enumerate(parse_gaf(gaf)):
        path_sequence = get_path(nodes, line.path)
        ref = path_sequence[line.path_start:line.path_end]
        query = fastafile.fetch(line.query_name, line.query_start, line.query_end)
        wfa_alignment(line, ref, query)

        #if cnt == 1000:
        #    print(cnt)
        #    break
        #gaf_lines.append(line)
    
    
    """for line in gaf_lines:
        path_sequence = get_path(nodes, line.path)
        ref = path_sequence[line.path_start:line.path_end]
        query = fastafile.fetch(line.query_name, line.query_start, line.query_end)

        wfa_alignment(line, ref, query)"""

    fastafile.close()


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf', metavar='GAF', help='GAF File')
    arg('graph', metavar='GFA', help='Input GFA file')
    arg('fasta', metavar='FASTA', help='Input FASTA file')
    #arg('-o', '--output', default=sys.stdout,
    #    help='Output GFA file. If omitted, use standard output.')
 
def main(args):
    run(**vars(args))
