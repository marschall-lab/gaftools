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


def run(gaf, graph, fasta, output=sys.stdout):
    timers = StageTimer()
    
    realign_gaf(gaf, graph, fasta, output)
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


complement = str.maketrans('ACGT', 'TGCA')

class Node:
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
    edges = defaultdict(list)

    for nr, line in enumerate(open(gfa_filename)):
        fields = line.split('\t')
        if fields[0] == 'S':
            name = fields[1]
            tags = dict(parse_tag(s) for s in fields[3:])
            sequence = None
            if with_sequence and (fields[2] != '*'):
                sequence = fields[2]
            nodes[name] = Node(name,tags,sequence)
        elif fields[0] == 'L':
            from_node = fields[1]
            from_dir = fields[2]
            to_node = fields[3]
            to_dir = fields[4]
            overlap = fields[5]
            tags = dict(parse_tag(s) for s in fields[6:])
            e = Edge(from_node,from_dir,to_node,to_dir,overlap, tags)
            edges[(from_node,to_node)].append(e)

    return nodes, edges


GafLine = namedtuple("GafLine", "query_name query_length query_start query_end strand path path_length path_start path_end residue_matches alignment_block_length mapping_quality is_primary cigar")

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
            cigar = [k for k in fields if k.startswith("cg:Z:")][0][5:].strip()
        )
    return


def get_path(nodes, path):
    l = []
    for s in re.findall('[><][^><]+', path):
        node = nodes[s[1:]]
        #print(node.name, len(node.sequence))
        if s[0] == '>':
            l.append(node.sequence)
        elif s[0] == '<':
            l.append(node.sequence[::-1].translate(complement))
        else:
            assert False
    return ''.join(l)


def realign_gaf(gaf, graph, fasta, output):
    """
        Uses pyWFA (https://github.com/kcleal/pywfa)
    """
    
    #print('Reading', graph)
    fastafile = pysam.FastaFile(fasta)
    nodes, edges = parse_gfa(graph, with_sequence=True)
    
    f = open(output, "w")
    
    for cnt, gaf_line in enumerate(parse_gaf(gaf)):
        #print('RESULT GAF', gaf_line.query_name, gaf_line.query_start, gaf_line.query_end, gaf_line.path, gaf_line.strand)

        path_sequence = get_path(nodes, gaf_line.path)
        ref = path_sequence[gaf_line.path_start:gaf_line.path_end]
        query = fastafile.fetch(gaf_line.query_name, gaf_line.query_start, gaf_line.query_end)
        
        aligner = WavefrontAligner(ref)
        res = aligner(query, clip_cigar=False)
        
        match_new = 0
        cigar_len_new = 0
        for op_type, op_len in res.cigartuples:
            if op_type == 0:
                match_new += op_len
            cigar_len_new += op_len
        
        cigar = aligner.cigarstring.replace("M", "=")
        
        #Write the alignment back to the GAF
        f.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" %(gaf_line.query_name, gaf_line.query_length, 
                     gaf_line.query_start, gaf_line.query_end, gaf_line.strand, gaf_line.path, gaf_line.path_length, 
                     gaf_line.path_start, gaf_line.path_end, match_new, cigar_len_new, gaf_line.mapping_quality))
        
        if gaf_line.is_primary:
            f.write("\t%s" %gaf_line.is_primary[0])
        f.write("\tcg:Z:%s\n" %cigar)
        
        """if cnt == 10000:
            print(cnt)
            break"""
    
    f.close()
    fastafile.close()


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf', metavar='GAF', help='GAF File')
    arg('graph', metavar='GFA', help='Input GFA file')
    arg('fasta', metavar='FASTA', help='Input FASTA file')
    arg('-o', '--output', default=sys.stdout,
        help='Output GFA file. If omitted, use standard output.')
 
# fmt: on
def validate(args, parser):
    """if not args.gaf:
        parser.error("GAF file has to be provided")
    if args.graph:
        parser.error("Graph file (GFA) has to be provided")
    if args.fasta:
        parser.error("FASTA file has to be provided")
    """
def main(args):
    run(**vars(args))
