import re
import logging

from gaftools import __version__
from collections import namedtuple, defaultdict
from gaftools.utils import is_file_gzipped

logger = logging.getLogger(__name__)

complement = str.maketrans('ACGT', 'TGCA')

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


GafLine = namedtuple("GafLine", "query_name query_length query_start query_end strand path path_length path_start path_end residue_matches alignment_block_length mapping_quality is_primary")

def parse_gaf(filename):
    for line in open(filename):
        fields = line.split('\t')
        print(line)
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
            #cigar = [k for k in fields if k.startswith("cg:Z:")][-1][5:].strip()
        )
    return


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



def gfa_sort_basic(gfa_path):
    '''This function sorts the given gfa file based on the contig name and start position within the
    contig. Note that it only sorts S lines and leaves the others.
    This can be called from the command line or from another funtion by providing "True" to the
    return_list argument.
    '''

    import glob
    from heapq import merge
    import functools
    import gzip
    import os


    gfa_lines = []
    path = "part*.gfa"
    chunk_size = 250000
    chunk_id = 1

    if is_file_gzipped(gfa_path):
       open_gfa = gzip.open
       is_gzipped = True
    else:
        open_gfa = open
    
    with open_gfa(gfa_path, 'rt') as gfa_file:
        f_out = open('part_{}.gfa'.format(chunk_id), 'w')
        
        for line_num, mapping in enumerate(gfa_file, 1):
            val = mapping.rstrip().split('\t')
            gfa_lines.append(val)
        
            if not line_num % chunk_size:
                gfa_lines.sort(key=functools.cmp_to_key(compare_gfa_lines))
                
                for line_count, line in enumerate(gfa_lines):
                    f_out.write('\t'.join(line) + '\n') 
            
                logger.info('INFO: Splitting %d' %chunk_id)
                f_out.close()
                gfa_lines = []
                chunk_id += 1
                f_out = open('part_{}.gfa'.format(chunk_id), 'w')


        if gfa_lines:
            logger.info('INFO: Splitting %d' %chunk_id)
            gfa_lines.sort(key=functools.cmp_to_key(compare_gfa_lines))
            for line_count, line in enumerate(gfa_lines):
                f_out.write('\t'.join(line) + '\n') 
            f_out.close()
            gfa_lines = []

    chunks = []
    for filename in glob.glob(path):
        chunks += [open(filename, 'r')]
     
    gfa_s = []
    tmp = merge(*chunks, key=functools.cmp_to_key(compare_gfa_merge))
    for i in tmp:
        if i[0] == "S":
            gfa_s.append(i.rstrip().split('\t'))
        
    for part_file in glob.glob(path):
        if os.path.isfile(part_file):
            os.remove(part_file)

    return gfa_s
 

def compare_gfa_lines(ln1, ln2):
    
    if not ln1[0] == "S" and not ln2[0] == "S":
        #If both are not S lines, then leave it
        return -1
    elif ln1[0] == "S" and not ln2[0] == "S":
        return -1
    elif not ln1[0] == "S" and ln2[0] == "S":
        return 1

    chr1 = [k for k in ln1 if k.startswith("SN:Z:")][0][5:]
    chr2 = [k for k in ln2 if k.startswith("SN:Z:")][0][5:]
    start1 = int([k for k in ln1 if k.startswith("SO:i:")][0][5:])
    start2 = int([k for k in ln2 if k.startswith("SO:i:")][0][5:])

    if chr1 == chr2:
        if start1 < start2:
            return -1
        else:
            return 1
    if chr1 < chr2:
        return -1
    else:
        return 1


def compare_gfa_merge(ln1, ln2):

    ln1 = ln1.rstrip().split('\t')
    ln2 = ln2.rstrip().split('\t')
    return compare_gfa_lines(ln1, ln2)
    
