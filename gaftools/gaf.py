import re
import gzip
import logging
import gaftools.utils as utils
from gaftools import __version__
from collections import defaultdict

logger = logging.getLogger(__name__)

complement = str.maketrans('ACGT', 'TGCA')

class Alignment:
    def __init__(self, query_name, query_length, query_start, query_end, strand, path, path_length, path_start,
                 path_end, residue_matches, alignment_block_length, mapping_quality, is_primary,
                 cigar, cigar_length=None, score=None, tags=None):

        self.query_name = query_name
        self.query_length = query_length
        self.query_start = query_start
        self.query_end = query_end
        self.strand = strand
        self.path = path
        self.path_length = path_length
        self.path_start = path_start
        self.path_end = path_end
        self.residue_matches = residue_matches
        self.alignment_block_length = alignment_block_length
        self.cigar_length = cigar_length
        self.score = score
        self.cigar = cigar
        self.mapping_quality = mapping_quality
        self.is_primary = is_primary
        self.duplicate = False
        self.tags = tags


class Read:
    def __init__(self, query_name, read_length, map_ratio, seq_identity):
        self.rname = query_name
        self.length = read_length
        self.aln_count = 1
        self.highest_map_ratio = map_ratio
        self.total_map_ratio = map_ratio
        self.highest_seq_identity = seq_identity
        self.total_seq_identity = seq_identity


def parse_gaf(filename):
   
    gz_flag = False
    if utils.is_file_gzipped(filename):
        gaf_file = gzip.open(filename,"r")
        gz_flag = True
    else:
        gaf_file = open(filename,"r") 

    for line in gaf_file:
        if not gz_flag:
            fields = line.rstrip().split('\t')
        else:
            fields = line.decode("utf-8").rstrip().split('\t')
        
        #If the query name has spaces (e.g., GraphAligner), we get rid of the segment after the space
        query_name = fields[0].split(' ')[0]
        query_length = int(fields[1])
        query_start = int(fields[2])
        query_end = int(fields[3])
        strand = fields[4]
        path = fields[5]
        path_length = int(fields[6])
        path_start = int(fields[7])
        path_end = int(fields[8])
        residue_matches = int(fields[9])
        alignment_block_length = int(fields[10])
        mapping_quality = int(fields[11])
        is_primary = True
        cigar = ""

        #Check if there are additional tags
        tags = {}
        for k in fields:
            if re.match("[A-Za-z][A-Za-z0-9]:[AifZHB]:[A-Za-z0-9]+", k):
                pattern = re.findall(r"([A-Za-z][A-Za-z0-9]:[AifZHB]:)[A-Za-z0-9]+", k)[0]
                
                if pattern == "cg:Z:":
                    val = re.findall(r"[A-Za-z][A-Za-z0-9]:[AifZHB]:([A-Za-z0-9=]+)", k)[0]
                    cigar = val
                else:
                    val = re.findall(r"[A-Za-z][A-Za-z0-9]:[AifZHB]:([A-Za-z0-9]+)", k)[0]
                    if pattern not in tags:
                        tags[pattern] = val

                    if pattern == "tp:A" and (val != "P" or val != "p"):
                        is_primary = False
        
        yield Alignment(query_name, query_length, query_start, query_end, strand, path,
                        path_length, path_start, path_end, residue_matches, alignment_block_length,
                        mapping_quality, is_primary, cigar, tags=tags)
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



