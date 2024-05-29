import re
import glob
from heapq import merge
import functools
import gzip
import os

complement = str.maketrans('ACGT', 'TGCA')
tag_regex = r"^[A-Za-z][A-Za-z][:][AifZHB][:][ !-~]*$"

types_regex = {
    "A": r"^[!-~]$",
    "i": r"^[-+]?[0-9]+$",
    "f": r"^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$",
    "Z": r"^[ !-~]*$",
    "H": r"^([0-9A-F][0-9A-F])*$",
    "B": r"^[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)*$"
    }


def rev_comp(seq):
    return seq[::-1].translate(complement)


def is_correct_tag(tag):
    # first check if tag follows the scheme two_letters:{AifZHB}:value
    if not re.match(tag_regex, tag):
        return False
    name, tag_type, value = tag.split(":")
    if not re.match(types_regex[tag_type], value):
        return False
    return True


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'


def search_intervals(intervals, query_start, query_end, start, end):
    '''Given the start-end coordinates in the GFA file for the given contig (SO, SO+LN), it
    searches for the given (query_start, query_end) matches. (query_start, query_end) is the start
    and end location of a mapping in the gaf file.
    '''

    if start <= end:    
        mid = start + (end - start) // 2
        if query_end <= intervals[mid].start:
            return search_intervals(intervals, query_start, query_end, start, mid - 1)
        elif query_start >= intervals[mid].end:
            return search_intervals(intervals, query_start, query_end, mid + 1, end)
        else:
            return start, end

    return -1, -1


def reverse_cigar(cg):
    import itertools
    cg = cg[5:]
    all_cigars = ["".join(x) for _, x in itertools.groupby(cg, key=str.isdigit)]
    new_cigar = "cg:Z:"
    for i in range(len(all_cigars), 0, -2):
        new_cigar += str(all_cigars[i-2]) + str(all_cigars[i-1])
    return new_cigar


def is_unstable(line):
    g = list(filter(None, re.split('(>)|(<)', line.rstrip().split('\t')[5])))
    
    if len(g) == 1:
        return False
    
    if ":" in g[1]:
        return False
    
    return True

def gfa_sort_basic(gfa_path):
    '''This function sorts the given gfa file based on the contig name and start position within the
    contig. Note that it only sorts S lines and leaves the others.
    This can be called from the command line or from another funtion by providing "True" to the
    return_list argument.
    '''

    gfa_lines = []
    path = "part*.gfa"
    chunk_size = 250000
    chunk_id = 1

    if is_file_gzipped(gfa_path):
       open_gfa = gzip.open
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
            
                f_out.close()
                gfa_lines = []
                chunk_id += 1
                f_out = open('part_{}.gfa'.format(chunk_id), 'w')


        if gfa_lines:
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