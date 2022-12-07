"""
Sorting GAF and GFA files
"""

import logging
import sys
import platform

from gaftools import __version__
from gaftools.cli import log_memory_usage
from gaftools.cli import CommandLineError


logger = logging.getLogger(__name__)

def run(gaf_file, gfa_file, output=sys.stdout):
    if gaf_file and gfa_file:
        print("Error: Please input only the GAF or the GFA file")
        print()
        return
    elif not gaf_file and not gfa_file:
        print("Error: Please input a GAF or a GFA file to be sorted")
        print()
        return

    if (gaf_file):
        gaf_sort(gaf_file, output)
    else:
        gfa_sort(gfa_file, output)


def gfa_sort(gfa_path, out_path = None, return_list = False):
    '''This function sorts the given gfa file based on the contig name and start position within the
    contig. Note that it only sorts S lines and leaves the others.
    This can be called from the command line or from another funtion by providing "True" to
    return_list argument.
    '''

    import functools
    import gzip

    gfa_lines = []

    gz_flag = gfa_path[-2:] == "gz"
    if gz_flag:
        gfa_file = gzip.open(gfa_path,"r")
    else:
        gfa_file = open(gfa_path,"r")
    for line_count, line in enumerate(gfa_file):
        gfa_lines.append(line)
    #print("Read", line_count, "lines")
    gfa_file.close()
    gfa_s = []
    gfa_other = []
    for line in gfa_lines:
        val = line.rstrip().split('\t')
        if val[0] == "S":
            gfa_s.append(val)
        else:
            gfa_other.append(val)
    
    gfa_s.sort(key=functools.cmp_to_key(compare_gfa))
    
    if return_list:
       return gfa_s

    elif out_path:
        with open(out_path, "w", encoding='utf-8') as out_file:
            for line_count, line in enumerate(gfa_s):
                out_file.write('\t'.join(line) + '\n')
        
            for line_count2, line in enumerate(gfa_other):
                out_file.write('\t'.join(line) + '\n')

        #print("Wrote", line_count + line_count2, "lines to", out_file)
    
    else:
        for line_count, line in enumerate(gfa_s):
             print('\t'.join(line) + '\n')
        
        for line_count2, line in enumerate(gfa_other):
            print('\t'.join(line) + '\n')
 
    return True


def gaf_sort(gaf_path, out_path = None):
    '''This function sorts a gaf file (mappings to a pangenome graph) in sstable coordinate system based on 1)Contig name 2)Start
    position of the contig's mapping loci.
    '''

    import functools
    import gzip

    gaf_lines = []
    
    if is_file_gzipped(gaf_path):
       open_gaf = gzip.open
       is_gzipped = True
    else:
        open_gaf = open

    with open_gaf(gaf_path, "rt") as gaf_file:
        for line_count, mapping in enumerate(gaf_file):
            val = mapping.rstrip().split('\t')
            gaf_lines.append(val)

    gaf_lines.sort(key=functools.cmp_to_key(compare_gaf))

    if out_path:
        with open(out_path, "w") as out_file:
            for line_count, line in enumerate(gaf_lines):
                out_file.write('\t'.join(line) + '\n') 
    else:
        for line_count, line in enumerate(gaf_lines):
             print('\t'.join(line) + '\n')

    return True



def compare_gfa(ln1, ln2):
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



def compare_gaf(ln1, ln2):

    if ln1[5][0] == ">" or ln1[5][0] == "<":
        chr1 = ln1[5][1:].lower()
    else:
        chr1 = ln1[5].lower()


    if ln2[5][0] == ">" or ln2[5][0] == "<":
        chr2 = ln2[5][1:].lower()
    else:
        chr2 = ln2[5].lower()

    start1 = ln1[7]
    start2 = ln2[7]

    if chr1 == chr2:
        if start1 < start2:
            return -1
        else:
            return 1
    if chr1 < chr2:
        return -1
    else:
        return 1


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'


def add_arguments(parser):
    arg = parser.add_argument

    arg('--gaf', dest = 'gaf_file', metavar='GAF', default = None, required = False, help='GAF File whose coordinates have to be changed')
    arg('--gfa', dest = 'gfa_file', metavar='GFA', default = None, required = False, help='Input GFA file to conver the coordinates')
    arg('-o', '--output', default=sys.stdout, help='Output GAF file. If omitted, use standard output.')


def validate(args, parser):
    return True


def main(args):
    run(**vars(args))
