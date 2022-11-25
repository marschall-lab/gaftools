"""
Sort GAF File
"""

import sys

def run(gaf_path, output=sys.stdout):
    gaf_sort(gaf_path, output)


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
    #print("Read", line_count, "lines")

    gaf_lines.sort(key=functools.cmp_to_key(compare))

    if out_path:
        with open(out_path, "w") as out_file:
            for line_count, line in enumerate(gaf_lines):
                out_file.write('\t'.join(line) + '\n') 
    else:
        for line_count, line in enumerate(gaf_lines):
             print('\t'.join(line) + '\n')
     
    return True


def compare(ln1, ln2):
    
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
    # Positional arguments
    arg('gaf_path', metavar='GAF', help='Input GAF file to be sorted')
    arg('-o', '--output', default=sys.stdout, help='The name of the output file (Write to stdout by default).')


def validate(args, parser):
    return True


def main(args):
    run(**vars(args))


