"""
Sort GFA File
"""

import sys

def run(gfa_path, output=sys.stdout):
    gfa_sort(gfa_path, output)

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
    
    gfa_s.sort(key=functools.cmp_to_key(compare))
    
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


def compare(ln1, ln2):
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


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gfa_path', metavar='GFA', help='Input GFA file to conver the coordinates')
    arg('-o', '--output', default=sys.stdout,
        help='Output sorted GFA file. If omitted, use standard output.')
    
# fmt: on
def validate(args, parser):
    return True

def main(args):
    run(**vars(args))
