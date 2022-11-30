"""
Statistics of a GAF File
"""

import sys

def run(gaf_path, cigar_stat, output=sys.stdout):
    gaf_stat(gaf_path, cigar_stat)


def gaf_stat(gaf_path, cigar_stat=False):
    '''This function outputs some statistics of a GAF file. If you run with "--cigar" option, then
    cigar related statistics are also output. This increases the running time more than 10 folds
    because I need to iterate through the cigar of each alignment making it run in O(n^2)
    '''

    import gzip
    import itertools

    gaf_lines = []
    
    if is_file_gzipped(gaf_path):
       open_gaf = gzip.open
       is_gzipped = True
    else:
        open_gaf = open
    
    read_names = set()
    total_aligned_bases = 0
    total_mapq = 0
    total_primary = 0
    total_secondary = 0

    if cigar_stat:
        total_del = 0
        total_ins = 0
        total_x = 0
        total_del_large = 0
        total_ins_large = 0
        total_x_large = 0
        total_match = 0
        total_match_large = 0
        total_perfect = 0
        
    with open_gaf(gaf_path, "rt") as gaf_file:
        for mapping in gaf_file:
            val = mapping.rstrip().split('\t')
            gaf_lines.append(val)

    #print("Done reading")
    for line_count, val in enumerate(gaf_lines):
        hashed_readname = hash(val[0])
        read_names.add(hashed_readname)
        total_aligned_bases += int(val[9])
        total_mapq += int(val[11])
        is_primary = [k for k in val if k.startswith("tp:A:")][0][5:]
        if is_primary != "P":
            total_secondary += 1
        else:
            total_primary += 1
        
        if cigar_stat:
            '''Cigar string analysis'''
            cigar = [k for k in val if k.startswith("cg:Z:")][0][5:]
            all_cigars = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
            if len(all_cigars) == 2:
                total_perfect += 1
            #print(all_cigars)
            for cnt in range(0, len(all_cigars)-1, 2):
                if all_cigars[cnt+1] == "D":
                    total_del += 1
                    if int(all_cigars[cnt]) >= 50:
                        total_del_large += 1
                elif all_cigars[cnt+1] == "I":
                    total_ins += 1
                    if int(all_cigars[cnt]) >= 50:
                        total_ins_large += 1
                elif all_cigars[cnt+1] == "X":
                    total_x += 1
                    if int(all_cigars[cnt]) >= 50:
                        total_x_large += 1
                elif all_cigars[cnt+1] == "=":
                    total_match += 1
                    if int(all_cigars[cnt]) >= 50:
                         total_match_large += 1
    
    print()
    print("Total alignments:", line_count + 1)
    print("\tPrimary:", total_primary)
    print("\tSecondary:", total_secondary)
    print("Reads with at least one alignment:", len(read_names))
    print("Total aligned bases:", str(total_aligned_bases))
    print("Average mapping quality:", round((total_mapq/line_count), 1))
    
    if cigar_stat:
        print("Cigar string statistics:\n\tTotal deletion regions: %d (%d >50bps)\n\tTotal insertion regions: %d (%d >50bps)\n\tTotal substitution regions: %d (%d >50bps)\n\tTotal match regions: %d (%d >50bps)" % 
          (total_del, total_del_large, total_ins, total_ins_large, total_x, total_x_large,
           total_match, total_match_large))
        
        print("Total perfect alignments (exact match):", total_perfect)
    
    print()
    
    return True


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'


def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_path', metavar='GAF', help='Input GAF file')
    arg('--cigar', dest='cigar_stat', default=False, action='store_true', help='Outputs cigar related statistics (takes a long)')


def validate(args, parser):
    return True


def main(args):
    run(**vars(args))

