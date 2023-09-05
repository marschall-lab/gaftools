import re

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
    a = line.rstrip().split('\t')[5]
    g = list(filter(None, re.split('(>)|(<)', a)))
    if len(g) == 1:
        return False
    if ":" in g[1]:
        return False
    return True
