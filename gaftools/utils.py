import re
import tracemalloc
import linecache
import sys
from pysam import libcbgzf

DEFAULT_CHROMOSOME = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
    "chrM",
]

complement = str.maketrans("ACGT", "TGCA")
tag_regex = r"^[A-Za-z][A-Za-z][:][AifZHB][:][ !-~]*$"

types_regex = {
    "A": r"^[!-~]$",
    "i": r"^[-+]?[0-9]+$",
    "f": r"^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$",
    "Z": r"^[ !-~]*$",
    "H": r"^([0-9A-F][0-9A-F])*$",
    "B": r"^[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)*$",
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
        return inp.read(2) == b"\x1f\x8b"


def search_intervals(intervals, query_start, query_end, start, end):
    """Given the start-end coordinates in the GFA file for the given contig (SO, SO+LN), it
    searches for the given (query_start, query_end) matches. (query_start, query_end) is the start
    and end location of a mapping in the gaf file.
    """

    if start <= end:
        mid = start + (end - start) // 2
        if query_end <= int(intervals[mid].tags["SO"][1]):
            return search_intervals(intervals, query_start, query_end, start, mid - 1)
        elif query_start >= int(intervals[mid].tags["SO"][1]) + int(intervals[mid].tags["LN"][1]):
            return search_intervals(intervals, query_start, query_end, mid + 1, end)
        else:
            return start, end

    return -1, -1


def reverse_cigar(cg):
    import itertools

    all_cigars = ["".join(x) for _, x in itertools.groupby(cg, key=str.isdigit)]
    new_cigar = ""
    for i in range(len(all_cigars), 0, -2):
        new_cigar += str(all_cigars[i - 2]) + str(all_cigars[i - 1])
    return new_cigar


def display_top(snapshot, key_type="lineno", limit=3):
    snapshot = snapshot.filter_traces(
        (
            tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
            tracemalloc.Filter(False, "<unknown>"),
        )
    )
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = "/".join(frame.filename.split("/")[-2:])
        print("#%s: %s:%s: %.1f KiB" % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print("    %s" % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))


class FileWriter:
    def __init__(self, name: str):
        self.name = name
        self.writer = None
        self.compress = False
        if name is None:
            self.writer = sys.stdout
            self.compress = False
        elif name.endswith(".gz"):
            self.writer = libcbgzf.BGZFile(name, "wb")
            self.compress = True
        else:
            self.writer = open(name, "w")
            self.compress = False

    def write(self, line):
        self.writer.write(str.encode(line) if self.compress else line)

    def tell(self):
        if self.name is not None:
            return self.writer.tell()
        else:
            raise RuntimeError("Cannot get the file pointer position for standard output")

    def close(self):
        if self.name is not None:
            self.writer.close()
