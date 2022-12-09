import gzip
import logging
from collections import defaultdict
from typing import Optional, DefaultDict

import pyfaidx
from dataclasses import dataclass


class InvalidRegion(Exception):
    pass


def detect_file_format(path):
    """
    Detect file format and return 'BAM', 'CRAM', 'VCF' or None. None indicates an
    unrecognized file format.

    'VCF' is returned for both uncompressed and compressed VCFs (.vcf and .vcf.gz).
    """
    with open(path, "rb") as f:
        first_bytes = f.read(16)
        if first_bytes.startswith(b"CRAM"):
            return "CRAM"
        if first_bytes.startswith(b"##fileformat=VCF"):
            return "VCF"

    gzip_header = b"\037\213"
    if first_bytes.startswith(gzip_header):
        with gzip.GzipFile(path, "rb") as f:
            first_bytes = f.read(16)
            if first_bytes.startswith(b"BAM\1"):
                return "BAM"
            elif first_bytes.startswith(b"##fileformat=VCF"):
                return "VCF"

    return None

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

@dataclass
class Region:
    chromosome: str
    start: int
    end: Optional[int]

    def __repr__(self):
        return f'Region("{self.chromosome}", {self.start}, {self.end})'

    @staticmethod
    def parse(spec: str):
        """
        >>> Region.parse("chr1")
        Region("chr1", 0, None)
        >>> Region.parse("chr1:")
        Region("chr1", 0, None)
        >>> Region.parse("chr1:101")
        Region("chr1", 100, None)
        >>> Region.parse("chr1:101-")
        Region("chr1", 100, None)
        >>> Region.parse("chr1:101-200")
        Region("chr1", 100, 200)
        >>> Region.parse("chr1:101:200")  # for backwards compatibility
        Region("chr1", 100, 200)
        """
        parts = spec.split(":", maxsplit=1)
        chromosome = parts[0]
        if len(parts) == 1 or not parts[1]:
            start, end = 0, None
        else:
            try:
                sep = ":" if ":" in parts[1] else "-"
                start_end = parts[1].split(sep, maxsplit=1)
                start = int(start_end[0]) - 1
                if len(start_end) == 1 or not start_end[1]:
                    end = None
                else:
                    end = int(start_end[1])
                    if end <= start:
                        raise InvalidRegion("end is before start in specified region")
            except ValueError:
                raise InvalidRegion("Region must be specified as chrom[:start[-end]])") from None
        return Region(chromosome, start, end)

