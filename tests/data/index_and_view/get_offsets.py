import argparse
from pysam import libcbgzf
import gaftools.utils as utils


def run(fname):
    if utils.is_file_gzipped(fname):
        gaf_file = libcbgzf.BGZFile(fname, "rb")
        print("File is gzipped")
    else:
        gaf_file = open(fname, "rt")
    count = 0
    while True:
        offset = gaf_file.tell()
        mapping = gaf_file.readline()
        if not mapping:
            break
        print(f"Line {count} -> {offset}")
        count += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="get_offsets.py")
    parser.add_argument("-file", required=True, help="File name.")

    options = parser.parse_args()
    run(options.file)
