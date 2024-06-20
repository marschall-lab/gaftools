"""
Tests for gaf.py.
"""

from gaftools.gaf import GAF


def test_parse_gaf_tags():
    # gaftools realign alignments.gaf smallgraph.gfa reads.fa

    gaf_file = open("tests/data/alignments-graphaligner.gaf", "r")
    gaf = GAF("tests/data/alignments-graphaligner.gaf")
    gaf_file2 = list(gaf.read_file())

    for cnt, line in enumerate(gaf_file):
        fields = line.rstrip().split("\t")

        for t in gaf_file2[cnt].tags.keys():
            matching = list(filter(lambda x: t in x, fields))
            assert len(matching) > 0

    gaf_file.close()
