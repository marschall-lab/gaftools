"""
Tests for 'gaftools sort'
"""

from gaftools.cli.sort import run_sort
from gaftools.cli.order_gfa import run_order_gfa
from gaftools.cli import CommandLineError
from gaftools.gaf import GAF
import pytest


def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split("\t"))

    return [parse_line(l) for l in open(filename)]


@pytest.mark.parametrize("gaf_file", ["smallgraph.gaf", "smallgraph.gaf.gz"])
def test_sort_smallgraph(tmp_path, gaf_file):
    input_gaf = f"tests/data/sort/{gaf_file}"
    input_gfa = "tests/data/order_gfa/smallgraph-ordered.gfa"
    output = str(tmp_path) + "/output.gaf"
    run_sort(gfa=input_gfa, gaf=input_gaf, outgaf=output)
    gaf = GAF(output)
    gaf_lines = list(gaf.read_file())
    assert gaf_lines[0].query_name == "read_s3"
    assert gaf_lines[0].tags["bo:i:"] == "3"
    assert gaf_lines[0].tags["sn:Z:"] == "chr1"
    assert gaf_lines[1].query_name == "read_s464827_revcomp"
    assert gaf_lines[1].tags["bo:i:"] == "3"
    assert gaf_lines[1].tags["sn:Z:"] == "unknown"
    assert gaf_lines[2].query_name == "read_s4_s5_100_100"
    assert gaf_lines[2].tags["bo:i:"] == "4"
    assert gaf_lines[2].tags["sn:Z:"] == "chr1"
    assert gaf_lines[3].query_name == "read_s4_s6_100_100"
    assert gaf_lines[3].tags["bo:i:"] == "4"
    assert gaf_lines[3].tags["sn:Z:"] == "chr1"
    assert gaf_lines[4].query_name == "read_s4_s5_95_105_revcomp"
    assert gaf_lines[4].tags["bo:i:"] == "4"
    assert gaf_lines[4].tags["sn:Z:"] == "chr1"
    assert gaf_lines[5].query_name == "read_s4_s5_90_110"
    assert gaf_lines[5].tags["bo:i:"] == "4"
    assert gaf_lines[5].tags["sn:Z:"] == "chr1"
    assert gaf_lines[6].query_name == "read_s7"
    assert gaf_lines[6].tags["bo:i:"] == "7"
    assert gaf_lines[6].tags["sn:Z:"] == "chr1"
    assert gaf_lines[7].query_name == "read_s8_s9"
    assert gaf_lines[7].tags["bo:i:"] == "8"
    assert gaf_lines[7].tags["sn:Z:"] == "chr1"


@pytest.mark.parametrize(
    "gaf_file",
    [
        "graphaligner.gaf",
        "graphaligner-stable.gaf",
        "graphaligner.gaf.gz",
        "graphaligner-stable.gaf.gz",
    ],
)
def test_sort_customgraph(tmp_path, gaf_file):
    input_gaf = f"tests/data/index_and_view/{gaf_file}"
    input_gfa = "tests/data/gfa2rgfa/reference-graph.gfa"
    ordered_gfa_dir = str(tmp_path) + "/ordered_gfa/"
    run_order_gfa(gfa_filename=input_gfa, outdir=ordered_gfa_dir)
    ordered_gfa = ordered_gfa_dir + "reference-graph-complete.gfa"
    output = str(tmp_path) + "/output.gaf"
    if "stable" in gaf_file:
        try:
            run_sort(gfa=ordered_gfa, gaf=input_gaf, outgaf=output)
        except CommandLineError:
            # expected
            return
    run_sort(gfa=ordered_gfa, gaf=input_gaf, outgaf=output)
    gaf = GAF(output)
    gaf_lines = list(gaf.read_file())
    sorted_query_names = [
        "read_s1_s2_s3_s4",
        "read_s1_s14_s4",
        "read_s1_s14_s3",
        "read_s1_s14_s3_reversed",
        "read_s3_s4_s28_s6",
        "read_s15_s16_s6_s22",
        "read_s15_s16_s6_s22_reversed",
        "read_s28_s6_s27",
        "read_s27_s32_s8",
        "read_s27_s32_s8_reversed",
        "read_s27_s32_s8",
        "read_s27_s32_s8_reversed",
        "read_s12_s13",
        "read_s18",
        "read_s18_s12",
        "read_s20",
        "read_s20_s21_s18",
        "read_s24",
        "read_s29",
    ]
    sorted_sn_tags = [
        "REF#0#CONTIG1",
        "REF#0#CONTIG1",
        "REF#0#CONTIG1",
        "REF#0#CONTIG1",
        "REF#0#CONTIG1",
        "REF#0#CONTIG1",
        "REF#0#CONTIG1",
        "unknown",
        "unknown",
        "unknown",
        "REF#0#CONTIG1",
        "REF#0#CONTIG1",
        "REF#0#CONTIG1",
        "unknown",
        "REF#0#CONTIG1",
        "unknown",
        "unknown",
        "unknown",
        "unknown",
    ]
    sorted_bo_tags = [
        "0",
        "0",
        "0",
        "0",
        "1",
        "3",
        "3",
        "5",
        "5",
        "5",
        "5",
        "5",
        "7",
        "7",
        "7",
        "7",
        "7",
        "7",
        "7",
    ]
    assert len(gaf_lines) == len(sorted_query_names)
    for i in range(len(gaf_lines)):
        assert gaf_lines[i].query_name == sorted_query_names[i]
        assert gaf_lines[i].tags["sn:Z:"] == sorted_sn_tags[i]
        assert gaf_lines[i].tags["bo:i:"] == sorted_bo_tags[i]
