"""
Test for gaftools index.
"""

import pickle
from gaftools.cli.index import run


# parsing the index file which is a pickled dictionary
def parse_output(filename):
    data = None
    with open(filename, "rb") as handle:
        data = pickle.load(handle)
    return data


# Tests


# testing index of gaf with stable coordinates
def test_index_stable(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gvi"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck.gaf.gvi"
    run(gaf_path=input_gaf, gfa_path=input_gfa, output=output)
    output_dict = parse_output(output)
    truth_dict = parse_output(truth)
    assert output_dict == truth_dict


# testing index of gaf with stable coordinates
def test_index_unstable(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gvi"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf.gvi"
    run(gaf_path=input_gaf, gfa_path=input_gfa, output=output)
    output_dict = parse_output(output)
    truth_dict = parse_output(truth)
    assert output_dict == truth_dict


# testing index of gzipped gaf with stable coordinates
def test_index_stable_gzipped(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gvi"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck.gaf.gz.gvi"
    run(gaf_path=input_gaf, gfa_path=input_gfa, output=output)
    output_dict = parse_output(output)
    truth_dict = parse_output(truth)
    assert output_dict == truth_dict


# testing index of gzipped gaf with stable coordinates
def test_index_unstable_gzipped(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gvi"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf.gz.gvi"
    run(gaf_path=input_gaf, gfa_path=input_gfa, output=output)
    output_dict = parse_output(output)
    truth_dict = parse_output(truth)
    assert output_dict == truth_dict
