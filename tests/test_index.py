"""
Test for gaftools index.
"""

import pickle
from gaftools.cli.index import run
import pytest


# parsing the index file which is a pickled dictionary
def parse_output(filename):
    data = None
    with open(filename, "rb") as handle:
        data = pickle.load(handle)
    return data


@pytest.mark.parametrize(
    "gaf_file",
    [
        "customgraph.gaf",
        "customgraph.gaf.gz",
        "customgraph-stable.gaf",
        "customgraph-stable.gaf.gz",
    ],
)
def test_index(tmp_path, gaf_file):
    input_gaf = f"tests/data/index_and_view/{gaf_file}"
    input_gfa = "tests/data/gfa2rgfa/reference-graph.gfa"
    output = str(tmp_path) + "/output.gvi"
    truth = (
        "tests/data/index_and_view/customgraph-stable.gvi"
        if "stable" in gaf_file
        else "tests/data/index_and_view/customgraph-unstable.gvi"
    )
    run(gaf_path=input_gaf, gfa_path=input_gfa, output=output)
    output_dict = parse_output(output)
    truth_dict = parse_output(truth)
    assert len(output_dict.keys()) == len(truth_dict.keys())
    for id in output_dict.keys():
        assert id in truth_dict.keys()
        output_offsets = sorted(output_dict[id])
        truth_offsets = sorted(truth_dict[id])
        assert len(output_offsets) == len(truth_offsets)
        for i in range(len(output_offsets)):
            assert truth_offsets[i] == output_offsets[i]
