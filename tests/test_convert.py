"""
Test for gaftools convert.
"""

from gaftools.cli.convert import stable_to_unstable, unstable_to_stable

def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split('\t'))
    return [parse_line(l) for l in open(filename)]

def test_stable_to_unstable(tmp_path):
    input_gaf = 'tests/data/alignments-minigraph-stable-conversioncheck.gaf'
    input_gfa = 'tests/data/smallgraph.gfa'
    output = str(tmp_path) + '/output.gaf'
    truth = 'tests/data/alignments-minigraph-unstable-conversioncheck.gaf'
    stable_to_unstable(
        gaf_path=input_gaf,
        gfa_path=input_gfa,
        out_path=output
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert (output_lines[n] == truth_lines[n])

def test_unstable_to_stable(tmp_path):
    input_gaf = 'tests/data/alignments-minigraph-unstable-conversioncheck.gaf'
    input_gfa = 'tests/data/smallgraph.gfa'
    output = str(tmp_path) + '/output.gaf'
    truth = 'tests/data/alignments-minigraph-stable-conversioncheck.gaf'
    unstable_to_stable(
        gaf_path=input_gaf,
        gfa_path=input_gfa,
        out_path=output
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert (output_lines[n] == truth_lines[n])