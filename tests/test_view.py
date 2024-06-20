"""
Test for gaftools view.
"""

from gaftools.cli.view import run


def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split("\t"))

    return [parse_line(l) for l in open(filename)]


######################
# Stable to Unstable #
######################


# testing whole file conversion from stable to unstable
def test_stable_to_unstable(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="unstable")
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from stable to unstable with nodes specified
def test_stable_to_unstable_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s2.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="unstable", nodes=["s2"])
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from stable to unstable with regions specified
def test_stable_to_unstable_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s2.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="unstable",
        regions=["chr1:3300-3400"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from stable to unstable with non-reference nodes specified
def test_stable_to_unstable_non_reference_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s464827.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="unstable", nodes=["s464827"])
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from stable to unstable with non-reference regions specified
def test_stable_to_unstable_non_reference_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s464827.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="unstable",
        regions=["NA20129#1#JAHEPE010000248.1:4927-5000"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from stable to unstable with multiple nodes specified
def test_stable_to_unstable_multiple_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s2.gaf"
    run(
        gaf_path=input_gaf, gfa=input_gfa, output=output, format="unstable", nodes=["s464827", "s2"]
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from stable to unstable with multiple regions specified
def test_stable_to_unstable_multiple_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s2.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="unstable",
        regions=["NA20129#1#JAHEPE010000248.1:4927-5000", "chr1:3300-3400"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing whole gzipped file conversion from stable to unstable
def test_gzipped_stable_to_unstable(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="unstable")
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from stable to unstable with nodes specified
def test_gzipped_stable_to_unstable_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s2.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="unstable", nodes=["s2"])
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from stable to unstable with regions specified
def test_gzipped_stable_to_unstable_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s2.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="unstable",
        regions=["chr1:3300-3400"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from stable to unstable with non-reference nodes specified
def test_gzipped_stable_to_unstable_non_reference_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s464827.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="unstable", nodes=["s464827"])
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from stable to unstable with non-reference regions specified
def test_gzipped_stable_to_unstable_non_reference_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s464827.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="unstable",
        regions=["NA20129#1#JAHEPE010000248.1:4927-5000"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from stable to unstable with multiple nodes specified
def test_gzipped_stable_to_unstable_multiple_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s2.gaf"
    run(
        gaf_path=input_gaf, gfa=input_gfa, output=output, format="unstable", nodes=["s464827", "s2"]
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from stable to unstable with multiple regions specified
def test_gzipped_stable_to_unstable_multiple_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-stable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-unstable-conversioncheck-only-s2.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="unstable",
        regions=["NA20129#1#JAHEPE010000248.1:4927-5000", "chr1:3300-3400"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


######################
# Unstable to Stable #
######################


# testing whole file conversion from unstable to stable
def test_unstable_to_stable(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="stable")
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from unstable to stable with nodes specified
def test_unstable_to_stable_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s2.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="stable", nodes=["s2"])
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from unstable to stable with regions specified
def test_unstable_to_stable_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s2.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="stable",
        regions=["chr1:3300-3400"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from unstable to stable with non-reference nodes specified
def test_unstable_to_stable_non_reference_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s464827.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="stable", nodes=["s464827"])
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from unstable to stable with non-reference regions specified
def test_unstable_to_stable_non_reference_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s464827.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="stable",
        regions=["NA20129#1#JAHEPE010000248.1:4927-5000"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from unstable to stable with multiple nodes specified
def test_unstable_to_stable_multiple_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s2.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="stable", nodes=["s464827", "s2"])
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing file conversion from unstable to stable with multiple regions specified
def test_unstable_to_stable_multiple_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s2.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="stable",
        regions=["NA20129#1#JAHEPE010000248.1:4927-5000", "chr1:3300-3400"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing whole gzipped file conversion from unstable to stable
def test_gzipped_unstable_to_stable(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="stable")
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from unstable to stable with nodes specified
def test_gzipped_unstable_to_stable_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s2.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="stable", nodes=["s2"])
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from unstable to stable with regions specified
def test_gzipped_unstable_to_stable_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s2.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="stable",
        regions=["chr1:3300-3400"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from unstable to stable with non-reference nodes specified
def test_gzipped_unstable_to_stable_non_reference_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s464827.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="stable", nodes=["s464827"])
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from unstable to stable with non-reference regions specified
def test_gzipped_unstable_to_stable_non_reference_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s464827.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="stable",
        regions=["NA20129#1#JAHEPE010000248.1:4927-5000"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from unstable to stable with multiple nodes specified
def test_gzipped_unstable_to_stable_multiple_nodes(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s2.gaf"
    run(gaf_path=input_gaf, gfa=input_gfa, output=output, format="stable", nodes=["s464827", "s2"])
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# testing gzipped file conversion from unstable to stable with multiple regions specified
def test_gzipped_unstable_to_stable_multiple_regions(tmp_path):
    input_gaf = "tests/data/alignments-minigraph-unstable-conversioncheck.gaf.gz"
    input_gfa = "tests/data/smallgraph.gfa"
    output = str(tmp_path) + "/output.gaf"
    truth = "tests/data/alignments-minigraph-stable-conversioncheck-only-s2.gaf"
    run(
        gaf_path=input_gaf,
        gfa=input_gfa,
        output=output,
        format="stable",
        regions=["NA20129#1#JAHEPE010000248.1:4927-5000", "chr1:3300-3400"],
    )
    output_lines = parse_output(output)
    truth_lines = parse_output(truth)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
