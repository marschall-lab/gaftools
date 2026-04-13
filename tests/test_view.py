"""
Test for gaftools view.
"""

from gaftools.cli.view import run
from gaftools.cli import CommandLineError
from pysam import libcbgzf
import pytest


def parse_line(line):
    return tuple(s.strip() for s in line.split("\t"))


def parse_output(filename):
    return [parse_line(l) for l in open(filename)]


def parse_string_block(text):
    return [parse_line(l) for l in text.strip().split("\n")]


def parse_bgzip_output(filename):
    with libcbgzf.BGZFile(filename, "rb") as f:
        return [parse_line(line.decode("utf-8").strip()) for line in f]


def subset_output_lines(lines, indices):
    out = []
    for i in indices:
        out.append(lines[i])
    return out


# No conversion. Full file viewing.
@pytest.mark.parametrize(
    "gaf_file",
    [
        "graphaligner.gaf",
        "graphaligner.gaf.gz",
        "graphaligner-stable.gaf",
        "graphaligner-stable.gaf.gz",
    ],
)
def test_full_file_view(tmp_path, gaf_file):
    gaf = f"tests/data/index_and_view/{gaf_file}"
    output = str(tmp_path) + "/output.gaf"
    run(gaf_path=gaf, output=output)
    output_lines = parse_output(output)
    if gaf.endswith(".gz"):
        truth_lines = parse_bgzip_output(gaf)
    else:
        truth_lines = parse_output(gaf)
    assert len(output_lines) == len(truth_lines)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


# No conversion. Full file viewing in bgzip
@pytest.mark.parametrize(
    "gaf_file",
    [
        "graphaligner.gaf",
        "graphaligner.gaf.gz",
        "graphaligner-stable.gaf",
        "graphaligner-stable.gaf.gz",
    ],
)
def test_full_file_view_bgzip(tmp_path, gaf_file):
    gaf = f"tests/data/index_and_view/{gaf_file}"
    output = str(tmp_path) + "/output.gaf.gz"
    run(gaf_path=gaf, output=output)
    output_lines = parse_bgzip_output(output)
    if gaf.endswith(".gz"):
        truth_lines = parse_bgzip_output(gaf)
    else:
        truth_lines = parse_output(gaf)
    assert len(output_lines) == len(truth_lines)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]


@pytest.mark.parametrize(
    "gaf_file",
    [
        "graphaligner.gaf",
        "graphaligner.gaf.gz",
        # "graphaligner-stable.gaf",
        # "graphaligner-stable.gaf.gz",
    ],
)
@pytest.mark.parametrize("gfa_file", ["graph.gfa", "graph.gfa.gz", "reference-graph.gfa"])
def test_conversion_to_stable(tmp_path, gaf_file, gfa_file):
    gaf = f"tests/data/index_and_view/{gaf_file}"
    gfa = f"tests/data/gfa2rgfa/{gfa_file}"
    truth = "tests/data/index_and_view/graphaligner-stable.gaf"
    # non-bgzip output
    output = str(tmp_path) + "/output.gaf"
    if "reference" in gfa:
        run(gaf_path=gaf, format="stable", gfa=gfa, output=output)
        output_lines = parse_output(output)
        truth_lines = parse_output(truth)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]
    else:
        try:
            run(gaf_path=gaf, format="stable", gfa=gfa, output=output)
        except CommandLineError:
            # expected
            pass

    # bgzip output
    output = str(tmp_path) + "/output.gaf.gz"
    if "reference" in gfa:
        run(gaf_path=gaf, format="stable", gfa=gfa, output=output)
        output_lines = parse_bgzip_output(output)
        truth_lines = parse_output(truth)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]
    else:
        try:
            run(gaf_path=gaf, format="stable", gfa=gfa, output=output)
        except CommandLineError:
            # expected
            pass


@pytest.mark.parametrize(
    "gaf_file",
    [
        # "graphaligner.gaf",
        # "graphaligner.gaf.gz",
        "graphaligner-stable.gaf",
        "graphaligner-stable.gaf.gz",
    ],
)
@pytest.mark.parametrize("gfa_file", ["graph.gfa", "graph.gfa.gz", "reference-graph.gfa"])
def test_conversion_to_unstable(tmp_path, gaf_file, gfa_file):
    gaf = f"tests/data/index_and_view/{gaf_file}"
    gfa = f"tests/data/gfa2rgfa/{gfa_file}"
    truth = "tests/data/index_and_view/graphaligner.gaf"
    # non-bgzip output
    output = str(tmp_path) + "/output.gaf"
    if "reference" in gfa:
        run(gaf_path=gaf, format="unstable", gfa=gfa, output=output)
        output_lines = parse_output(output)
        truth_lines = parse_output(truth)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]
    else:
        try:
            run(gaf_path=gaf, format="unstable", gfa=gfa, output=output)
        except CommandLineError:
            # expected
            pass

    # bgzip output
    output = str(tmp_path) + "/output.gaf.gz"
    if "reference" in gfa:
        run(gaf_path=gaf, format="unstable", gfa=gfa, output=output)
        output_lines = parse_bgzip_output(output)
        truth_lines = parse_output(truth)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]
    else:
        try:
            run(gaf_path=gaf, format="unstable", gfa=gfa, output=output)
        except CommandLineError:
            # expected
            pass


# Here we only subset. No conversion
@pytest.mark.parametrize(
    "gaf_file",
    [
        "graphaligner.gaf",
        "graphaligner.gaf.gz",
        "graphaligner-stable.gaf",
        "graphaligner-stable.gaf.gz",
    ],
)
@pytest.mark.parametrize("gfa_file", ["graph.gfa", "graph.gfa.gz", "reference-graph.gfa"])
def test_node_selection(tmp_path, gaf_file, gfa_file):
    view_index_file = (
        "tests/data/index_and_view/view-index-stable.gvi"
        if "stable" in gaf_file
        else "tests/data/index_and_view/view-index-unstable.gvi"
    )
    gaf = f"tests/data/index_and_view/{gaf_file}"
    gfa = f"tests/data/gfa2rgfa/{gfa_file}"
    truth = (
        "tests/data/index_and_view/graphaligner-stable.gaf"
        if "stable" in gaf_file
        else "tests/data/index_and_view/graphaligner.gaf"
    )
    output = str(tmp_path) + "/output.gaf"

    ### Intersection Mode
    mode = "I"
    if True:
        # node not present in GAF
        selected_nodes = ["s5"]
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            nodes=selected_nodes,
            mode=mode,
            output=output,
        )
        output_lines = parse_output(output)
        assert len(output_lines) == 0

        selected_nodes = ["s1"]
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            nodes=selected_nodes,
            mode=mode,
            output=output,
        )
        truth_line_indices = [0, 1, 13, 15]
        truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
        output_lines = parse_output(output)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]

        selected_nodes = ["s1", "s3"]
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            nodes=selected_nodes,
            mode=mode,
            output=output,
        )
        truth_line_indices = [0, 13, 15]
        truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
        output_lines = parse_output(output)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]

        selected_nodes = ["s1", "s14", "s3", "s4"]
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            nodes=selected_nodes,
            mode=mode,
            output=output,
        )
        truth_line_indices = []
        truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
        output_lines = parse_output(output)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]

    ### Union Mode
    mode = "U"
    if True:
        # node not present in GAF
        selected_nodes = ["s5"]
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            nodes=selected_nodes,
            mode=mode,
            output=output,
        )
        output_lines = parse_output(output)
        assert len(output_lines) == 0

        selected_nodes = ["s1"]
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            nodes=selected_nodes,
            mode=mode,
            output=output,
        )
        truth_line_indices = [0, 1, 13, 15]
        truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
        output_lines = parse_output(output)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]

        selected_nodes = ["s1", "s3"]
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            nodes=selected_nodes,
            mode=mode,
            output=output,
        )
        truth_line_indices = [0, 1, 2, 13, 15]
        truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
        output_lines = parse_output(output)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]

        selected_nodes = ["s1", "s14", "s3", "s4"]
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            nodes=selected_nodes,
            mode=mode,
            output=output,
        )
        truth_line_indices = [0, 1, 2, 13, 15]
        truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
        output_lines = parse_output(output)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]


# Here we only subset and convert
@pytest.mark.parametrize(
    "gaf_file",
    [
        "graphaligner.gaf",
        "graphaligner.gaf.gz",
        "graphaligner-stable.gaf",
        "graphaligner-stable.gaf.gz",
    ],
)
@pytest.mark.parametrize("gfa_file", ["graph.gfa", "graph.gfa.gz", "reference-graph.gfa"])
def test_node_selection_withconversion(tmp_path, gaf_file, gfa_file):
    view_index_file = (
        "tests/data/index_and_view/view-index-stable.gvi"
        if "stable" in gaf_file
        else "tests/data/index_and_view/view-index-unstable.gvi"
    )
    gaf = f"tests/data/index_and_view/{gaf_file}"
    gfa = f"tests/data/gfa2rgfa/{gfa_file}"
    truth = (
        "tests/data/index_and_view/graphaligner.gaf"
        if "stable" in gaf_file
        else "tests/data/index_and_view/graphaligner-stable.gaf"
    )
    format = "unstable" if "stable" in gaf else "stable"
    output = str(tmp_path) + "/output.gaf"

    ### Intersection Mode
    mode = "I"
    if True:
        # node not present in GAF
        selected_nodes = ["s5"]
        if "reference" in gfa:
            run(
                gaf_path=gaf,
                gfa=gfa,
                index=view_index_file,
                nodes=selected_nodes,
                mode=mode,
                format=format,
                output=output,
            )
            output_lines = parse_output(output)
            assert len(output_lines) == 0
        else:
            try:
                run(
                    gaf_path=gaf,
                    gfa=gfa,
                    index=view_index_file,
                    nodes=selected_nodes,
                    mode=mode,
                    format=format,
                    output=output,
                )
            except CommandLineError:
                # expected
                pass

        selected_nodes = ["s1"]
        if "reference" in gfa:
            run(
                gaf_path=gaf,
                gfa=gfa,
                index=view_index_file,
                nodes=selected_nodes,
                mode=mode,
                format=format,
                output=output,
            )
            truth_line_indices = [0, 1, 13, 15]
            truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
            output_lines = parse_output(output)
            assert len(output_lines) == len(truth_lines)
            for n in range(len(output_lines)):
                assert output_lines[n] == truth_lines[n]
        else:
            try:
                run(
                    gaf_path=gaf,
                    gfa=gfa,
                    index=view_index_file,
                    nodes=selected_nodes,
                    mode=mode,
                    format=format,
                    output=output,
                )
            except CommandLineError:
                # expected
                pass

        selected_nodes = ["s1", "s3"]
        if "reference" in gfa:
            run(
                gaf_path=gaf,
                gfa=gfa,
                index=view_index_file,
                nodes=selected_nodes,
                mode=mode,
                format=format,
                output=output,
            )
            truth_line_indices = [0, 13, 15]
            truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
            output_lines = parse_output(output)
            assert len(output_lines) == len(truth_lines)
            for n in range(len(output_lines)):
                assert output_lines[n] == truth_lines[n]
        else:
            try:
                run(
                    gaf_path=gaf,
                    gfa=gfa,
                    index=view_index_file,
                    nodes=selected_nodes,
                    mode=mode,
                    format=format,
                    output=output,
                )
            except CommandLineError:
                # expected
                pass

        selected_nodes = ["s1", "s14", "s3", "s4"]
        if "reference" in gfa:
            run(
                gaf_path=gaf,
                gfa=gfa,
                index=view_index_file,
                nodes=selected_nodes,
                mode=mode,
                format=format,
                output=output,
            )
            truth_line_indices = []
            truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
            output_lines = parse_output(output)
            assert len(output_lines) == len(truth_lines)
            for n in range(len(output_lines)):
                assert output_lines[n] == truth_lines[n]
        else:
            try:
                run(
                    gaf_path=gaf,
                    gfa=gfa,
                    index=view_index_file,
                    nodes=selected_nodes,
                    mode=mode,
                    format=format,
                    output=output,
                )
            except CommandLineError:
                # expected
                pass

    ### Union Mode
    mode = "U"
    if True:
        # node not present in GAF
        selected_nodes = ["s5"]
        if "reference" in gfa:
            run(
                gaf_path=gaf,
                gfa=gfa,
                index=view_index_file,
                nodes=selected_nodes,
                mode=mode,
                format=format,
                output=output,
            )
            output_lines = parse_output(output)
            assert len(output_lines) == 0
        else:
            try:
                run(
                    gaf_path=gaf,
                    gfa=gfa,
                    index=view_index_file,
                    nodes=selected_nodes,
                    mode=mode,
                    format=format,
                    output=output,
                )
            except CommandLineError:
                # expected
                pass

        selected_nodes = ["s1"]
        if "reference" in gfa:
            run(
                gaf_path=gaf,
                gfa=gfa,
                index=view_index_file,
                nodes=selected_nodes,
                mode=mode,
                format=format,
                output=output,
            )
            truth_line_indices = [0, 1, 13, 15]
            truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
            output_lines = parse_output(output)
            assert len(output_lines) == len(truth_lines)
            for n in range(len(output_lines)):
                assert output_lines[n] == truth_lines[n]
        else:
            try:
                run(
                    gaf_path=gaf,
                    gfa=gfa,
                    index=view_index_file,
                    nodes=selected_nodes,
                    mode=mode,
                    format=format,
                    output=output,
                )
            except CommandLineError:
                # expected
                pass

        selected_nodes = ["s1", "s3"]
        if "reference" in gfa:
            run(
                gaf_path=gaf,
                gfa=gfa,
                index=view_index_file,
                nodes=selected_nodes,
                mode=mode,
                format=format,
                output=output,
            )
            truth_line_indices = [0, 1, 2, 13, 15]
            truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
            output_lines = parse_output(output)
            assert len(output_lines) == len(truth_lines)
            for n in range(len(output_lines)):
                assert output_lines[n] == truth_lines[n]
        else:
            try:
                run(
                    gaf_path=gaf,
                    gfa=gfa,
                    index=view_index_file,
                    nodes=selected_nodes,
                    mode=mode,
                    format=format,
                    output=output,
                )
            except CommandLineError:
                # expected
                pass

        selected_nodes = ["s1", "s14", "s3", "s4"]
        if "reference" in gfa:
            run(
                gaf_path=gaf,
                gfa=gfa,
                index=view_index_file,
                nodes=selected_nodes,
                mode=mode,
                format=format,
                output=output,
            )
            truth_line_indices = [0, 1, 2, 13, 15]
            truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
            output_lines = parse_output(output)
            assert len(output_lines) == len(truth_lines)
            for n in range(len(output_lines)):
                assert output_lines[n] == truth_lines[n]
        else:
            try:
                run(
                    gaf_path=gaf,
                    gfa=gfa,
                    index=view_index_file,
                    nodes=selected_nodes,
                    mode=mode,
                    format=format,
                    output=output,
                )
            except CommandLineError:
                # expected
                pass
