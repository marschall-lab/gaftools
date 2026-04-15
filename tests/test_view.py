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
@pytest.mark.parametrize(
    "gfa_file", ["customgraph.gfa", "customgraph.gfa.gz", "gfa2rgfa/reference-graph.gfa"]
)
def test_conversion_to_stable(tmp_path, gaf_file, gfa_file):
    gaf = f"tests/data/index_and_view/{gaf_file}"
    gfa = f"tests/data/{gfa_file}"
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
            assert False
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
            assert False
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
@pytest.mark.parametrize(
    "gfa_file", ["customgraph.gfa", "customgraph.gfa.gz", "gfa2rgfa/reference-graph.gfa"]
)
def test_conversion_to_unstable(tmp_path, gaf_file, gfa_file):
    gaf = f"tests/data/index_and_view/{gaf_file}"
    gfa = f"tests/data/{gfa_file}"
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
            assert False
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
            assert False
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
@pytest.mark.parametrize(
    "gfa_file", ["customgraph.gfa", "customgraph.gfa.gz", "gfa2rgfa/reference-graph.gfa"]
)
def test_node_selection(tmp_path, gaf_file, gfa_file):
    view_index_file = (
        "tests/data/index_and_view/view-index-stable.gvi"
        if "stable" in gaf_file
        else "tests/data/index_and_view/view-index-unstable.gvi"
    )
    gaf = f"tests/data/index_and_view/{gaf_file}"
    gfa = f"tests/data/{gfa_file}"
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
@pytest.mark.parametrize(
    "gfa_file", ["customgraph.gfa", "customgraph.gfa.gz", "gfa2rgfa/reference-graph.gfa"]
)
def test_node_selection_withconversion(tmp_path, gaf_file, gfa_file):
    view_index_file = (
        "tests/data/index_and_view/view-index-stable.gvi"
        if "stable" in gaf_file
        else "tests/data/index_and_view/view-index-unstable.gvi"
    )
    gaf = f"tests/data/index_and_view/{gaf_file}"
    gfa = f"tests/data/{gfa_file}"
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
                assert False
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
                assert False
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
                assert False
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
                assert False
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
                assert False
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
                assert False
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
                assert False
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
                assert False
            except CommandLineError:
                # expected
                pass


# the node selection tests ran fine. Just need to check if the conversion happens correctly.
@pytest.mark.parametrize("index_file", ["view-index-stable.gvi", "view-index-unstable.gvi"])
def test_region2node(index_file):
    from gaftools.cli.view import get_unstable
    import pickle

    index = f"tests/data/index_and_view/{index_file}"
    ind = None
    with open(index, "rb") as tmp:
        ind = pickle.load(tmp)
    # Note: Absent nodes don't show up
    # testing boundaries of the regions
    selected_regions = ["REF#0#CONTIG1:60-65"]
    assert get_unstable(selected_regions, ind) == []
    selected_regions = ["REF#0#CONTIG1:59-65"]
    assert get_unstable(selected_regions, ind) == ["s4"]
    selected_regions = ["REF#0#CONTIG1:59-90"]
    assert get_unstable(selected_regions, ind) == ["s4"]
    selected_regions = ["REF#0#CONTIG1:59-91"]
    assert set(get_unstable(selected_regions, ind)) == set(["s4", "s6"])
    selected_regions = ["BAR#1#ASSM1:170-200"]
    assert get_unstable(selected_regions, ind) == []
    selected_regions = ["BAR#1#ASSM1:170-201"]
    assert get_unstable(selected_regions, ind) == ["s27"]
    selected_regions = ["BAR#1#ASSM1:40-201"]
    assert get_unstable(selected_regions, ind) == ["s27"]
    selected_regions = ["BAR#1#ASSM1:39-201"]
    assert set(get_unstable(selected_regions, ind)) == set(["s24", "s27"])
    selected_regions = ["BAR#1#ASSM1:0-10"]
    assert get_unstable(selected_regions, ind) == []
    selected_regions = ["BAR#1#ASSM1:0-11"]
    assert get_unstable(selected_regions, ind) == ["s24"]
    selected_regions = ["BAR#1#ASSM1:10-12"]
    assert get_unstable(selected_regions, ind) == ["s24"]
    selected_regions = ["REF#0#CONTIG1:59-91", "BAR#1#ASSM1:39-201"]
    assert set(get_unstable(selected_regions, ind)) == set(["s4", "s6", "s24", "s27"])
    selected_regions = ["REF#0#CONTIG1:60-95", "REF#0#CONTIG1:70-91"]
    assert get_unstable(selected_regions, ind) == ["s6"]
    selected_regions = ["REF#0#CONTIG1:5-30", "REF#0#CONTIG1:70-91"]
    assert set(get_unstable(selected_regions, ind)) == set(["s1", "s2", "s6"])
    selected_regions = ["FOO#1#ASSM1:60-85", "BAR#1#ASSM1:39-201"]
    assert set(get_unstable(selected_regions, ind)) == set(["s15", "s24", "s27"])
    selected_regions = ["FOO#1#ASSM1:60-86", "BAR#1#ASSM1:39-201"]
    assert set(get_unstable(selected_regions, ind)) == set(["s15", "s16", "s24", "s27"])


# Region conversion was testing above but making sure that these nodes are processed properly
@pytest.mark.parametrize(
    "gaf_file",
    [
        "graphaligner.gaf",
        "graphaligner.gaf.gz",
        "graphaligner-stable.gaf",
        "graphaligner-stable.gaf.gz",
    ],
)
@pytest.mark.parametrize(
    "gfa_file", ["customgraph.gfa", "customgraph.gfa.gz", "gfa2rgfa/reference-graph.gfa"]
)
def test_region_selection(tmp_path, gaf_file, gfa_file):
    view_index_file = (
        "tests/data/index_and_view/view-index-stable.gvi"
        if "stable" in gaf_file
        else "tests/data/index_and_view/view-index-unstable.gvi"
    )
    gaf = f"tests/data/index_and_view/{gaf_file}"
    gfa = f"tests/data/{gfa_file}"
    truth = (
        "tests/data/index_and_view/graphaligner-stable.gaf"
        if "stable" in gaf_file
        else "tests/data/index_and_view/graphaligner.gaf"
    )
    output = str(tmp_path) + "/output.gaf"

    ### Intersection Mode
    mode = "I"
    if True:
        # this region is not present in GAF file
        selected_regions = ["REF#0#CONTIG1:60-65"]
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            regions=selected_regions,
            mode=mode,
            output=output,
        )
        truth_line_indices = []
        truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
        output_lines = parse_output(output)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]

        # this region is not present in GAF file
        selected_regions = ["BAR#1#ASSM1:170-200"]
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            regions=selected_regions,
            mode=mode,
            output=output,
        )
        truth_line_indices = []
        truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
        output_lines = parse_output(output)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]

        selected_regions = ["REF#0#CONTIG1:0-6"]  # only s1
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            regions=selected_regions,
            mode=mode,
            output=output,
        )
        truth_line_indices = [0, 1, 13, 15]
        truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
        output_lines = parse_output(output)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]

        selected_regions = ["REF#0#CONTIG1:0-6", "FOO#1#ASSM1:11-20"]  # s1 and s14
        run(
            gaf_path=gaf,
            gfa=gfa,  # since we do not convert coordinates, the gfa does not matter.
            index=view_index_file,
            regions=selected_regions,
            mode=mode,
            output=output,
        )
        truth_line_indices = [1, 13, 15]
        truth_lines = subset_output_lines(parse_output(truth), truth_line_indices)
        output_lines = parse_output(output)
        assert len(output_lines) == len(truth_lines)
        for n in range(len(output_lines)):
            assert output_lines[n] == truth_lines[n]
