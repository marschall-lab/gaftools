"""
Test the GFA to rGFA conversion script.
"""

import pysam.libcbgzf
from gaftools.cli.gfa2rgfa import run
import pysam
import pytest

### Parsing functions


# parse the whole file as a sorted list of tuples
def parse_output(filename, gzipped=False):
    def parse_line(line):
        return tuple(s.strip() for s in sorted(line.split("\t")))

    if gzipped:
        return [parse_line(l.decode("utf-8")) for l in pysam.libcbgzf.BGZFile(filename, "rb")]
    else:
        return [parse_line(l) for l in open(filename)]


# parse the whole file as a gfa (tags are separated)
def parse_gfa(filename, gzipped=False):
    def parse_s_line(line):
        fields = line.strip().split("\t")
        id = fields[1]
        seq = fields[2]
        tags = {}
        for field in fields[3:]:
            tag, ftype, value = field.split(":")
            tags[tag] = (ftype, value)
        sorted_tags = list(tuple((k, tags[k])) for k in sorted(tags.keys()))
        return tuple((id, seq, *sorted_tags))

    lines = []
    if gzipped:
        for l in pysam.libcbgzf.BGZFile(filename, "rb"):
            l = l.decode("utf-8")
            if l.startswith("S"):
                lines.append(parse_s_line(l))
                continue
            lines.append(tuple(s.strip() for s in sorted(l.split("\t"))))
    else:
        for l in open(filename):
            if l.startswith("S"):
                lines.append(parse_s_line(l))
                continue
            lines.append(tuple(s.strip() for s in sorted(l.split("\t"))))

    return lines


# parses the coordinates from rGFA and checks against sequence of the coordinate in assembly
def parse_coordinates(rgfa, seqfile, gzipped=False):
    def parse_seqfile(seqfile):
        assemblies = {}
        for l in open(seqfile):
            assembly, path = l.strip().split("\t")
            if "." in assembly:
                assembly = tuple(assembly.split("."))
            else:
                # for GRCh38 and CHM13
                assembly = (assembly, "0")
            assemblies[assembly] = pysam.FastaFile(path)
        return assemblies

    assemblies = parse_seqfile(seqfile)
    if gzipped:
        rgfa = pysam.libcbgzf.BGZFile(rgfa, "rb")
    else:
        rgfa = open(rgfa)
    for l in rgfa:
        if gzipped:
            l = l.decode("utf-8")
        if l.startswith("S"):
            fields = l.strip().split("\t")
            seq = fields[2]
            tags = {}
            for field in fields[3:]:
                tag, ftype, value = field.split(":")
                tags[tag] = (ftype, value)
            assert "SN" in tags
            assert "SO" in tags
            assert "SR" in tags
            assembly = tuple(tags["SN"][1].split("#"))
            if assembly[0] == "unknown":
                # this node was not identified to come from an assembly.
                # cannot perform this test.
                continue
            if len(assembly) == 2:
                assm_name = (assembly[0], "0")
                contig_name = assembly[1]
            elif len(assembly) == 3:
                assm_name = (assembly[0], assembly[1])
                contig_name = "#".join(assembly)
            start = int(tags["SO"][1])
            length = len(seq)
            end = start + length
            assert assemblies[assm_name].fetch(contig_name, start, end) == seq
    rgfa.close()


### Test functions


# testing GFA conversion to rGFA.
# output is uncompressed.
@pytest.mark.parametrize(
    "gfa_file",
    [
        "customgraph.gfa",
        "customgraph.gfa.gz",
        "gfa2rgfa/graph-invertednodes.gfa",
        "gfa2rgfa/graph-partial-tagged.gfa",
        "gfa2rgfa/graph-partial-tagged-alternateSN.gfa",
    ],
)
def test_standard_conversion_uncompressed(tmp_path, gfa_file):
    input_gfa = f"tests/data/{gfa_file}"
    truth_rgfa = "tests/data/gfa2rgfa/reference-graph.gfa"
    seqfile = "tests/data/gfa2rgfa/samples.seqfile"  # for parsing coordinates
    output = str(tmp_path) + "/output-rgfa.gfa"

    # providing REF as reference name explicitly
    run(gfa=input_gfa, reference_name="REF", output=output)
    output_lines = parse_gfa(output, gzipped=True)
    truth_lines = parse_gfa(truth_rgfa)
    # checking for tag exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    output_lines = parse_output(output, gzipped=True)
    truth_lines = parse_output(truth_rgfa)
    # checking for exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    # checking for coordinates
    parse_coordinates(output, seqfile, gzipped=True)

    # not providing REF as reference name explicitly. has to be detected from header.
    run(gfa=input_gfa, output=output)
    output_lines = parse_gfa(output, gzipped=True)
    truth_lines = parse_gfa(truth_rgfa)
    # checking for tag exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    output_lines = parse_output(output, gzipped=True)
    truth_lines = parse_output(truth_rgfa)
    # checking for exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    # checking for coordinates
    parse_coordinates(output, seqfile, gzipped=True)


# testing GFA conversion to rGFA.
# output is compressed.
@pytest.mark.parametrize(
    "gfa_file",
    [
        "customgraph.gfa",
        "customgraph.gfa.gz",
        "gfa2rgfa/graph-invertednodes.gfa",
        "gfa2rgfa/graph-partial-tagged.gfa",
        "gfa2rgfa/graph-partial-tagged-alternateSN.gfa",
    ],
)
def test_standard_conversion_compressed(tmp_path, gfa_file):
    input_gfa = f"tests/data/{gfa_file}"
    truth_rgfa = "tests/data/gfa2rgfa/reference-graph.gfa"
    seqfile = "tests/data/gfa2rgfa/samples.seqfile"  # for parsing coordinates
    output = str(tmp_path) + "/output-rgfa.gfa.gz"

    # providing REF as reference name explicitly
    run(gfa=input_gfa, reference_name="REF", output=output)
    output_lines = parse_gfa(output, gzipped=True)
    truth_lines = parse_gfa(truth_rgfa)
    # checking for tag exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    output_lines = parse_output(output, gzipped=True)
    truth_lines = parse_output(truth_rgfa)
    # checking for exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    # checking for coordinates
    parse_coordinates(output, seqfile, gzipped=True)

    # not providing REF as reference name explicitly. has to be detected from header.
    run(gfa=input_gfa, output=output)
    output_lines = parse_gfa(output, gzipped=True)
    truth_lines = parse_gfa(truth_rgfa)
    # checking for tag exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    output_lines = parse_output(output, gzipped=True)
    truth_lines = parse_output(truth_rgfa)
    # checking for exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    # checking for coordinates
    parse_coordinates(output, seqfile, gzipped=True)


# testing GFA conversion to rGFA using full seqfile
@pytest.mark.parametrize(
    "gfa_file",
    [
        "customgraph.gfa",
        "gfa2rgfa/graph-partial-tagged.gfa",
        "gfa2rgfa/graph-wrong-Wline-order.gfa",
    ],
)
def test_full_seqfile_conversion(tmp_path, gfa_file):
    input_gfa = f"tests/data/{gfa_file}"
    truth_rgfa = "tests/data/gfa2rgfa/reference-graph-seqfile.gfa"
    seqfile = "tests/data/gfa2rgfa/samples.seqfile"
    output = str(tmp_path) + "/output-rgfa.gfa"

    run(gfa=input_gfa, seqfile=seqfile, output=output)
    output_lines = parse_gfa(output, gzipped=True)
    truth_lines = parse_gfa(truth_rgfa)
    # checking for tag exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    output_lines = parse_output(output, gzipped=True)
    truth_lines = parse_output(truth_rgfa)
    # checking for exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    # checking for coordinates
    parse_coordinates(output, seqfile, gzipped=True)


# testing GFA conversion to rGFA using partial seqfile
# Note: cannot use graph-invertednodes.gfa since s24 is not in partial set of assemblies. So it is never flipped back
@pytest.mark.parametrize(
    "gfa_file",
    [
        "customgraph.gfa",
        "gfa2rgfa/graph-partial-tagged.gfa",
        "gfa2rgfa/graph-wrong-Wline-order.gfa",
    ],
)
def test_partial_seqfile_conversion(tmp_path, gfa_file):
    input_gfa = f"tests/data/{gfa_file}"
    truth_rgfa = "tests/data/gfa2rgfa/reference-graph-partial-seqfile.gfa"
    partial_seqfile = "tests/data/gfa2rgfa/samples-partial.seqfile"
    seqfile = "tests/data/gfa2rgfa/samples.seqfile"
    output = str(tmp_path) + "/output-rgfa.gfa"

    run(gfa=input_gfa, seqfile=partial_seqfile, output=output)
    output_lines = parse_gfa(output, gzipped=True)
    truth_lines = parse_gfa(truth_rgfa)
    # checking for tag exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    output_lines = parse_output(output, gzipped=True)
    truth_lines = parse_output(truth_rgfa)
    # checking for exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    # checking for coordinates
    parse_coordinates(output, seqfile, gzipped=True)


# testing GFA conversion to rGFA with reference override to FOO#2
@pytest.mark.parametrize(
    "gfa_file",
    [
        "customgraph.gfa",
        "customgraph.gfa.gz",
        "gfa2rgfa/graph-invertednodes.gfa",
        "gfa2rgfa/graph-partial-tagged.gfa",
        "gfa2rgfa/graph-partial-tagged-alternateSN.gfa",
    ],
)
def test_conversion_reference_override(tmp_path, gfa_file):
    input_gfa = f"tests/data/{gfa_file}"
    truth_rgfa = "tests/data/gfa2rgfa/reference-graph-FOO#2.gfa"
    seqfile = "tests/data/gfa2rgfa/samples.seqfile"  # for parsing coordinates
    output = str(tmp_path) + "/output-rgfa.gfa"

    # providing REF as reference name explicitly
    run(gfa=input_gfa, reference_name="FOO#2", override_reference=True, output=output)
    output_lines = parse_gfa(output, gzipped=True)
    truth_lines = parse_gfa(truth_rgfa)
    # checking for tag exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    output_lines = parse_output(output, gzipped=True)
    truth_lines = parse_output(truth_rgfa)
    # checking for exactness
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
    # checking for coordinates
    parse_coordinates(output, seqfile, gzipped=True)
