"""
Test the GFA to rGFA conversion script.
"""

from gaftools.cli.gfa2rgfa import run

### Parsing functions


# parse the whole file for checking exactness
def parse_gfa(filename):
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
    for l in open(filename):
        if l.startswith("S"):
            lines.append(parse_s_line(l))
            continue
        lines.append(tuple(s.strip() for s in l.split("\t")))

    return lines


# parses the coordinates from rGFA and checks against sequence of the coordinate in assembly
def parse_coordinates():
    pass


### Test functions


def test_gfa2rgfa(tmp_path):
    input_gfa = "tests/data/graph-conversioncheck-gfa.gfa"
    seqfile = "tests/data/graph-conversioncheck-samples.seqfile"
    truth_rgfa = "tests/data/graph-conversioncheck-rgfa.gfa"
    output = str(tmp_path) + "/output-rgfa.gfa"
    run(gfa=input_gfa, reference_name="REF", reference_tagged=False, seqfile=seqfile, output=output)
    output_lines = parse_gfa(output)
    truth_lines = parse_gfa(truth_rgfa)
    for n in range(len(output_lines)):
        assert output_lines[n] == truth_lines[n]
