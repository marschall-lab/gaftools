"""
Tests for 'gaftools realign'
"""

import gzip
import pytest

from gaftools.gaf import GAF
from gaftools.cli.realign import run_realign


def write_fastq_from_fasta(fasta_path, fastq_path):
    # Making a fake FASTQ from the test FASTA to test that gaftools is reading both
    # The quality here doesn't matter, so it's just a placeholder
    with open(fasta_path, "r") as fasta_file, open(fastq_path, "w") as fastq_file:
        name = None
        sequence_parts = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequence = "".join(sequence_parts)
                    fastq_file.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")
                    print(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}")
                name = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line)
        if name is not None:
            sequence = "".join(sequence_parts)
            fastq_file.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            print(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}")


# Helper to validate the same test data through gzipped input paths.
def write_gzip_file(src, dest):
    with open(src, "rb") as src_file, gzip.open(dest, "wb") as dest_file:
        dest_file.write(src_file.read())


# this tells pytest to run the test for all these file formats basically
@pytest.mark.parametrize("reads_file", ["reads.fa", "reads.fastq", "reads.fa.gz", "reads.fastq.gz"])
@pytest.mark.parametrize("gaf_file", ["graphaligner.gaf", "graphaligner.gaf.gz"])
def test_realign(tmp_path, reads_file, gaf_file):
    input_gaf = f"tests/data/realign/{gaf_file}"
    output_gaf = str(tmp_path) + "/output.gaf"
    reads_path = f"tests/data/realign/{reads_file}"

    """
    # testing different formats
    if reads_format == "fastq":
        reads_path = str(tmp_path / "reads.fastq")
        write_fastq_from_fasta("tests/data/realign/reads.fa", reads_path)
    elif reads_format == "fasta.gz":
        reads_path = str(tmp_path / "reads.fa.gz")
        write_gzip_file("tests/data/realign/reads.fa", reads_path)
    elif reads_format == "fastq.gz":
        reads_fastq = str(tmp_path / "reads.fastq")
        write_fastq_from_fasta("tests/data/realign/reads.fa", reads_fastq)
        reads_path = str(tmp_path / "reads.fastq.gz")
        write_gzip_file(reads_fastq, reads_path)
    """

    run_realign(
        gaf=input_gaf,
        graph="tests/data/smallgraph.gfa",
        fasta=reads_path,
        output=output_gaf,
        cores=1,
    )

    gaf = GAF(output_gaf)
    gaf_lines = list(gaf.read_file())
    assert len(gaf_lines) == 2

    assert gaf_lines[0].query_name == "read_s8_s9"
    assert gaf_lines[0].query_length == 398
    assert gaf_lines[0].query_start == 0
    assert gaf_lines[0].query_end == 398
    assert gaf_lines[0].strand == "+"
    assert gaf_lines[0].path == ">s8>s9"
    assert gaf_lines[0].path_length == 10109
    assert gaf_lines[0].path_start == 6166
    assert gaf_lines[0].path_end == 6564
    assert gaf_lines[0].residue_matches == 398
    assert gaf_lines[0].alignment_block_length == 398
    assert gaf_lines[0].mapping_quality == 60
    assert gaf_lines[0].cigar == "398="
    # assert gaf_lines[0].is_primary

    assert gaf_lines[1].query_name == "read_s8_s9_deletion15"
    assert gaf_lines[1].query_length == 383
    assert gaf_lines[1].query_start == 0
    assert gaf_lines[1].query_end == 383
    assert gaf_lines[1].strand == "+"
    assert gaf_lines[1].path == ">s8>s9"
    assert gaf_lines[1].path_length == 10109
    assert gaf_lines[1].path_start == 6166
    assert gaf_lines[1].path_end == 6564
    assert gaf_lines[1].residue_matches == 383
    assert gaf_lines[1].alignment_block_length == 398
    assert gaf_lines[1].mapping_quality == 60
    assert gaf_lines[1].cigar == "189=15D194="
    # assert gaf_lines[1].is_primary
