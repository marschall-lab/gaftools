"""
Statistics of a GAF File
"""

import sys
import logging
import itertools
from gaftools.gaf import GAF, Read
from gaftools.cli import log_memory_usage
from gaftools.timer import StageTimer

logger = logging.getLogger(__name__)


def run_stat(
    gaf_path,
    cigar_stat=False,
    output=None,
):
    """This function outputs some statistics of a GAF file. If you run with "--cigar" option, then
    cigar related statistics are also output. This increases the running time more than 10 folds
    because I need to iterate through the cigar of each alignment making it run in O(n^2)
    """

    timers = StageTimer()

    if output is None:
        output = sys.stdout
    else:
        output = open(output, "w")

    # read_names = set()
    total_aligned_bases = 0
    total_mapq = 0
    total_primary = 0
    total_secondary = 0

    if cigar_stat:
        total_del = 0
        total_ins = 0
        total_x = 0
        total_del_large = 0
        total_ins_large = 0
        total_x_large = 0
        total_match = 0
        total_match_large = 0
        total_perfect = 0

    reads = {}
    gaf_file = GAF(gaf_path)
    for alignment_count, mapping in enumerate(gaf_file.read_file(), 1):
        # hashed_readname = hash(mapping.query_name)
        # read_names.add(hashed_readname)

        if not (mapping.is_primary) or (mapping.mapping_quality <= 0):
            total_secondary += 1
            continue

        map_ratio = float(mapping.query_end - mapping.query_start) / (mapping.query_length)
        seq_identity = float(mapping.residue_matches) / mapping.alignment_block_length

        if mapping.query_name not in reads:
            reads[mapping.query_name] = Read(
                mapping.query_name, mapping.query_length, map_ratio, seq_identity
            )
        else:
            reads[mapping.query_name].aln_count += 1
            # reads[mapping.query_name].total_map_ratio += map_ratio
            # reads[mapping.query_name].total_seq_identity += seq_identity

            if reads[mapping.query_name].highest_map_ratio < map_ratio:
                reads[mapping.query_name].highest_map_ratio = map_ratio

            if reads[mapping.query_name].highest_seq_identity < seq_identity:
                reads[mapping.query_name].highest_seq_identity = seq_identity

        total_aligned_bases += mapping.residue_matches
        total_mapq += mapping.mapping_quality
        total_primary += 1

        if cigar_stat:
            """Cigar string analysis"""
            cigar = mapping.cigar
            all_cigars = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
            if len(all_cigars) == 2:
                total_perfect += 1
            for cnt in range(0, len(all_cigars) - 1, 2):
                if all_cigars[cnt + 1] == "D":
                    total_del += 1
                    if int(all_cigars[cnt]) >= 50:
                        total_del_large += 1
                elif all_cigars[cnt + 1] == "I":
                    total_ins += 1
                    if int(all_cigars[cnt]) >= 50:
                        total_ins_large += 1
                elif all_cigars[cnt + 1] == "X":
                    total_x += 1
                    if int(all_cigars[cnt]) >= 50:
                        total_x_large += 1
                elif all_cigars[cnt + 1] == "=":
                    total_match += 1
                    if int(all_cigars[cnt]) >= 50:
                        total_match_large += 1
    gaf_file.close()

    # avg_total_seq_identity = 0.0
    # avg_total_map_ratio = 0.0
    avg_highest_seq_identity = 0.0
    avg_highest_map_ratio = 0.0

    for k, v in reads.items():
        # avg_total_seq_identity += float(v.total_seq_identity) / v.aln_count
        # avg_total_map_ratio += float(v.total_map_ratio) / v.aln_count
        avg_highest_seq_identity += float(v.highest_seq_identity)
        avg_highest_map_ratio += float(v.highest_map_ratio)

    # avg_total_seq_identity /= len(reads)
    # avg_total_map_ratio /= len(reads)
    avg_highest_seq_identity /= len(reads)
    avg_highest_map_ratio /= len(reads)
    print()
    print("Total alignments:", alignment_count, file=output)
    print("\tPrimary:", total_primary, file=output)
    print("\tSecondary:", total_secondary, file=output)
    print("Reads with at least one alignment:", len(reads), file=output)
    print("Total aligned bases:", str(total_aligned_bases), file=output)
    print("Average mapping quality:", round((total_mapq / alignment_count), 1), file=output)
    # print("Average total sequence identity:", round(avg_total_seq_identity, 2), file=output)
    print("Average highest sequence identity:", round(avg_highest_seq_identity, 3), file=output)
    # print("Average total map ratio:", round(avg_total_map_ratio,2), file=output)
    print("Average highest map ratio:", round(avg_highest_map_ratio, 3), file=output)

    if cigar_stat:
        print(
            "Cigar string statistics:\n\tTotal deletion regions: %d (%d >50bps)\n\tTotal insertion regions: %d (%d >50bps)\n\tTotal substitution regions: %d (%d >50bps)\n\tTotal match regions: %d (%d >50bps)"
            % (
                total_del,
                total_del_large,
                total_ins,
                total_ins_large,
                total_x,
                total_x_large,
                total_match,
                total_match_large,
            ),
            file=output,
        )

        print("Total perfect alignments (exact match):", total_perfect, file=output)

    print()
    print(
        "* Numbers are based on primary alignments and the ones with >0 mapping quality",
        file=output,
    )
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg("gaf_path", metavar="GAF", help="Input GAF file (can be bgzip-compressed)")
    arg("-o", "--output", default=None, help="Output file. If omitted, use standard output.")
    arg(
        "--cigar",
        dest="cigar_stat",
        default=False,
        action="store_true",
        help="Outputs cigar related statistics (requires more time)",
    )


def validate(args, parser):
    return True


def main(args):
    run_stat(**vars(args))
