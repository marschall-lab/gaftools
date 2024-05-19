"""
Realign GAF file using wavefront alignment algorithm (WFA)
"""

import logging
import sys
import pysam
import gaftools.gaf
import multiprocessing as mp
from gaftools import __version__
from gaftools.cli import log_memory_usage
from gaftools.cli import CommandLineError
from gaftools.timer import StageTimer
from gaftools.gfa import GFA
from pywfa.align import WavefrontAligner, cigartuples_to_str

logger = logging.getLogger(__name__)


def run_realign(
        gaf,
        graph,
        fasta,
        output=None,
        ext=False,
        cores=1
):
    timers = StageTimer()

    if output is None:
        output = sys.stdout
    else:
        output = open(output, 'w')

    if cores > mp.cpu_count():
        logger.warning('Number of cores requested is greater than the total number of CPU cores.')
        cores = min(mp.cpu_count() - 1, cores)

    realign_gaf(gaf, graph, fasta, output, ext, cores)
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


def overlap_ratio(x_start, x_end, y_start, y_end):
    overlap = max(0, min(x_end, y_end) - max(x_start, y_start))
    total_length = x_end - x_start + y_end - y_start
    x_length = x_end - x_start
    y_length = y_end - y_start

    return max(2 * (overlap / total_length), (overlap / x_length), (overlap / y_length))


def filter_duplicates(aln):
    import functools

    for k in aln.keys():
        aln[k].sort(key=functools.cmp_to_key(gaftools.gaf.compare_aln))

    for read_name, mappings in aln.items():
        if len(mappings) == 1:
            continue
        for cnt, line in enumerate(mappings):
            if line.duplicate == True:
                continue
            for cnt2, line2 in enumerate(mappings):
                if cnt == cnt2 or line2.duplicate == True:
                    continue

                sim = overlap_ratio(line.query_start, line.query_end, line2.query_start, line2.query_end)

                if sim > 0.75:
                    if line.score > line2.score:
                        mappings[cnt2].duplicate = True
                        # print("......Duplicate")
                    else:
                        mappings[cnt].duplicate = True
                        # print("......Duplicate, breaking")
                        break


def write_alignments(aln, output):
    # Write the alignment back to the GAF
    for read_name, mappings in aln.items():
        for gaf_line in mappings:
            if not gaf_line.duplicate:
                # Write the alignment back to the GAF
                output.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" % (gaf_line.query_name,
                                                                                 gaf_line.query_length,
                                                                                 gaf_line.query_start,
                                                                                 gaf_line.query_end, gaf_line.strand,
                                                                                 gaf_line.path, gaf_line.path_length,
                                                                                 gaf_line.path_start, gaf_line.path_end,
                                                                                 match, cigar_len,
                                                                                 gaf_line.mapping_quality))

                for k in gaf_line.tags.keys():
                    output.write("\t%s%s" % (k, gaf_line.tags[k]))

                # print(gaf_line.cigar)
                # output.write("\tcg:Z:%s\n" %gaf_line.cigar)


def wfa_alignment(aln, gaf_line, ref, query, path_start, extended, queue):
    # If the sequence is too large, the running time and memory requirement increases too much
    # Thus, we write the original alignment back to the GAF instead of realignment
    if gaf_line.query_end - gaf_line.query_start > 60000:
        # output.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" %(gaf_line.query_name,
        #                 gaf_line.query_length, gaf_line.query_start, gaf_line.query_end, gaf_line.strand,
        #                 gaf_line.path, gaf_line.path_length, gaf_line.path_start, gaf_line.path_end,
        #                 gaf_line.residue_matches, gaf_line.alignment_block_length, gaf_line.mapping_quality))
        out_string = "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" % (gaf_line.query_name,
                                                                         gaf_line.query_length, gaf_line.query_start,
                                                                         gaf_line.query_end, gaf_line.strand,
                                                                         gaf_line.path, gaf_line.path_length,
                                                                         gaf_line.path_start, gaf_line.path_end,
                                                                         gaf_line.residue_matches,
                                                                         gaf_line.alignment_block_length,
                                                                         gaf_line.mapping_quality)

        for k in gaf_line.tags.keys():
            # output.write("\t%s%s"%(k, gaf_line.tags[k]))
            out_string += "\t%s%s" % (k, gaf_line.tags[k])
        queue.put(out_string + "\n")
        # output.write("\n")
        # return

    else:
        aligner = WavefrontAligner(ref)
        if extended:
            res = aligner(query, clip_cigar=True, min_aligned_bases_left=30, min_aligned_bases_right=30)
        else:
            res = aligner(query, clip_cigar=False)

        match, mismatch, cigar_len, ins, deletion, soft_clip = 0, 0, 0, 0, 0, 0
        cigar = ""

        for op_type, op_len in res.cigartuples:
            if op_type == 0:
                match += op_len
                cigar += str(op_len) + "="
            elif op_type == 1:
                ins += op_len
                cigar += str(op_len) + "I"
            elif op_type == 2:
                deletion += op_len
                cigar += str(op_len) + "D"
            elif op_type == 4:
                soft_clip += op_len
            elif op_type == 8:
                mismatch += op_len
                cigar += str(op_len) + "X"
            else:
                assert False
            cigar_len += op_len

        if extended:
            if match < 30:
                # return
                queue.put("NA")

            if gaf_line.query_name in aln:
                # todo maybe I can add this to the queue as a tuple and when I return from the queue I check the item
                #     and if it's a string I write it, if it's a tuple I add it to aln
                #     aln seems to be a dictionary of key:query_name and value [alignment(s)]
                aln[gaf_line.query_name].append(gaftools.gaf.Alignment("", gaf_line.query_length, res.text_start,
                                                                       res.text_end, gaf_line.strand, gaf_line.path,
                                                                       gaf_line.path_length,
                                                                       path_start + res.pattern_start,
                                                                       path_start + res.pattern_end, match, 0, 0, False,
                                                                       cigar,
                                                                       cigar_len, res.score))
            else:
                aln[gaf_line.query_name] = [gaftools.gaf.Alignment("", gaf_line.query_length, res.text_start,
                                                                   res.text_end, gaf_line.strand, gaf_line.path,
                                                                   gaf_line.path_length, path_start + res.pattern_start,
                                                                   path_start + res.pattern_end, match, 0, 0, False,
                                                                   cigar, cigar_len, res.score)]
        else:
            cigar = aligner.cigarstring.replace("M", "=")

            # Write the alignment back to the GAF
            # output.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" %(gaf_line.query_name,
            #                 gaf_line.query_length, gaf_line.query_start, gaf_line.query_end, gaf_line.strand,
            #                 gaf_line.path, gaf_line.path_length, gaf_line.path_start, gaf_line.path_end,
            #                 match, cigar_len, gaf_line.mapping_quality))
            out_string = "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" % (gaf_line.query_name,
                                                                             gaf_line.query_length,
                                                                             gaf_line.query_start, gaf_line.query_end,
                                                                             gaf_line.strand,
                                                                             gaf_line.path, gaf_line.path_length,
                                                                             gaf_line.path_start, gaf_line.path_end,
                                                                             match, cigar_len, gaf_line.mapping_quality)
            # queue.put(out_string + "\n")

            for k in gaf_line.tags.keys():
                # output.write("\t%s%s"%(k, gaf_line.tags[k]))
                out_string += "\t%s%s" % (k, gaf_line.tags[k])
            queue.put(out_string + "\n")
            # output.write("\n")
    queue.put(None)


def realign_gaf(gaf, graph, fasta, output, extended, cores=1):
    """
    Uses pyWFA (https://github.com/kcleal/pywfa)
    parallelized on the level of alignment because in python cores
    """

    fastafile = pysam.FastaFile(fasta)
    graph_obj = GFA(graph)
    processes = []
    queue = mp.Queue()

    aln = {}
    for cnt, line in enumerate(gaftools.gaf.parse_gaf(gaf)):

        path_sequence = graph_obj.extract_path(line.path)

        if extended:
            extension_start = line.query_start
            extension_end = line.query_length - line.query_end

            # Also add 10% of the extension to each side
            if extension_start > 0:
                extension_start += extension_start // 10
            if extension_end > 0:
                extension_end += extension_end // 10

            path_start = line.path_start - extension_start
            path_end = line.path_end + extension_end

            if path_start < 0:
                path_start = 0
            if path_end > line.path_length:
                path_end = line.path_length

            ref = path_sequence[path_start:path_end]
            query = fastafile.fetch(line.query_name)

            # making a process
            processes.append(
                mp.Process(target=wfa_alignment, args=(aln, line, ref, query, path_start, extended, queue,)))

            # wfa_alignment(aln, line, ref, query, path_start, output, extended)
        else:
            ref = path_sequence[line.path_start:line.path_end]
            query = fastafile.fetch(line.query_name, line.query_start, line.query_end)
            processes.append(mp.Process(target=wfa_alignment, args=([], line, ref, query, 0, False, queue,)))
            # wfa_alignment([], line, ref, query, 0, output, False)
        if len(processes) == cores:
            # start running the processes and checking when they're done before starting the next batch
            for p in processes:
                p.start()
            n_sentinals = 0
            while n_sentinals != len(processes):  # exits when all processes finished
                out_string = queue.get()
                if out_string is None:
                    n_sentinals += 1
                elif out_string == "NA":
                    pass
                else:
                    output.write(out_string)

            for p in processes:
                p.join()
            processes = []
            queue = mp.Queue()
            # todo: There's a problem with storing a lot of alignments in aln and then written separately later
            #       I need to think about this more

    # leftover alignments
    if len(processes) != 0:
        for p in processes:
            p.start()
        n_sentinals = 0
        while n_sentinals != len(processes):  # exits when all processes finished
            out_string = queue.get()
            if out_string is None:
                n_sentinals += 1
            elif out_string == "NA":
                pass
            else:
                output.write(out_string)

        for p in processes:
            p.join()

    fastafile.close()

    if extended:
        filter_duplicates(aln)
        write_alignments(aln, output)


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf', metavar='GAF', help='GAF File')
    arg('graph', metavar='GFA', help='Input GFA file')
    arg('fasta', metavar='FASTA', help='Input FASTA file')
    arg('-o', '--output', default=None, help='Output file. If omitted, use standard output.')
    arg("--ext", action='store_true', help="Extend the aligned path (use for Minigraph alignments).")
    arg("-c", "--cores", metavar="CORES", default=1, type=int, help="Number of cores to use for alignments.")


def main(args):
    run_realign(**vars(args))
