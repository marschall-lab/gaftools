"""
Realign GAF file using wavefront alignment algorithm (WFA)
"""

import logging
import sys
import pysam
import queue
import multiprocessing as mp
from dataclasses import dataclass, field
from gaftools.cli import log_memory_usage
from gaftools.timer import StageTimer
from gaftools.gaf import GAF
from gaftools.gfa import GFA
from pywfa.align import WavefrontAligner


logger = logging.getLogger(__name__)


@dataclass(order=True)
class PriorityAlignment:
    """
    Simple data class to store the sequences and their priorities so the GAF output matches the GAf input
    """

    priority: int
    seq: str
    field(compare=False)


def peak(prior_queue):
    """
    Queue and PriorityQueue do not have a peak function
    Note: for PriorityQueue, only the item at 0 is the smallest value
    that doesn't guarantee that item 1 and 2 are correctly sorted by priority
    """
    if not prior_queue.empty():
        return prior_queue.queue[0]
    else:
        return None


def all_are_alive(processes):
    """
    Check if all processes are alive
    """
    for p in processes:
        if not p.is_alive():
            return False
    return True


def one_is_alive(processes):
    """
    Check if at least one process is alive
    """
    for p in processes:
        if p.is_alive():
            return True
    return False


def all_exited(processes):
    for p in processes:
        if p.exitcode != 0:
            return False
    return True


def run_realign(gaf, graph, fasta, output=None, cores=1):
    timers = StageTimer()

    if output is None:
        output = sys.stdout
    else:
        output = open(output, "w")

    if cores > mp.cpu_count():
        logger.warning("Number of cores requested is greater than the total number of CPU cores.")
        cores = min(mp.cpu_count() - 1, cores)

    # realign_gaf(gaf, graph, fasta, output, ext, cores)
    realign_gaf(gaf, graph, fasta, output, cores)
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage(include_children=True)
    logger.info("Total time:                                  %9.2f s", total_time)


def wfa_alignment(seq_batch, qu):
    """
    Realigns the sequence using the wavefront alignment algorithm.
    """
    # if the sequence is too large, the running time and memory requirement increases too much
    # Thus, we write the original alignment back to the GAF instead of realignment
    for gaf_line, ref, query, prior_counter in seq_batch:
        if gaf_line.query_end - gaf_line.query_start > 60_000:
            out_string = (
                f"{gaf_line.query_name}\t{gaf_line.query_length}\t{gaf_line.query_start}\t"
                f"{gaf_line.query_end}\t{gaf_line.strand}\t{gaf_line.path}\t{gaf_line.path_length}\t"
                f"{gaf_line.path_start}\t{gaf_line.path_end}\t{gaf_line.residue_matches}"
                f"\t{gaf_line.alignment_block_length}\t{gaf_line.mapping_quality}"
            )
            for k in gaf_line.tags.keys():
                out_string += f"\t{k}{gaf_line.tags[k]}"
            qu.put(PriorityAlignment(prior_counter, out_string + "\n"))
            # queue.put(out_string + "\n")

        else:
            aligner = WavefrontAligner(ref)
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

            out_string = (
                f"{gaf_line.query_name}\t{gaf_line.query_length}\t{gaf_line.query_start}\t"
                f"{gaf_line.query_end}\t{gaf_line.strand}\t{gaf_line.path}\t{gaf_line.path_length}\t"
                f"{gaf_line.path_start}\t{gaf_line.path_end}\t{match}"
                f"\t{cigar_len}\t{gaf_line.mapping_quality}"
            )
            cigar = aligner.cigarstring.replace("M", "=")
            gaf_line.tags["cg:Z:"] = cigar
            # queue.put(out_string + "\n")
            for k in gaf_line.tags.keys():
                # output.write("\t%s%s"%(k, gaf_line.tags[k]))
                out_string += f"\t{k}{gaf_line.tags[k]}"
                # out_string += "\t%s%s" % (k, gaf_line.tags[k])
            qu.put(PriorityAlignment(prior_counter, out_string + "\n"))

    qu.put(None)  # sentinel for finished process


def realign_gaf(gaf, graph, fasta, output, cores=1):
    """
    Uses pyWFA (https://github.com/kcleal/pywfa)
    parallelized on the level of alignment
    """
    # tracemalloc.start()

    fastafile = pysam.FastaFile(fasta)
    step_timer = StageTimer()
    step_timer.start("read_gfa")
    graph_obj = GFA(graph)
    step_timer.stop("read_gfa")
    logger.info(f"\n Finished loading {graph}. Took: {step_timer.elapsed('read_gfa')}")
    processes = []
    align_queue = mp.Queue()

    seq_batch = []
    batch_size = 1000
    gaf_file = GAF(gaf)
    priority_counter = 0
    for line in gaf_file.read_file():
        path_sequence = graph_obj.extract_path(line.path)
        ref = path_sequence[line.path_start : line.path_end]
        query = fastafile.fetch(line.query_name, line.query_start, line.query_end)

        seq_batch.append((line, ref, query, priority_counter))
        priority_counter += 1
        if len(seq_batch) != batch_size:
            continue
        else:
            # logger.info(f"{time.time()} finished preparing one batch and adding a process for it")
            processes.append(
                mp.Process(
                    target=wfa_alignment,
                    args=(
                        seq_batch,
                        align_queue,
                    ),
                )
            )
            seq_batch = []

        if len(processes) == cores:
            p_queue = queue.PriorityQueue()
            # logger.info(f"{time.time()} Running the prepared batches of processes")
            for p in processes:
                p.start()
            n_sentinels = 0
            while n_sentinels != len(processes):  # exits when all processes finished
                try:
                    out_string_obj = align_queue.get(timeout=0.5)
                except queue.Empty:  # queue throws Empty exception after timeout
                    # check if all threads are still alive
                    if one_is_alive(processes):
                        continue
                    else:
                        # all_exited returns false if one exited with non-zero code
                        if not all_exited(processes):
                            logger.error(
                                "One of the processes had a none-zero exit code. One reason could be that one of the processes consumed too much memory and was killed"
                            )
                            sys.exit(1)
                if out_string_obj is None:  # sentinel counter to count finished processes
                    n_sentinels += 1
                else:  # priority queue to keep the output order same as input order
                    p_queue.put(out_string_obj)
                    # output.write(out_string_obj)

            for p in processes:
                p.join()
            queue_len = len(p_queue.queue)
            for _ in range(queue_len):
                output.write(p_queue.get().seq)
            processes = []
            align_queue = mp.Queue()
            p_queue = queue.PriorityQueue()
    gaf_file.close()

    if len(seq_batch) > 0:  # leftover alignments to re-align
        processes.append(
            mp.Process(
                target=wfa_alignment,
                args=(
                    seq_batch,
                    align_queue,
                ),
            )
        )
    # leftover batches
    if len(processes) != 0:
        p_queue = queue.PriorityQueue()
        for p in processes:
            p.start()
        n_sentinels = 0
        while n_sentinels != len(processes):  # exits when all processes finished
            try:
                out_string_obj = align_queue.get(timeout=0.1)
            except queue.Empty:
                # check if all threads are still alive
                if one_is_alive(processes):
                    continue
                else:
                    # all_exited returns false if one exited with non-zero code
                    if not all_exited(processes):
                        logger.error("One of the processes had a none-zero exit code")
                        sys.exit(1)
            if out_string_obj is None:
                n_sentinels += 1
            else:
                p_queue.put(out_string_obj)
                # output.write(out_string_obj)
        for p in processes:
            p.join()
        queue_len = len(p_queue.queue)
        for _ in range(queue_len):
            output.write(p_queue.get().seq)
    logger.info("Done!")

    fastafile.close()


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf', metavar='GAF', 
        help='Input GAF file (can be bgzip-compressed)')
    arg('graph', metavar='rGFA', 
        help='reference rGFA file')
    arg('fasta', metavar='FASTA', 
        help='Input FASTA file of the read')
    arg('-o', '--output', default=None, 
        help='Output GAF file. If omitted, use standard output.')
    arg("-c", "--cores", metavar="CORES", default=1, type=int, 
        help="Number of cores to use for alignments.")


def main(args):
    run_realign(**vars(args))
