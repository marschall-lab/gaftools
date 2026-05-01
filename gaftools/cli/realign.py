"""
Realign a GAF file using wavefront alignment algorithm (WFA).
"""

import logging
import queue
import gzip
import tempfile
import multiprocessing as mp
from typing import Optional
from dataclasses import dataclass, field
from gaftools.cli import log_memory_usage
from gaftools.errors import (
    CommandLineError,
    ReadNotFoundError,
    IncorrectReadFormatError,
    IncorrectGafFormatError,
)
from gaftools.timer import StageTimer
from gaftools.gaf import GAF
from gaftools.gfa import GFA
from gaftools.utils import FileWriter, is_file_gzipped
import pyfastx
from pywfa.align import WavefrontAligner

logger = logging.getLogger(__name__)
timers = StageTimer()


@dataclass(order=True)
class PriorityAlignment:
    """
    Simple data class to store the sequences and their priorities so the GAF output matches the GAf input
    """

    priority: int
    seq: str
    field(compare=False)


@dataclass
class ReadAccessor:
    """
    Random-access wrapper for FASTA/FASTQ reads using pyfastx instead of pysam.
    """

    reader: object
    tempdir: Optional[tempfile.TemporaryDirectory] = None

    def fetch(self, query_name, query_start, query_end):
        # Avoid relying on `query_name in reader` semantics, which can vary by pyfastx build.
        try:
            record = self.reader[query_name]
        except (KeyError, IndexError, RuntimeError):
            raise ReadNotFoundError(
                f"Read '{query_name}' not found in input reads file. This is a known issue. "
                "gaftools currently uses pyfastx which has issues with FASTQ read fetching. "
                "Consider converting to FASTA to avoid this issue."
            )
        seq = record.seq if hasattr(record, "seq") else str(record)
        return seq[query_start:query_end]

    def close(self):
        # instead of the pysam FastaFile close
        if self.tempdir is not None:
            self.tempdir.cleanup()


def detect_reads_format(reads_path):
    # Peek at the first non-empty record line to decide which pyfastx reader to use.
    opener = gzip.open if is_file_gzipped(reads_path) else open
    with opener(reads_path, "rt") as reads_file:
        for line in reads_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                return "fasta"
            if line.startswith("@"):
                return "fastq"
            break
    raise IncorrectReadFormatError(
        f"Could not determine read format for '{reads_path}'. Expected FASTA or FASTQ."
    )


def open_reads(reads_path):
    reads_format = detect_reads_format(reads_path)
    # Keep pyfastx index files in a temp directory instead of next to the input reads.
    tempdir = tempfile.TemporaryDirectory(prefix="gaftools-realign-pyfastx-")
    index_path = f"{tempdir.name}/reads.fxi"
    if reads_format == "fasta":
        return ReadAccessor(pyfastx.Fasta(reads_path, index_file=index_path), tempdir=tempdir)
    return ReadAccessor(pyfastx.Fastq(reads_path, index_file=index_path), tempdir=tempdir)


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
    writer = FileWriter(output)
    if cores > mp.cpu_count():
        logger.warning("Number of cores requested is greater than the total number of CPU cores.")
        cores = min(mp.cpu_count() - 1, cores)
    # realign_gaf(gaf, graph, fasta, output, ext, cores)
    realign_gaf(gaf, graph, fasta, writer, cores)
    writer.close()
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage(include_children=True)
    logger.info("Time to read GFA:                            %9.2f s", timers.elapsed("read_gfa"))
    logger.info("Time to read, realign, and write GAF:        %9.2f s", timers.elapsed("realign"))
    logger.info("Time spent on rest:                          %9.2f s", total_time - timers.sum())
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


def realign_gaf(gaf, graph, fasta, writer, cores=1):
    """
    Uses pyWFA (https://github.com/kcleal/pywfa)
    parallelized on the level of alignment
    """
    # tracemalloc.start()

    logger.debug(f"Using {cores} cores for realignment")
    logger.debug(f"Reading reads file: {fasta}")
    reads = open_reads(fasta)  # this used to be pysam.FastaFile
    logger.info("Complete initializing reads file.")
    try:
        logger.debug(f"Reading GFA: {graph}")
        with timers("read_gfa"):
            graph_obj = GFA(graph)
        logger.info("Complete initializing GFA file.")
        processes = []
        align_queue = mp.Queue()

        seq_batch = []
        batch_size = 1000
        logger.debug(f"Reading GAF: {gaf}")
        gaf_file = GAF(gaf)
        logger.info("Complete initializing GAF file.")
        priority_counter = 0
        logger.info("Reading GAF file and preparing batches for alignment.")
        with timers("realign"):
            logger.debug("Starting realignment of GAF file.")
            for line in gaf_file.read_file():
                if line.detect_path_format():
                    raise IncorrectGafFormatError(
                        "Detected stable coordinates in GAF. Realignment requires unstable coordinates. "
                        "Please convert the GAF with gaftools view."
                    )
                path_sequence = graph_obj.extract_path(line.path)
                ref = path_sequence[line.path_start : line.path_end]
                # realign keeps GAF coordinates, so slice only the mapped query interval.
                query = reads.fetch(
                    line.query_name, line.query_start, line.query_end
                )  # was also fetch but from pysam library

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
                                    raise CommandLineError(
                                        "One of the processes had a none-zero exit code. One reason could be that one of the processes consumed too much memory and was killed"
                                    )
                        if out_string_obj is None:  # sentinel counter to count finished processes
                            n_sentinels += 1
                        else:  # priority queue to keep the output order same as input order
                            p_queue.put(out_string_obj)
                            # output.write(out_string_obj)

                    for p in processes:
                        p.join()
                    queue_len = len(p_queue.queue)
                    for _ in range(queue_len):
                        writer.write(p_queue.get().seq)
                    processes = []
                    align_queue = mp.Queue()
                    p_queue = queue.PriorityQueue()
            gaf_file.close()
            logger.debug("Completed realignment of GAF file.")
            logger.debug("Re-aligning leftover alignments.")
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
                                raise CommandLineError(
                                    "One of the processes had a none-zero exit code"
                                )
                    if out_string_obj is None:
                        n_sentinels += 1
                    else:
                        p_queue.put(out_string_obj)
                        # output.write(out_string_obj)
                for p in processes:
                    p.join()
                queue_len = len(p_queue.queue)
                for _ in range(queue_len):
                    writer.write(p_queue.get().seq)
            logger.debug("Finished all re-alignments along with the leftover alignments.")
        logger.info("Realignment complete.")
    finally:
        reads.close()


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf', metavar='GAF',
        help='GAF file (can be bgzip-compressed)')
    arg('graph', metavar='GFA',
        help='GFA file (can be bgzip-compressed)')
    arg('fasta', metavar='READS',
        help='Reads file in FASTA or FASTQ format (optionally gzipped)')
    arg('-o', '--output', default=None,
        help='Output GAF file (bgzipped if the file ends with .gz). If omitted, use standard output.')
    arg("-c", "--cores", metavar="CORES", default=1, type=int,
        help="Number of cores to use for alignments.")


def main(args):
    run_realign(**vars(args))
