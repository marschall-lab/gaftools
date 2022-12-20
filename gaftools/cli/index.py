"""
Index the GAF File
"""

import logging
import re

from gaftools import __version__
from gaftools.timer import StageTimer
from gaftools.cli import log_memory_usage

class Node1:
    def __init__(self, node_id, start, end):
        self.node_id = node_id
        self.start = start
        self.end = end

logger = logging.getLogger(__name__)

def run(gaf_path, reference, output=None):
    
    import gzip
    import copy
    from gaftools.cli.sort import gfa_sort
    from gaftools.utils import is_file_gzipped
    import pickle
    from pysam import libcbgzf

    timers = StageTimer()
    if output == None:
        output = gaf_path+".gai"
    
    #Detecting if GAF has stable or unstable coordinate
    if is_file_gzipped(gaf_path):
        gaf_file = libcbgzf.BGZFile(gaf_path,"rb")
        line = gaf_file.readline().decode("utf-8")
    else:
        gaf_file = open(gaf_path,"rt")
        line = gaf_file.readline()
    unstable = detect_unstable(line)
    gaf_file.close()
    
    nodes = {}
    ref = {}
    ref_contig = []
    if not unstable:
        logger.info("Detected stable coordinates in the GAF file.")
        with timers("sort_gfa"):
            logger.info("Sorting GFA File")
            gfa_lines = gfa_sort(reference, None, True)
        contig_name = None
        with timers("store_contig_info"):
            logger.info("Storing Contig Informationa")
            for gfa_line in gfa_lines:
                tmp_contig_name = [k for k in gfa_line if k.startswith("SN:Z:")][0][5:]
                
                if tmp_contig_name != contig_name:
                    contig_name = copy.deepcopy(tmp_contig_name)
                    if contig_name not in ref:
                        ref[contig_name] = []

                start_pos = int([k for k in gfa_line if k.startswith("SO:i:")][0][5:])
                end_pos = int([k for k in gfa_line if k.startswith("LN:i:")][0][5:]) + start_pos
                tmp = Node1(gfa_line[1], start_pos, end_pos)
                rank = int([k for k in gfa_line if k.startswith("SR:i:")][0][5:])
                try:
                    rank = int([k for k in gfa_line if k.startswith("SR:i:")][0][5:])
                except IndexError:
                    logger.error("ERROR: No Rank present in the reference GFA File. Input rGFA file should have SR field.")
                    exit()
                ref[contig_name].append(tmp)
                nodes[gfa_line[1]] = (gfa_line[1], tmp_contig_name, start_pos, end_pos)
                if tmp_contig_name not in ref_contig:
                    if rank == 0:
                        ref_contig.append(contig_name)
                else:
                    assert (rank == 0)
    else:
        logger.info("Detected unstable coordinates in the GAF file.")
        gz_flag = reference[-2:] == "gz"
        if gz_flag:
            gfa_file = gzip.open(reference,"r")
        else:
            gfa_file = open(reference,"r")
        for gfa_line in gfa_file:
            if gz_flag:
                gfa_line = gfa_line.decode("utf-8")
            if gfa_line[0] != "S":
                break
            gfa_line = gfa_line.rstrip().split('\t')
            contig_name = [k for k in gfa_line if k.startswith("SN:Z:")][0][5:]
            start_pos = int([k for k in gfa_line if k.startswith("SO:i:")][0][5:])
            end_pos = int([k for k in gfa_line if k.startswith("LN:i:")][0][5:]) + start_pos
            rank = int([k for k in gfa_line if k.startswith("SR:i:")][0][5:])
            try:
                rank = int([k for k in gfa_line if k.startswith("SR:i:")][0][5:])
            except IndexError:
                logger.error("ERROR: No Rank present in the reference GFA File. Input rGFA file should have SR field.")
                exit()
            nodes[gfa_line[1]] = (gfa_line[1], contig_name, start_pos, end_pos)
            if contig_name not in ref_contig:
                if rank == 0:
                    ref_contig.append(contig_name)
            else:
                assert (rank == 0)
    
    if is_file_gzipped(gaf_path):
        logger.info("GAF file compression detected. BGZF compression needed for optimal performance. Generating appropriated index (using virtual offsets defined by the BGZF compression).")
        gaf_file = libcbgzf.BGZFile(gaf_path,"rb")
    else:
        logger.info("Uncompressed GAF file detected. Generating appropriate index (using offset values).")
        gaf_file = open(gaf_path,"rt")
    out_dict = {}
    logger.info("Indexing the file")
    offset = 0
    while True:
        offset = gaf_file.tell()
        mapping = gaf_file.readline()
        if not mapping:
            break
        try:
            val = mapping.rstrip().split('\t')
        except TypeError:
            val = mapping.decode("utf-8").rstrip().split('\t')
        if not unstable:
            with timers("convert_coord"):
                alignment = convert_coord(val, ref)
        else:
            alignment = list(re.split('>|<', val[5]))[1:]
        for a in alignment:
            try:
                out_dict[nodes[a]].append(offset)
            except KeyError:
                out_dict[nodes[a]] = [offset]
    out_dict["ref_contig"] = ref_contig
    
    gaf_file.close()
        
    with open(output, 'wb') as handle:
        with timers("write_file"):
            pickle.dump(out_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Time spent sorting GFA:                      %9.2f s", timers.elapsed("sort_gfa"))
    logger.info("Time spent storing contigs:                  %9.2f s", timers.elapsed("store_contig_info"))
    logger.info("Time spent converting coordinates:           %9.2f s", timers.elapsed("convert_coord"))
    logger.info("Total time:                                  %9.2f s", total_time)

def convert_coord(line, ref):
    
    from gaftools.cli.convert import search_intervals

    unstable_coord = []
    gaf_contigs = list(filter(None, re.split('(>)|(<)', line[5])))
    for nd in gaf_contigs:
        if nd == ">" or nd == "<":
            continue
        if ':' in nd and '-' in nd:
            tmp = nd.rstrip().split(':')
            query_contig_name = tmp[0]
            (query_start, query_end) = tmp[1].rstrip().split('-')
        else:
            query_start = line[7]
            query_end = line[8] 
            query_contig_name = nd
        
        #print(orient, query_contig_name, query_start, query_end)
        '''Find the matching nodes from the reference genome here'''
        start, end = search_intervals(ref[query_contig_name], int(query_start), int(query_end), 0, len(ref[query_contig_name]))

        for i in ref[query_contig_name][start:end+1]:
            cases = -1
            if i.start <= int(query_start) < i.end:
                cases = 1
            elif i.start < int(query_end) <= i.end:
                cases = 2
            elif int(query_start) < i.start < i.end < int(query_end):
                cases = 3
            
            if cases != -1:    
                unstable_coord.append(i.node_id)
    
    return unstable_coord

def detect_unstable(line):
    a = line.rstrip().split('\t')[5]
    g = list(filter(None, re.split('(>)|(<)', a)))
    if len(g) == 1:
        return False
    if ":" in g[1]:
        return False
    return True

# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_path', metavar='GAF', help='Input GAF file (can be gzip-compressed)')
    arg('reference', metavar='rGFA', help='Reference rGFA file has to be input.')
    arg('-o', '--output', default=None, help='Output Indexed GAF file. If omitted, use <GAF File>.gai.')
    
# fmt: on
def validate(args, parser):
    return True


def main(args):
    run(**vars(args))
