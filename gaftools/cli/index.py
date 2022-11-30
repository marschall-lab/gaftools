"""
Index the GAF File
"""

import logging
import sys
import platform
import re

from gaftools import __version__
from gaftools.timer import StageTimer
from gaftools.cli import log_memory_usage
from gaftools.cli import CommandLineError

class Node1:
    def __init__(self, node_id, start, end):
        self.node_id = node_id
        self.start = start
        self.end = end

logger = logging.getLogger(__name__)

def run(gaf_path, reference, output=None, unstable=False):
    
    import gzip
    import copy
    from gaftools.cli.sort_gfa import gfa_sort
    import pickle

    timers = StageTimer()
    if output == None:
        output = gaf_path+".gai"
    nodes = {}
    ref = {}
    if not unstable:
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
                ref[contig_name].append(tmp)
                nodes[gfa_line[1]] = (gfa_line[1], tmp_contig_name, start_pos, end_pos)
    else:
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
            nodes[gfa_line[1]] = (gfa_line[1], contig_name, start_pos, end_pos)
    
    if is_file_gzipped(gaf_path):
       open_gaf = gzip.open
    else:
        open_gaf = open
    out_dict = {}
    logger.info("Indexing the file")
    with open_gaf(gaf_path, "rt") as gaf_file:
        for line_count, mapping in enumerate(gaf_file):
            val = mapping.rstrip().split('\t')
            if not unstable:
                with timers("convert_coord"):
                    alignment = convert_coord(val, ref)
            else:
                alignment = list(re.split('>|<', val[5]))[1:]
            for a in alignment:
                try:
                    out_dict[nodes[a]].append(line_count)
                except KeyError:
                    out_dict[nodes[a]] = [line_count]
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

def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_path', metavar='GAF', help='Input GAF file (can be gzip-compressed)')
    arg('reference', metavar='GFA', help='Reference GFA file has to be input.')
    arg('-o', '--output', default=None, help='Output Indexed GAF file. If omitted, use <GAF File>.gai.')
    
    arg('--unstable', dest='unstable', default=False, action='store_true', help='Specify flag if GAF file has unstable coordinates')
    
    
# fmt: on


def validate(args, parser):
    return True


def main(args):
    run(**vars(args))
