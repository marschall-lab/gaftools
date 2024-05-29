"""
Index the GAF File
"""

import gzip
import copy
import re
import pickle
from pysam import libcbgzf
import logging

import gaftools.utils as utils
from gaftools.gaf import GAF
from gaftools import __version__
from gaftools.timer import StageTimer
from gaftools.cli import log_memory_usage, CommandLineError


class Node1:
    def __init__(self, node_id, start, end):
        self.node_id = node_id
        self.start = start
        self.end = end

logger = logging.getLogger(__name__)


def run(gaf_path, 
        gfa_path, 
        output=None
       ):

    timers = StageTimer()
    if output == None:
        output = gaf_path+".gvi"
    
    #Detecting if GAF has stable or unstable coordinate
    if utils.is_file_gzipped(gaf_path):
        gaf_file = libcbgzf.BGZFile(gaf_path,"rb")
        line = gaf_file.readline().decode("utf-8")
    else:
        gaf_file = open(gaf_path,"rt")
        line = gaf_file.readline()
    
    unstable = utils.is_unstable(line)
    gaf_file.close()
    
    nodes = {}
    ref = {}
    ref_contig = []
    if not unstable:
        # TODO: Need to remove this gfa sorting part. It can probably be incorporated into the GFA class
        with timers("sort_gfa"):
            logger.info("INFO: Sorting GFA File")
            gfa_lines = utils.gfa_sort_basic(gfa_path)
        contig_name = None
        
        with timers("store_contig_info"):
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
                    logger.error("No Rank present in the reference GFA File. Input rGFA file should have SR field.")
                    exit()
                ref[contig_name].append(tmp)
                nodes[gfa_line[1]] = (gfa_line[1], tmp_contig_name, start_pos, end_pos)
                if tmp_contig_name not in ref_contig:
                    if rank == 0:
                        ref_contig.append(contig_name)
                else:
                    assert (rank == 0)
    else:
        gz_flag = gfa_path[-2:] == "gz"
        if gz_flag:
            gfa_file = gzip.open(gfa_path,"r")
        else:
            gfa_file = open(gfa_path,"r")
        with timers("store_contig_info"):
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
                    raise CommandLineError("ERROR: No Rank present in the reference GFA File. Input rGFA file should have SR field.")
                nodes[gfa_line[1]] = (gfa_line[1], contig_name, start_pos, end_pos)
                if contig_name not in ref_contig:
                    if rank == 0:
                        ref_contig.append(contig_name)
                else:
                    assert (rank == 0)
    
    if utils.is_file_gzipped(gaf_path):
        gaf_file = libcbgzf.BGZFile(gaf_path,"rb")
    else:
        gaf_file = open(gaf_path,"rt")
    
    out_dict = {}
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
    logger.info("Time to sort gfa:                            %9.2f s", timers.elapsed('sort_gfa'))
    logger.info("Time to store contig info:                   %9.2f s", timers.elapsed('store_contig_info'))
    logger.info("Total time:                                  %9.2f s", total_time)


def convert_coord(line, ref):

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
        
        '''Find the matching nodes from the reference genome here'''
        start, end = utils.search_intervals(ref[query_contig_name], int(query_start), int(query_end), 0, len(ref[query_contig_name]))

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


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_path', metavar='GAF', help='Input GAF file (can be gzip-compressed)')
    arg('gfa_path', metavar='rGFA', help='Reference rGFA file has to be input.')
    arg('-o', '--output', default=None, help='Output Indexed GAF file. If omitted, use <GAF File>.gvi.')
    
# fmt: on
def validate(args, parser):
    return True


def main(args):
    run(**vars(args))
