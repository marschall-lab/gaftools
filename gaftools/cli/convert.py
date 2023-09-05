"""
Convert Coordinate Systems between the stable system and unstable system
"""

import logging
import sys
import gaftools.utils as utils
import gaftools.gaf

from gaftools import __version__
from gaftools.cli import log_memory_usage
from gaftools.cli import CommandLineError
from gaftools.timer import StageTimer
from gaftools.gaf import parse_gaf


logger = logging.getLogger(__name__)

class Node1:
    def __init__(self, node_id, start, end):
        self.node_id = node_id
        self.start = start
        self.end = end


class Node2:
    def __init__(self, contig_id, start, end):
        self.contig_id = contig_id
        self.start = start
        self.end = end
    
    def to_string(self, orient):
        return "%s%s:%d-%d"%(orient,self.contig_id,self.start,self.end)


def run(gaf_file, gfa_file, output=sys.stdout, unstable=False, stable=False):

    timers = StageTimer()
    assert (unstable != stable)
    if (unstable):
        stable_to_unstable(gaf_file, gfa_file, output)
    else:
        unstable_to_stable(gaf_file, gfa_file, output)
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


def merge_nodes(node1, node2, orient1, orient2):
    
    if (node1.contig_id != node2.contig_id) or (orient1 != orient2):
        return False
    if (orient1 == ">") and (node1.end != node2.start):
        return False
    if (orient1 == "<") and (node1.start != node2.end):
        return False
    if (orient1 == "<"):
        node = Node2(node1.contig_id, node2.start, node1.end)
    else:
        node = Node2(node1.contig_id, node1.start, node2.end)
    return [node, orient1]


def stable_to_unstable(gaf_path, gfa_path, out_path):
    '''This function converts a GAF file (mappings to a pangenome graph) into unstable coordinate. It does not expect sorted
    input however it strictly assumes that SO, LN and SN tags are available in the rGFA...
    '''

    import re
    import copy
    import gzip
    import itertools
    import gaftools.utils as utils
    from gaftools.cli.sort import gfa_sort_basic
    
    '''Needs to sort the gfa to use logn time binary search'''
    logger.info("INFO: Sorting the GFA file...")
    gfa_lines = gfa_sort_basic(gfa_path)
    print("Done")  

    '''We load the GFA into memory for fast execution. GFA is not very large
    so it does not seem to be a big issue... This creates a dictionary where each element is a
    contig that keeps the list of start and end locations with node name(S).
    '''
    logger.info("INFO: Loading the rGFA file into memory...")
    reference = {}    
    contig_name = None
    for gfa_line in gfa_lines:
        tmp_contig_name = [k for k in gfa_line if k.startswith("SN:Z:")][0][5:]
        
        if tmp_contig_name != contig_name:
            contig_name = copy.deepcopy(tmp_contig_name)
            if contig_name not in reference:
                reference[contig_name] = []

        start_pos = int([k for k in gfa_line if k.startswith("SO:i:")][0][5:])
        end_pos = int([k for k in gfa_line if k.startswith("LN:i:")][0][5:]) + start_pos
        tmp = Node1(gfa_line[1], start_pos, end_pos)
        reference[contig_name].append(tmp)
    
    gaf_unstable = open(out_path, "w")

    logger.info("INFO: Reading and converting the alignments...")
     
    line_count = 0

    for gaf_line in parse_gaf(gaf_path):
        
        gaf_contigs = list(filter(None, re.split('(>)|(<)', gaf_line.path)))
        assert len(gaf_contigs) >= 1

        unstable_coord = ""
        orient = None

        new_start = -1
        new_total = 0
        for nd in gaf_contigs:
            if nd == ">" or nd == "<":
                orient = nd
                continue
            if ':' in nd and '-' in nd:
                tmp = nd.rstrip().split(':')
                query_contig_name = tmp[0]
                (query_start, query_end) = tmp[1].rstrip().split('-')
                split_contig = True
            else:
                query_start = gaf_line.path_start
                query_end = gaf_line.path_end
                query_contig_name = nd
                split_contig = False
            if not orient:
                if gaf_line.strand == "+":
                    orient = ">"
                else:
                    orient = "<"
 
            '''Find the matching nodes from the reference genome here'''
            start, end = utils.search_intervals(reference[query_contig_name], query_start, query_end, 0, len(reference[query_contig_name]))
 
            nodes_tmp = []
            for i in reference[query_contig_name][start:end+1]:
                cases = -1
                if i.start <= int(query_start) < i.end:
                    cases = 1
                    if new_start == -1:
                        if split_contig:
                            new_start = int(query_start)
                        else:
                            new_start = int(query_start) - i.start
                elif i.start < int(query_end) <= i.end:
                    cases = 2
                elif int(query_start) < i.start < i.end < int(query_end):
                    cases = 3
                
                if cases != -1:
                    nodes_tmp.append(i.node_id)
                    new_total += (i.end - i.start)
            
            if orient == "<":
                for i in reversed(nodes_tmp):
                    unstable_coord += orient+i
            else:
                for i in nodes_tmp:
                    unstable_coord += orient+i


        if line_count != 0:
            gaf_unstable.write("\n")
        
        if gaf_line.strand == "-":
            new_end = new_total - new_start
            new_start = new_end - (gaf_line.path_end - gaf_line.path_start)
        else:
            new_end = new_start + (gaf_line.path_end - gaf_line.path_start)

        gaf_unstable.write("%s\t%s\t%s\t%s\t+\t%s\t%d\t%d\t%d\t%d\t%d\t%d"
                           %(gaf_line.query_name, gaf_line.query_length, gaf_line.query_start, gaf_line.query_end,
                            unstable_coord, new_total, new_start, new_end, gaf_line.residue_matches,
                            gaf_line.alignment_block_length, gaf_line.mapping_quality))
        

        line_count += 1
        """for i in gaf_line_elements[9:len(gaf_line_elements)-1]:
            gaf_unstable.write("\t%s"%i)"""

        for k in gaf_line.tags.keys():
            gaf_unstable.write("\t%s:%s"%(k,gaf_line.tags[k]))

        #Add cigar in reverse 
        if gaf_line.strand == "-":
            new_cigar = utils.reverse_cigar(gaf_line.cigar)
            gaf_unstable.write("\t%s" %new_cigar)
        else:
            gaf_unstable.write("\t%s" %gaf_line.cigar)

    gaf_file.close()
    gaf_unstable.close()


def unstable_to_stable(gaf_path, gfa_path, out_path):
    '''This function converts a gaf file (mappings to a pangenome graph). It does not need sorted
    input however it strictly assumes that SO, LN and SN tags are available...
    '''

    import re
    import gzip
    
    logger.info("INFO: Loading the rGFA file into memory")
    nodes = {}
    contig_len = {}
    ref_contig = []
    contig_name = None
    
    gz_flag = gfa_path[-2:] == "gz"
    if gz_flag:
        gfa_file = gzip.open(gfa_path,"r")
    else:
        gfa_file = open(gfa_path,"r")

    for gfa_line in gfa_file:
        if gz_flag:
            gfa_line = gfa_line.decode("utf-8")
        if gfa_line[0] != "S":
            break
        
        gfa_line = gfa_line.rstrip().split('\t')
        contig_name = [k for k in gfa_line if k.startswith("SN:Z:")][0][5:]
        start_pos = int([k for k in gfa_line if k.startswith("SO:i:")][0][5:])
        end_pos = int([k for k in gfa_line if k.startswith("LN:i:")][0][5:]) + start_pos
        try:
            rank = int([k for k in gfa_line if k.startswith("SR:i:")][0][5:])
        except IndexError:
            raise CommandLineError("No Rank present in the reference GFA File. Input rGFA file should have SR field.")
        if contig_name not in ref_contig:
            if rank == 0:
                ref_contig.append(contig_name)
        else:
            assert (rank == 0)

        tmp = Node2(contig_name, start_pos, end_pos)
        nodes[gfa_line[1]] = tmp

        try:
            contig_len[contig_name] += end_pos - start_pos
        except KeyError:
            contig_len[contig_name] = end_pos-start_pos
    gfa_file.close()
    

    gaf_stable = open(out_path, "w")
    logger.info("INFO: Reading and converting the alignments...")
    line_count = 0
    
    new_total = None
    new_start = None
    for gaf_line in parse_gaf(gaf_path):
        reverse_flag = False
        gaf_nodes = list(filter(None, re.split('(>)|(<)', gaf_line.path)))
        node_list = []
        stable_coord = ""
        orient = None
        
        for nd in gaf_nodes:
            if nd == ">" or nd == "<":
                orient = nd
                continue
            if not orient:
                orient = ">"
            node_list.append([nodes[nd], orient])
        out_node = [node_list[0]]

        for i in range(len(node_list)-1):
            n1 = out_node[-1][0]
            o1 = out_node[-1][1]
            n2 = node_list[i+1][0]
            o2 = node_list[i+1][1]
            node_merge = merge_nodes(n1, n2, o1, o2)
            if node_merge == False:
                stable_coord += n1.to_string(o1)
                out_node.append([n2, o2])
            else:
                out_node[-1] = node_merge
        if len(out_node) == 1 and out_node[0][0].contig_id in ref_contig:
            if out_node[0][1] == "<":
                reverse_flag = True
                gaf_line.strand = "-"
                new_start = out_node[0][0].start + gaf_line.path_length - gaf_line.path_end
            else:
                new_start = out_node[0][0].start + gaf_line.path_start
            stable_coord = out_node[0][0].contig_id
            new_total = contig_len[stable_coord]
        else:
            stable_coord += out_node[-1][0].to_string(out_node[-1][1])
            new_start = gaf_line.path_start
            new_total = gaf_line.path_length

        if line_count != 0:
            gaf_stable.write("\n")

        gaf_stable.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" %(gaf_line.query_name, 
                        gaf_line.query_length, gaf_line.query_start, gaf_line.query_end, gaf_line.strand, 
                        stable_coord, new_total, new_start, new_start + gaf_line.path_start - gaf_line.path_end, 
                        gaf_line.residue_matches, gaf_line.alignment_block_length, gaf_line.mapping_quality))

        line_count += 1

        for k in gaf_line.tags.keys():
            gaf_stable.write("\t%s:%s"%(k,gaf_line.tags[k]))
        
        #Add cigar in reverse 
        if reverse_flag:
            tmp_cigar = utils.reverse_cigar(gaf_line.cigar)
            gaf_stable.write("\t%s" %tmp_cigar)
        else:
            gaf_stable.write("\t%s" %gaf_line.cigar)

    gaf_file.close()
    gaf_stable.close()
    return True

# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_file', metavar='GAF', help='GAF File whose coordinates have to be changed')
    arg('gfa_file', metavar='rGFA', help='Input rGFA file to convert the coordinates')
    arg('-o', '--output', default=sys.stdout,
        help='Output GAF file. If omitted, use standard output.')
    arg('--unstable', dest='unstable', default=False, action='store_true',
        help='Convert to Unstable Coordinates')
    arg('--stable', dest='stable', default=False, action='store_true',
        help='Convert to Stable Coordinates')


# fmt: on
def validate(args, parser):
    if not args.unstable and not args.stable:
        parser.error("Either --unstable or --stable has to be provided")
    if args.unstable and args.stable:
        parser.error("Specify --unstable or --stable; both cannot be provided")


def main(args):
    run(**vars(args))
