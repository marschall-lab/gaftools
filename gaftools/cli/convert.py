"""
Convert Coordinate Systems between the stable system and unstable system
"""

import logging
import sys

from gaftools import __version__
from gaftools.cli import log_memory_usage
from gaftools.cli import CommandLineError


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

def run(
    gaf_file,
    gfa_file,
    output=sys.stdout,
    unstable=False,
    stable=False
):
    assert (unstable != stable)
    if (unstable):
        stable_to_unstable(gaf_file, gfa_file, output)
    else:
        unstable_to_stable(gaf_file, gfa_file, output)


def stable_to_unstable(gaf_path, gfa_path, out_path):
    '''This function converts a gaf file (mappings to a pangenome graph). It does not expect sorted
    input however it strictly assumes that SO, LN and SN tags are available...
    '''

    import re
    import copy
    from gaftools.cli.sort_gfa import gfa_sort
    import gzip
    
    '''Needs to sort the gfa to use logn time binary search'''
    print("Sorting the GFA file...")
    gfa_lines = gfa_sort(gfa_path, None, True)
    
    '''We load the reference genome into memory for fast execution. The reference is not very large
    so it does not seem to be a big issue... This creates a dictionary where each element is a
    contig which keeps the list of start and end locations with node name(S).
    '''
    #print("Loading the GFA file into memory")
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
    #print()
    
    gaf_unstable = open(out_path, "w")

    print("Reading the alignments...")
    gz_flag = gaf_path[-2:] == "gz"
    if gz_flag:
        gaf_file = gzip.open(gaf_path,"r")
    else:
        gaf_file = open(gaf_path,"r")
    
    line_count = 0
    for gaf_line in gaf_file:
        if gz_flag:
            gaf_line_elements = gaf_line.decode("utf-8").rstrip().split('\t')
        else:
            gaf_line_elements = gaf_line.rstrip().split('\t')
        gaf_contigs = list(filter(None, re.split('(>)|(<)', gaf_line_elements[5])))
        unstable_coord = ""
        orient = None

        new_start = -1
        new_total = 0
        for nd in gaf_contigs:
            if nd == ">" or nd == "<":
                orient = nd
                continue
            if ':' in nd:
                tmp = nd.rstrip().split(':')
                query_contig_name = tmp[0]
                (query_start, query_end) = tmp[1].rstrip().split('-')
            else:
                query_start = gaf_line_elements[7]
                query_end = gaf_line_elements[8] 
                query_contig_name = nd
            if not orient:
                orient = ">"

            #print(orient, query_contig_name, query_start, query_end)
            '''Find the matching nodes from the reference genome here'''
            start, end = search_intervals(reference[query_contig_name], int(query_start), int(query_end), 0, len(reference[query_contig_name]))

            for i in reference[query_contig_name][start:end+1]:
                cases = -1
                if i.start <= int(query_start) < i.end:
                    cases = 1
                    if new_start == -1:
                        new_start = int(query_start) - i.start
                elif i.start < int(query_end) <= i.end:
                    cases = 2
                elif int(query_start) < i.start < i.end < int(query_end):
                    cases = 3
                
                if cases != -1:    
                    unstable_coord += orient+i.node_id
                    new_total += (i.end - i.start)
        
        if line_count != 0:
            gaf_unstable.write("\n")
        gaf_unstable.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d" %(gaf_line_elements[0], gaf_line_elements[1], gaf_line_elements[2], 
                                                                    gaf_line_elements[3],
                                                                    gaf_line_elements[4],
                                                                    unstable_coord, new_total,
                                                                    new_start, new_start +
                                                                    int(gaf_line_elements[9])))
        line_count += 1
        for i in gaf_line_elements[9:len(gaf_line_elements)]:
            gaf_unstable.write("\t%s"%i)
    gaf_file.close()
    gaf_unstable.close()
    print("Done...")


def unstable_to_stable(gaf_path, gfa_path, out_path):
    '''This function converts a gaf file (mappings to a pangenome graph). It does not need sorted
    input however it strictly assumes that SO, LN and SN tags are available...
    '''

    import re
    import gzip
    
    print("Loading the GFA file into memory")
    nodes = {}
    contig_len = {}
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
        tmp = Node2(contig_name, start_pos, end_pos)
        nodes[gfa_line[1]] = tmp

        try:
            contig_len[contig_name] += end_pos-start_pos
        except KeyError:
            contig_len[contig_name] = end_pos-start_pos
    gfa_file.close()
    
    gaf_stable = open(out_path, "w")

    print("Reading the alignments...")
    line_count = 0
    gz_flag = gaf_path[-2:] == "gz"
    if gz_flag:
        gaf_file = gzip.open(gaf_path,"r")
    else:
        gaf_file = open(gaf_path,"r")
    
    new_total = None
    new_start = None
    for gaf_line in gaf_file:
        if gz_flag:
            gaf_line_elements = gaf_line.decode("utf-8").rstrip().split('\t')
        else:
            gaf_line_elements = gaf_line.rstrip().split('\t')
        gaf_nodes = list(filter(None, re.split('(>)|(<)', gaf_line_elements[5])))
        node_list = []
        stable_coord = ""
        orient = None
        for nd in gaf_nodes:
            if nd == ">" or nd == "<":
                orient = nd
                continue
            if not orient:
                orient = ">"
            node_list.append([nodes[nd],orient])
        for i in range(len(node_list)-1):
            n1=node_list[i][0]
            o1=node_list[i][1]
            n2=node_list[i+1][0]
            o2=node_list[i+1][1]
            node_merge=merge_nodes(n1,n2,o1,o2)
            if node_merge==False:
                stable_coord += n1.to_string(o1)
            else:
                node_list[i+1]=node_merge
        if len(node_list) == 1:
            stable_coord = node_list[0][0].contig_id
            new_total = contig_len[stable_coord]
            new_start = node_list[0][0].start + int(gaf_line_elements[7])
        else:
            stable_coord += node_list[-1][0].to_string(node_list[-1][1])
            new_start = int(gaf_line_elements[7])
            new_total = int(gaf_line_elements[6])
        
        if line_count != 0:
            gaf_stable.write("\n")
        gaf_stable.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d" %(gaf_line_elements[0], gaf_line_elements[1], gaf_line_elements[2], 
                                                                    gaf_line_elements[3],
                                                                    gaf_line_elements[4],
                                                                    stable_coord, new_total,
                                                                    new_start, new_start +
                                                                    int(gaf_line_elements[9])))
        line_count += 1
        for i in gaf_line_elements[9:len(gaf_line_elements)]:
            gaf_stable.write("\t%s"%i)
            
    gaf_file.close()
    gaf_stable.close()
    return True


def search_intervals(intervals, query_start, query_end, start, end):
    '''Given the start-end coordinates in the GFA file for the given contig (SO, SO+LN), it
    searches for the given (query_start, query_end) matches. (query_start, query_end) is the start
    and end location of a mapping in the gaf file.
    '''

    if start <= end:    
        mid = start + (end - start) // 2
        if query_end <= intervals[mid].start:
            return search_intervals(intervals, query_start, query_end, start, mid - 1)
        elif query_start >= intervals[mid].end:
            return search_intervals(intervals, query_start, query_end, mid + 1, end)
        else:
            return start, end

    return -1, -1

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

# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_file', metavar='GAF', help='GAF File whose coordinates have to be changed')
    arg('gfa_file', metavar='GFA', help='Input GFA file to conver the coordinates')
    arg('-o', '--output', default=sys.stdout,
        help='Output GAF file. If omitted, use standard output.')
    arg('--unstable', dest='unstable', default=False, action='store_true',
        help='Convert to Unstable Coordinates')
    arg('--stable', dest='stable', default=False, action='store_true',
        help='Convert to Stable Coordinates')
    
# fmt: on
def validate(args, parser):
    if not args.unstable and not args.stable:
        parser.error("Either one of --unstable or --stable has to be provided")
    if args.unstable and args.stable:
        parser.error("Specify one of --unstable or --stable. Both cannot be provided")

def main(args):
    run(**vars(args))
