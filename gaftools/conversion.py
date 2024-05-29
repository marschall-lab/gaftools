"""
Core scripts for converting coordinates.
"""
import logging
import re

import gaftools.utils as utils
<<<<<<< HEAD
=======
from gaftools.cli import CommandLineError
>>>>>>> f00e349 (reshuffling code and making minor changes to the GAF class. View command has been altered to include the convert function.)
from gaftools.gaf import GAF


logger = logging.getLogger(__name__)

class Node1:
    def __init__(self, node_id, start, end):
        self.node_id = node_id
        self.start = start
        self.end = end

# Stable Node representation of the Nodes.
# Used in the coversion from unstable to stable coordinate system.
class StableNode:
    def __init__(self, contig_id, start, end):
        self.contig_id = contig_id
        self.start = start
        self.end = end
    
    def to_string(self, orient):
        return "%s%s:%d-%d"%(orient,self.contig_id,self.start,self.end)


def merge_nodes(node1, node2, orient1, orient2):
    
    if (node1.contig_id != node2.contig_id) or (orient1 != orient2):
        return False
    if (orient1 == ">") and (node1.end != node2.start):
        return False
    if (orient1 == "<") and (node1.start != node2.end):
        return False
    if (orient1 == "<"):
        node = StableNode(node1.contig_id, node2.start, node1.end)
    else:
        node = StableNode(node1.contig_id, node1.start, node2.end)
    return [node, orient1]


def stable_to_unstable(gaf_path, reference):
    '''
    This function converts a GAF file (mappings to a pangenome graph) into unstable coordinate.
    It does not expect sorted input however it strictly assumes that SO, LN and SN tags are available in the rGFA.
    '''
    gaf_input = GAF(gaf_path)
    for gaf_line in gaf_input.read_file():
        yield to_unstable(gaf_line, reference)
    gaf_input.close()

def unstable_to_stable(gaf_path, nodes, ref_contig, contig_len):
    '''This function converts a gaf file (mappings to a pangenome graph). It does not need sorted
    input however it strictly assumes that SO, LN and SN tags are available...
    '''
    gaf_input = GAF(gaf_path)
    for gaf_line in gaf_input.read_file():
        yield to_stable(gaf_line, nodes, ref_contig, contig_len)
    gaf_input.close()


# separate function for converting lines to unstable coordinate
def to_unstable(gaf_line, reference):
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
        start, end = utils.search_intervals(reference[query_contig_name], int(query_start), int(query_end), 0, len(reference[query_contig_name]))
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
    
    if gaf_line.strand == "-":
        if split_contig:
            new_total = gaf_line.path_length
            new_start = new_total - gaf_line.path_end
            new_end = new_total - gaf_line.path_start
        else:
            new_end = new_total - new_start
            new_start = new_end - (gaf_line.path_end - gaf_line.path_start)
    else:
        if split_contig:
            new_total = gaf_line.path_length
            new_start = gaf_line.path_start
            new_end = gaf_line.path_end
        else:
            new_end = new_start + (gaf_line.path_end - gaf_line.path_start)

    new_line = "%s\t%s\t%s\t%s\t+\t%s\t%d\t%d\t%d\t%d\t%d\t%d"%(gaf_line.query_name, gaf_line.query_length, gaf_line.query_start, gaf_line.query_end,
                        unstable_coord, new_total, new_start, new_end, gaf_line.residue_matches,
                        gaf_line.alignment_block_length, gaf_line.mapping_quality)
    
    #Add cigar in reverse 
    if gaf_line.strand == "-":
        new_cigar = utils.reverse_cigar(gaf_line.cigar)
        gaf_line.tags['cg:Z:'] = new_cigar

    for k in gaf_line.tags.keys():
        new_line += "\t%s%s"%(k,gaf_line.tags[k])

    return new_line


def to_stable(gaf_line, nodes, ref_contig, contig_len):
    reverse_flag = False
    new_total = None
    new_start = None
    gaf_nodes = list(filter(None, re.split('(>)|(<)', gaf_line.path)))
    node_list = []
    stable_coord = ""
    orient = None
    new_line=""
    
    for nd in gaf_nodes:
        if nd == ">" or nd == "<":
            orient = nd
            continue
        if not orient:
            orient = ">"
        node_list.append([nodes[nd], orient])
    out_node = [node_list[0]]

    for i in range(len(node_list) - 1):
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

    new_line+="%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" %(gaf_line.query_name, 
                    gaf_line.query_length, gaf_line.query_start, gaf_line.query_end, gaf_line.strand, 
                    stable_coord, new_total, new_start, new_start + gaf_line.path_end - gaf_line.path_start, 
                    gaf_line.residue_matches, gaf_line.alignment_block_length, gaf_line.mapping_quality)
    
    #Add cigar in reverse 
    if reverse_flag:
        new_cigar = utils.reverse_cigar(gaf_line.cigar)
        gaf_line.tags['cg:Z:'] = new_cigar
        
    # adding tags in the sequence it was found
    for k in gaf_line.tags.keys():
        new_line+="\t%s%s"%(k,gaf_line.tags[k])

<<<<<<< HEAD
    return new_line
=======
    return new_line 


# Required for converting stable coordinates to unstable coordinates
# TODO: This needs to go into the GFA object.
def making_reference_object(gfa_path):
    '''Needs to sort the gfa to use logn time binary search'''
    
    gfa_lines = utils.gfa_sort_basic(gfa_path)

    '''We load the GFA into memory for fast execution. GFA is not very large
    so it does not seem to be a big issue... This creates a dictionary where each element is a
    contig that keeps the list of start and end locations with node name(S).
    '''
    logger.info("INFO: Loading the rGFA file into memory...")
    # the reference object needs to be given as a input here and not made here
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
    
    return reference

# Required for converting unstable coordinates to stable coordinates
# TODO: This needs to go into the GFA object.
def read_gfa_unstable_to_stable(gfa_path):
    # reading GFA file
    
    import gzip

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
            continue
        
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

    return nodes, contig_len, ref_contig
>>>>>>> f00e349 (reshuffling code and making minor changes to the GAF class. View command has been altered to include the convert function.)
