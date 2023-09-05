"""
View the GAF File based on parameters
"""
import logging
import sys
import platform

from gaftools import __version__
from gaftools.cli import log_memory_usage
from gaftools.cli import CommandLineError
from gaftools.timer import StageTimer


logger = logging.getLogger(__name__)

def run(gaf_path, 
        index=None, 
        output=sys.stdout,
        node=None,
        only_alignment=False,
        remove_cigar=False,
        remove_read_id=False,
        show_node_id=False,
        full_alignment=False
    ):
    
    from gaftools.utils import is_file_gzipped, search_intervals
    import pickle
    from pysam import libcbgzf

    timers = StageTimer()
    if output != sys.stdout:
        out = open(output, 'w')
    else:
        out = output

    if is_file_gzipped(gaf_path):
        logger.info("INFO: Compressed GAF file detected. The index provided must have been built on the compressed file.")
        gaf_file = libcbgzf.BGZFile(gaf_path,"rb")
    else:
        logger.info("INFO: Uncompressed GAF file detected. The index provided must have been built on the uncompressed file.")
        gaf_file = open(gaf_path,"rt")
    
    if len(node) != 0:
        if index == None:
            index = gaf_path+".gai"
        ind = None
        with open(index, 'rb') as tmp:
            ind = pickle.load(tmp)
        ind_key = sorted(list(ind.keys()), key = lambda x: (x[1], x[2]))
        ind_dict = {}
        for i in ind_key:
            ind_dict[i[0]] = i
        isRegion = detect_format(node)
        if isRegion:
            node = get_unstable(node, ind)
        node_id = [n[0] for n in node]
        if (len(node_id)==1):
            logger.info("INFO: One node ID recovered from the list of regions/nodes given. Output will contain entire alignments which contain that node.")
        elif full_alignment:
            logger.info("INFO: Multiple node IDs recovered from the list of regions/nodes given. Output will contain entire alignment since --full-alignment flag has been given.")
        else:
            logger.info("INFO: Multiple node IDs recovered from the list of regions/nodes given. Output will contain parts of the alignment which from the first node given to the last node given.")
        offsets=ind[node[0]]
        for nd in node[1:]:
            offsets = list(set(offsets) & set(ind[nd]))
        offsets.sort()
        c = 0
        for ofs in offsets:
            gaf_file.seek(ofs)
            mapping = gaf_file.readline()
            try:
                val = mapping.rstrip().split('\t')
            except TypeError:
                val = mapping.decode("utf-8").rstrip().split('\t')
            out_str = ""
            if not remove_read_id:
                out_str += "%s\t"%(val[0])
            if only_alignment:
                result = change_format(val, show_node_id, node_id, ind_key, ind_dict, full_alignment)
                out_str += result+"\n"
            else:
                for n, fd in enumerate(val[1:]):
                    if fd[:3] == "cg:" and remove_cigar:
                        continue
                    if n == 4:
                        fd = change_format(val, show_node_id, node_id, ind_key, ind_dict, full_alignment)
                    out_str += fd+"\t"
                out_str = out_str.strip("\t")
                out_str += "\n"
            out.write(out_str)        
                
        
    else:
        logger.info("INFO: No Nodes specified.")
        for ofs in offsets:
            gaf_file.seek(ofs)
            mapping = gaf_file.readline()
            val = mapping.rstrip().split('\t')
            out_str = ""
            for n,i in enumerate(val):
                if remove_read_id and n == 0:
                    continue
                if i[:3] == "cg:" and remove_cigar:
                    continue
                out_str += i+"\t"
            out_str = out_str.strip("\t")
            out_str += "\n"
            out.write(out_str)
    
    gaf_file.close()
    if output != sys.stdout:
        out.close()
    
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)

def change_format(a, show_node_id, node_id, ind, ind_dict, fa):
    import re
    x = list(filter(None, re.split('(>)|(<)', a[5])))
    isStable = False
    if len(x) == 1:
        isStable = True
    elif ":" in x[1]:
        isStable = True
    if isStable:
        x = convert_to_unstable(a, ind)
    tmp = None
    out_str = ""
    if len(node_id) == 1:
        
        if show_node_id:
            for nd in x:
                out_str += nd
        else:
            tmp = [x[0], list(ind_dict[x[1]])]
            for i,nd in enumerate(x):
                if i == 0 or i == 1:
                    continue 
                if nd == ">" or nd == "<":
                    orient = nd
                    continue
                n = ind_dict[nd]
                if orient == tmp[-2] and n[1] == tmp[-1][1]:
                    if orient == ">" and n[2] == tmp[-1][3]:
                        tmp[-1][3] = n[3]
                        continue
                    if orient == "<" and n[3] == tmp[-1][2]:
                        tmp[-1][2] = n[2]
                        continue
                tmp.append(orient)
                tmp.append(list(ind_dict[nd]))
            for nd in tmp:
                if nd == ">" or nd == "<":
                    out_str += nd
                else:
                    out_str += "%s:%d-%d"%(nd[1],nd[2],nd[3])
        return out_str
    if not fa:
        orient=None
        tmp = []
        toStart = False
        end = None
        for nd in x:
            if nd == ">" or nd == "<":
                orient = nd
                continue
            if nd not in node_id and not toStart:
                continue
            if nd in node_id:
                toStart=True
                tmp.extend([orient, list(ind_dict[nd])])
                end = len(tmp)
                continue
            tmp.extend([orient, list(ind_dict[nd])])
        result = [tmp[0], tmp[1]]
        orient=None
        for i,nd in enumerate(tmp[:end]):
            if i == 0 or i == 1:
                continue 
            if nd == ">" or nd == "<":
                orient = nd
                continue
            n = nd
            if orient == result[-2] and n[1] == result[-1][1]:
                if orient == ">" and n[2] == result[-1][3]:
                    result[-1][3] = n[3]
                    continue
                if orient == "<" and n[3] == result[-1][2]:
                    result[-1][2] = n[2]
                    continue
            result.append(orient)
            result.append(nd)
        for nd in result:
            if nd == ">" or nd == "<":
                out_str += nd
            else:
                out_str += "%s:%d-%d"%(nd[1],nd[2],nd[3])
    else:
        if show_node_id:
            for nd in x:
                out_str += nd
        else:
            tmp = [x[0], list(ind_dict[x[1]])]
            for i,nd in enumerate(x):
                if i == 0 or i == 1:
                    continue 
                if nd == ">" or nd == "<":
                    orient = nd
                    continue
                n = ind_dict[nd]
                if orient == tmp[-2] and n[1] == tmp[-1][1]:
                    if orient == ">" and n[2] == tmp[-1][3]:
                        tmp[-1][3] = n[3]
                        continue
                    if orient == "<" and n[3] == tmp[-1][2]:
                        tmp[-1][2] = n[2]
                        continue
                tmp.append(orient)
                tmp.append(list(ind_dict[nd]))
            for nd in tmp:
                if nd == ">" or nd == "<":
                    out_str += nd
                else:
                    out_str += "%s:%d-%d"%(nd[1],nd[2],nd[3])
    return out_str


def convert_to_unstable(x, ind):
    
    import copy
    from gaftools.cli.convert import Node1, search_intervals
    import re
    
    reference = {}    
    contig_name = None
    for n in ind:
        tmp_contig_name = n[1]
        
        if tmp_contig_name != contig_name:
            contig_name = copy.deepcopy(tmp_contig_name)
            if contig_name not in reference:
                reference[contig_name] = []

        start_pos = n[2]
        end_pos = n[3]
        tmp = Node1(n[0], start_pos, end_pos)
        reference[contig_name].append(tmp)
    gaf_contigs = list(filter(None, re.split('(>)|(<)', x[5])))
    unstable_coord = []
    orient = None

    for nd in gaf_contigs:
        if nd == ">" or nd == "<":
            orient = nd
            continue
        if ':' in nd and '-' in nd:
            tmp = nd.rstrip().split(':')
            query_contig_name = tmp[0]
            (query_start, query_end) = tmp[1].rstrip().split('-')
        else:
            query_start = x[7]
            query_end = x[8] 
            query_contig_name = nd
        if not orient:
            orient = ">"

        start, end = search_intervals(reference[query_contig_name], int(query_start), int(query_end), 0, len(reference[query_contig_name]))
        tmp = []
        for i in reference[query_contig_name][start:end+1]:
            cases = -1
            if i.start <= int(query_start) < i.end:
                cases = 1
            elif i.start < int(query_end) <= i.end:
                cases = 2
            elif int(query_start) < i.start < i.end < int(query_end):
                cases = 3
            if cases != -1:
                if orient == "<":
                    tmp.insert(0, orient)
                    tmp.insert(1, i.node_id)
                elif orient == ">":
                    tmp.append(orient)
                    tmp.append(i.node_id)    
        for t in tmp:
            unstable_coord.append(t)
            
    return unstable_coord
            

def detect_format(node):
    """
    Detect whether node IDs or regions are given 
    """

    if all(":" in n for n in node):
        return True
    if all(":" not in n for n in node):
        return False
    
    logger.warning("All the input nodes should be of the same format type: node ids or regions")


def get_unstable(nodes, index):
    """ Takes the regions and returns the node IDs """

    contig = [x.split(":")[0] for x in nodes]
    node_dict = {}
    start = [x.split(":")[1].split("-")[0] for x in nodes]
    end = [x.split(":")[1].split("-")[-1] for x in nodes]
    
    result = []
    for n, c in enumerate(contig):
        
        try:
            node_list = node_dict[c]
        except KeyError:
            node_list = list(filter(lambda x: (x[1] == c),list(index.keys())))
            node_list.sort(key = lambda x: x[2])
            node_dict[c] = node_list
        
        node = search([contig[n], start[n], end[n]], node_list)
        if len(node) > 1:
            logger.info("INFO: Region %s spans multiple nodes.\nThe nodes are:"%(nodes[n]))
            for n in node:
                logger.info("INFO: %s\t%s\t%d\t%d"%(n[0],n[1],n[2],n[3]))
        
        result.extend(node)
    
    return result    

def search(node, node_list):
    """ Find the unstable node id from the region """

    s = 0
    e = len(node_list)-1
    q_s = int(node[1])
    q_e = int(node[2])
    while (s != e):
        m = int((s+e)/2)
        if (q_s >= node_list[m][2]) and (q_s < node_list[m][3]):
            pos = m
            break
        elif (q_s >= node_list[m][3]):
            s = m+1
        else:
            e = m-1
        pos = s
    
    result = [node_list[pos]]
    while True:
        if q_e < node_list[pos][3]:
            break
        pos += 1
        result.append(node_list[pos])

    return result
    
# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('gaf_path', metavar='GAF', help='Input GAF file (can be gzip-compressed)')
    arg('-i', '--index', default=None, help='Path to GAF Index file. If not provided, it is assumed to be in the same directory as GAF file with the same name and .gaf.gai extension')
    arg('-o', '--output', default=sys.stdout, help='Output Indexed GAF file. If omitted, use <GAF File>.gai.')
    arg('--node', dest='node', metavar='NODE', default=[], action='append', help='Specify nodes to filter alignments. Instead of node ID, regions can also be specified. Can be used multiple times.\n When multiple nodes are specified, output contains partial alignment between first and last node. Entire alignment can be shown with --full-alignment flag.')
    arg('--only-alignment', dest='only_alignment', default=False, action='store_true', help='Show alignments which contain the list of nodes.')
    arg('--full-alignment', dest='full_alignment', default=False, action='store_true', help='Show the entire alignment with multiple nodes given.')
    arg('--remove-cigar', dest='remove_cigar', default=False, action='store_true', help='Option to remove CIGAR strings (if --only-alignment has not been chosen).')
    arg('--remove-read-id', dest='remove_read_id', default=False, action='store_true', help='Option to remove read IDs from output.')
    arg('--show-node-id', dest='show_node_id', default=False, action='store_true', help='Show list of nodes as node IDs. Option only available with --only-alignment (Default: Show in human readable form.)')
    
    
# fmt: on


def validate(args, parser):
    import os
    if not args.index and not os.path.exists(args.gaf_path+".gai"):
        parser.error("Could not find the index of the GAF file. Please provide one with --index.")
    if args.only_alignment and args.remove_cigar:
        parser.error("CIGAR is not shown with the --only-alignment flag. Specify one.")
    if args.node == [] and args.full_alignment:
        parser.error("The --full-alignment flag can only be given along with the --node option.")

def main(args):
    run(**vars(args))
