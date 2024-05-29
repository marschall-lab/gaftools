"""
View the GAF File based on parameters
"""
import logging
import pickle
import os
import copy
import re
import sys

from gaftools import __version__
from gaftools.cli import log_memory_usage, CommandLineError
from gaftools.utils import search_intervals
from gaftools.timer import StageTimer
from gaftools.gaf import GAF
from gaftools.conversion import Node1, stable_to_unstable, unstable_to_stable, making_reference_object, read_gfa_unstable_to_stable, to_stable, to_unstable


logger = logging.getLogger(__name__)

def run(gaf_path,
        gfa=None,
        output=None,
        index=None, 
        nodes=[],
        regions=[],
        format=None
    ):
    
    timers = StageTimer()
    
    if output == None:
        writer = sys.stdout
    else:
        writer = open(output, 'w')
    # Need to detect the format of the input gaf
    gaf=GAF(gaf_path)
    gaf_format=None
    # checking format in the first 10 lines.
    for i, gaf_line in enumerate(gaf.read_file()):
        if i == 10:
            break
        if i == 0:
            gaf_format=gaf_line.detect_path_format()
        assert gaf_format==gaf_line.detect_path_format()
    gaf.close()

    # if format is given, prepare some objects for use later
    if format:
        if format == 'stable':
            if gaf_format == True:
                raise CommandLineError('Input GAF already has stable coordinates. Please remove the --format stable option')
            # TODO: Need to incorporate these objects into the GFA class
            gfa_nodes, contig_len, ref_contig = read_gfa_unstable_to_stable(gfa)
        else:
            assert format == 'unstable'
            if gaf_format == False:
                raise CommandLineError('Input GAF already has unstable coordinates. Please remove the --format unstable option')
            # TODO: Need to incorporate these objects into the GFA class
            reference = making_reference_object(gfa)        

    # now find out what lines to view and how to view
    if len(nodes) != 0 or len(regions) != 0:
        if index == None:
            index = gaf_path+".gvi"
            if not os.path.exists(index):
                raise CommandLineError("No index found. Please provide the path to the index or create one with gaftools index.")

        ind = None
        with open(index, 'rb') as tmp:
            ind = pickle.load(tmp)
        
        ind_key = sorted(list(ind.keys()), key = lambda x: (x[1], x[2]))
        ind_dict = {}
        for i in ind_key:
            ind_dict[i[0]] = i
        
        if regions:
            assert nodes == []
            nodes = get_unstable(regions, ind)
        offsets=ind[ind_dict[nodes[0]]]
        for nd in nodes[1:]:
            # extracting all the lines that touches at least one of the nodes
            offsets = list(set(offsets) | set(ind[ind_dict[nd]]))
        offsets.sort()
        if len(offsets) == 0:
            raise CommandLineError("No alignments found for the given nodes/regions")
        gaf = GAF(gaf_path)
        # if format specified, have to make the changes.
        if format:
            if format == 'stable':
                for ofs in offsets:
                    line = gaf.read_line(ofs)
                    print(to_stable(line, gfa_nodes, ref_contig, contig_len), file=writer)
            else:
                assert format == 'unstable'
                for ofs in offsets:
                    line = gaf.read_line(ofs)
                    print(to_unstable(line, reference), file=writer)
        # if no format given, then just print the selected lines
        else:
            for ofs in offsets:
                line = gaf.read_line(ofs)
                print(line, file=writer)
        gaf.close()
    else:
        # No nodes or regions indicates the entire file will be viewed
        # converting the format of the entire file
        if format:
            if format == 'stable':
                for line in unstable_to_stable(gaf_path, gfa_nodes, ref_contig, contig_len):
                    print(line, file=writer)
            else:
                for line in stable_to_unstable(gaf_path, reference):
                    print(line, file=writer)
        else:
            # No format also given. So just need to print the file.
            gaf = GAF(gaf_path)
            for line in gaf.file:
                if gaf.gz_flag:
                    print(line.decode("utf-8").rstrip(), file=writer)
                else:
                    print(line.rstrip(), file=writer)

    
    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Total time:                                  %9.2f s", total_time)


def change_format(a, show_node_id, node_id, ind, ind_dict, fa):
    
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


def get_unstable(regions, index):
    """ Takes the regions and returns the node IDs """

    contig = [x.split(":")[0] for x in regions]
    node_dict = {}
    start = [x.split(":")[1].split("-")[0] for x in regions]
    end = [x.split(":")[1].split("-")[-1] for x in regions]
    
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
            logger.info("INFO: Region %s spans multiple nodes.\nThe nodes are:"%(node[n]))
            for n in node:
                logger.info("INFO: %s\t%s\t%d\t%d"%(n[0],n[1],n[2],n[3]))
        
        result.append(node[0][0])
    
    return result    

def search(node, node_list):
    """ Find the unstable node id from the region """
    
    s = 0
    pos = 0
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
    # if there is only one node for the entire contig (case for non-reference nodes)
    # then the above loop is not executed and we extract the only node with pos=0
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
    arg('-g', '--gfa', dest='gfa', metavar='GFA', default=None, help='Input GFA file (can be gzip-compressed). Required when converting from one coordinate system to another.')
    arg('-o', '--output', dest='output', metavar='OUTPUT', default=None, help='Output file. Default is stdout.')
    arg('-i', '--index', default=None, help='Path to GAF Index file. This index is created using gaftools index. If path is not provided, it is assumed to be in the same directory as GAF file with the same name and .gvi extension (default location of the index script).\nIt is required if nodes or regions are specified.')
    arg('-n', '--node', dest='nodes', metavar='NODE', default=[], action='append', help='Nodes to search for. Multiple can be provided.')
    arg('-r', '--region', dest='regions', metavar='REGION', default=[], action='append', help='Regions to search for. Multiple can be provided.')
    arg('-f', '--format', dest='format', metavar='FORMAT', help='format of output path (unstable | stable)')
    
# fmt: on

def validate(args, parser):
    if args.format and (args.format not in ['unstable', 'stable']):
        parser.error("--format only accepts unstable or stable as input.")
    if args.nodes and args.regions:
        parser.error("provide either of the --regions and --nodes options and not both.")
    if args.format and not args.gfa:
        parser.error("GFA file has to be provided along with --format.")

def main(args):
    run(**vars(args))
