class Ref:
    def __init__(self, node_id, start, end):
        self.node_id = node_id
        self.start = start
        self.end = end


def main():
    import argparse

    parser = argparse.ArgumentParser(description="gaftools stable-to-unstable")
    parser.add_argument("-a", "--alignments", dest="gaf_file",
                        help="The gaf file with stable coordinate system to be converted to unstable.") 
    parser.add_argument("-g", "--genome", dest="gfa_file",
                        help="Pangenome graph.")
    parser.add_argument("-o", "--out", dest="out_file",
                        help="The name of the output file (Write to stdout by default).")

    args = parser.parse_args()
    if not args.gfa_file:
        print("Please input the gfa file to be converted using --genome (or -g)...")
        return
    if not args.gaf_file:
        print("Please input the gaf file to be converted using --alignments (or -a)...")
        return
    if not args.out_file:
        print("Please input the output file path using --out (or -o)...")
        return

    stable_to_unstable(args.gfa_file, args.gaf_file, args.out_file)


def stable_to_unstable(gfa_path, gaf_path, out_path):
    '''This function converts a gaf file (mappings to a pangenome graph). It does not expect sorted
    input however it strictly assumes that SO, LN and SN tags are available...
    '''

    import re
    import copy
    from sort_gfa import gfa_sort
    
    '''Needs to sort the gfa to use logn time binary search'''
    print("Sorting the GFA file")
    gfa_lines = gfa_sort(gfa_path, None, True)
    
    '''We load the reference genome into memory for fast execution. The reference is not very large
    so it does not seem to be a big issue... This creates a dictionary where each element is a
    contig which keeps the list of start and end locations with node name(S).
    '''
    print("Loading the GFA file into memory")
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
        tmp = Ref(gfa_line[1], start_pos, end_pos)
        reference[contig_name].append(tmp)
    print()
    
    gaf_unstable = open(out_path, "w")

    print("Reading the alignments...")
    line_count = 0
    with open(gaf_path, "r") as gaf_file:
        for gaf_line in gaf_file:
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
                    
                    #if query_contig_name == "NA20129#2#JAHEPD010000261.1":
                    #    print(i.start, i.end, i.node_id, cases)

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

    gaf_unstable.close()
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
    

if __name__ == "__main__":
    main()

