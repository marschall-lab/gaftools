"""
GFA sorting function
"""

import logging

from gaftools import __version__
from gaftools.utils import is_file_gzipped

logger = logging.getLogger(__name__)

def gfa_sort(gfa_path, out_path = None, return_list = True):
    '''This function sorts the given gfa file based on the contig name and start position within the
    contig. Note that it only sorts S lines and leaves the others.
    This can be called from the command line or from another funtion by providing "True" to the
    return_list argument.
    '''

    import glob
    from heapq import merge
    import functools
    import gzip
    import os


    gfa_lines = []
    path = "part*.gfa"
    chunk_size = 250000
    chunk_id = 1

    if is_file_gzipped(gfa_path):
       open_gfa = gzip.open
       is_gzipped = True
    else:
        open_gfa = open
    
    with open_gfa(gfa_path, 'rt') as gfa_file:
        f_out = open('part_{}.gfa'.format(chunk_id), 'w')
        
        for line_num, mapping in enumerate(gfa_file, 1):
            val = mapping.rstrip().split('\t')
            gfa_lines.append(val)
        
            if not line_num % chunk_size:
                gfa_lines.sort(key=functools.cmp_to_key(compare_gfa))
                
                for line_count, line in enumerate(gfa_lines):
                    f_out.write('\t'.join(line) + '\n') 
            
                logger.info('INFO: Splitting %d' %chunk_id)
                f_out.close()
                gfa_lines = []
                chunk_id += 1
                f_out = open('part_{}.gfa'.format(chunk_id), 'w')


        if gfa_lines:
            logger.info('INFO: Splitting %d' %chunk_id)
            gfa_lines.sort(key=functools.cmp_to_key(compare_gfa))
            for line_count, line in enumerate(gfa_lines):
                f_out.write('\t'.join(line) + '\n') 
            f_out.close()
            gfa_lines = []

    chunks = []
    for filename in glob.glob(path):
        chunks += [open(filename, 'r')]
   
    if return_list:
        gfa_s = []
        tmp = merge(*chunks, key=functools.cmp_to_key(compare_gfa2))
        for i in tmp:
            #print(i)
            if i[0] == "S":
                #print(i.rstrip())
                gfa_s.append(i.rstrip().split('\t'))
        
        for part_file in glob.glob(path):
            if os.path.isfile(part_file):
                os.remove(part_file)

        return gfa_s

    with open(out_path, 'w') as f_out:
        f_out.writelines(merge(*chunks, key=functools.cmp_to_key(compare_gfa2)))
    
    for part_file in glob.glob(path):
        if os.path.isfile(part_file):
            os.remove(part_file)


def compare_gfa(ln1, ln2):
    
    if not ln1[0] == "S" and not ln2[0] == "S":
        #If both are not S lines, then leave it
        return -1
    elif ln1[0] == "S" and not ln2[0] == "S":
        return -1
    elif not ln1[0] == "S" and ln2[0] == "S":
        return 1

    chr1 = [k for k in ln1 if k.startswith("SN:Z:")][0][5:]
    chr2 = [k for k in ln2 if k.startswith("SN:Z:")][0][5:]
    start1 = int([k for k in ln1 if k.startswith("SO:i:")][0][5:])
    start2 = int([k for k in ln2 if k.startswith("SO:i:")][0][5:])

    if chr1 == chr2:
        if start1 < start2:
            return -1
        else:
            return 1
    if chr1 < chr2:
        return -1
    else:
        return 1


def compare_gfa2(ln1, ln2):

    ln1 = ln1.rstrip().split('\t')
    ln2 = ln2.rstrip().split('\t')
    #print(ln1, ln2)
    if not ln1[0] == "S" and not ln2[0] == "S":
        #If both are not S lines, then leave it
        return -1
    elif ln1[0] == "S" and not ln2[0] == "S":
        return -1
    elif not ln1[0] == "S" and ln2[0] == "S":
        return 1

    chr1 = [k for k in ln1 if k.startswith("SN:Z:")][0][5:]
    chr2 = [k for k in ln2 if k.startswith("SN:Z:")][0][5:]
    start1 = int([k for k in ln1 if k.startswith("SO:i:")][0][5:])
    start2 = int([k for k in ln2 if k.startswith("SO:i:")][0][5:])

    if chr1 == chr2:
        if start1 < start2:
            return -1
        else:
            return 1
    if chr1 < chr2:
        return -1
    else:
        return 1
