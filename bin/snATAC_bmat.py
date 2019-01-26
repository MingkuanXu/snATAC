#!/usr/bin/env python

"""
create binary accessible matrix
Created by Rongxin Fang
"""

import sys
import numpy as np
from operator import itemgetter
import pybedtools
import gzip
import bz2
from itertools import islice

# process one million reads at a time
CHUNK_SIZE=1000000
mat = np.zeros((0, 0))

magic_dict = {
    "\x1f\x8b\x08": "gz",
    "\x42\x5a\x68": "bz2",
    "\x50\x4b\x03\x04": "zip"
    }

max_len = max(len(x) for x in magic_dict)

def file_type(filename):
    with open(filename) as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return "txt"

def open_file(fname):
    if file_type(fname) == "gz":
        fin = gzip.open(fname, 'rb')
    elif file_type(fname) == "bz2":
        fin = bz2.BZ2File(fname, 'r')
    elif file_type(fname) == "txt":
        fin = open(fname, 'r')
    return fin

def next_n_lines(file_opened, N):
    """
    Read N lines at one time.
    """
    return [x for x in islice(file_opened, N)]

def chunkIt(seq, num):
    """
    seperate the list into every two pair
    i.g. [1,2,3,4,5,6,7,8] -> [[1,2], [3,4], [5,6], [7,8]]
    every two read as a file
    """
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

def main():
    from argparse import ArgumentParser
    # parameters
    # parser = ArgumentParser(description='generate binary accessibility matrix')
    # parser.add_argument('-i', '--input', help='bed file contains read', required=True)
    # parser.add_argument('-x', '--barcode', help='rows: file contains selected cell barcode', required=True)
    # parser.add_argument('-y', '--peak', help='columns: file contains selected peaks', required=True)
    # parser.add_argument('-o', '--output', help='output file', required=True)
    # options = parser.parse_args()
    # input parsing
    read_bed = "../data/p56.rep1.bed_100000"
    # options.input
    peak_bed = "../data/p56.rep1.ygi"
    # options.peak
    barcode_txt = "../data/p56.rep1.xgi"
    # options.barcode
    output =  "result1.txt"
    #options.output
    # read peaks and reads
    peaks = pybedtools.BedTool(peak_bed)
    reads = pybedtools.BedTool(read_bed)

    # find overlap
    # ov = peaks.intersect(reads, wa=True,  wb=True)

    regions = {}
    barcodes = {}

    i = 0; j = 0;
    with open(peak_bed) as fin:
        for line in fin:
            regions['\t'.join(line.split()[:3])] = i
            i += 1

    with open(barcode_txt) as fin:
        for line in fin:
            cur_barcode = line.strip().split()[0]
            barcodes[cur_barcode] = j
            j += 1

    global mat
    mat = np.zeros((j, i))

    peaks = pybedtools.BedTool(peak_bed)


    counter = 0
    fin = open_file(read_bed)
    while True:
        chunk = next_n_lines(fin, CHUNK_SIZE)
        if(len(chunk) == 0): break
        update_matrix_with_chunk(chunk,regions,barcodes,peaks)
    fin.close()

    # convert the matrix to a binary matrix
    # mat[ np.where( mat > 1 ) ] = 1

    np.savetxt(output, mat, delimiter='\t', fmt="%d")



def update_matrix_with_chunk(chunk,regions,barcodes,peaks):
    reads = pybedtools.BedTool(chunk)
    for line in peaks.intersect(reads, wa=True, wb=True):
        elems = str(line).split()
        cur_region = '\t'.join(elems[:3])
        cur_barcode = elems[6]
        if cur_barcode not in barcodes: continue
        if cur_region not in regions: sys.exit("error(main): region not in the list")
        mat[barcodes[cur_barcode], regions[cur_region]] = 1


if __name__ == '__main__':
    main()
