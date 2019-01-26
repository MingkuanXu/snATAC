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
import threading
from time import ctime, sleep

# process one million reads at a time
CHUNK_SIZE = 1000000
# CHUNK_SIZE = 10000
coordinates = set()
MAX_THREAD = 120

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

def find_coordinates(chunk,peaks,regions,barcodes):
    # print("Handling the chunk at %s" % ctime())
    reads = pybedtools.BedTool(chunk)
    for line in peaks.intersect(reads, wa=True, wb=True):
        elems = str(line).split()
        cur_region = '\t'.join(elems[:3])
        cur_barcode = elems[6]
        if cur_barcode not in barcodes: continue
        if cur_region not in regions: sys.exit("error(main): region not in the list")
        coordinates.add((barcodes[cur_barcode], regions[cur_region]))

def create_barcode_hash(barcode_txt):
    """
    Create a hash table for the barcodes
    """
    barcodes = {}
    i = 0
    with open(barcode_txt) as fin:
        for line in fin:
            cur_barcode = line.strip().split()[0]
            barcodes[cur_barcode] = i
            i += 1
    return barcodes

def create_region_hash(peak_bed):
    """
    Create a hash table for the regions
    """
    regions = {}
    i = 0;
    with open(peak_bed) as fin:
        for line in fin:
            regions['\t'.join(line.split()[:3])] = i
            i += 1
    return regions

def update_matrix_with_coordinates(coordinates,i,j):
    mat = np.zeros((j, i))
    for each in coordinates:
        mat[each[0], each[1]] = 1
    return mat

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
    output =  "result.txt"
    #options.output

    # read peaks and reads
    peaks = pybedtools.BedTool(peak_bed)
    reads = pybedtools.BedTool(read_bed)

    # create hash tables for peaks and barcodes
    barcodes = create_barcode_hash(barcode_txt)
    regions = create_region_hash(peak_bed)

    fin = open_file(read_bed)

    while True:
        # chunk1 = next_n_lines(fin, CHUNK_SIZE)
        chunk1 = next_n_lines(fin, CHUNK_SIZE)
        if(len(chunk1) == 0): break
        # find_coordinates(chunk1,peaks,regions,barcodes)
        t1 = threading.Thread(target = find_coordinates,args=(chunk1,peaks,regions,barcodes,))
        t1.start()
        if(threading.activeCount()>MAX_THREAD):
            t1.join()
        # print("Current threading is %d" % threading.activeCount())

    fin.close()
    t1.join()

    matrix = update_matrix_with_coordinates(coordinates,len(regions),len(barcodes))

    np.savetxt(output, matrix, delimiter='\t', fmt="%d")

if __name__ == '__main__':
    main()
