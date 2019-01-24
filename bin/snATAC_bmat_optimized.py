#!/usr/bin/env python

"""
create binary accessible matrix
Created by Rongxin Fang
"""

import sys
import numpy as np
from operator import itemgetter
import gzip
import bz2
from itertools import islice

# process one million reads at a time
CHUNK_SIZE=1000000

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


def find_regions(peak_bed):
    regions = []
    with open(peak_bed) as fin:
        for line in fin:
            regions.append('\t'.join(line.split()[:3]))
    return regions

def find_barcode(barcode_txt):
    barcodes = []
    with open(barcode_txt) as fin:
        for line in fin:
            barcodes.append(line.strip().split()[0])
    return barcodes


def main():
    '''
    from argparse import ArgumentParser
    # parameters
    parser = ArgumentParser(description='generate binary accessibility matrix')
    parser.add_argument('-i', '--input', help='bed file contains read', required=True)
    parser.add_argument('-x', '--barcode', help='rows: file contains selected cell barcode', required=True)
    parser.add_argument('-y', '--peak', help='columns: file contains selected peaks', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    options = parser.parse_args()
    '''
    # input parsing
    '''
    read_bed = options.input
    peak_bed = options.peak
    barcode_txt = options.barcode
    output = options.output
    '''

    read_bed = "../data/p56.rep1.bed_100000"
    # options.input
    peak_bed = "../data/p56.rep1.ygi"
    # options.peak
    barcode_txt = "../data/p56.rep1.xgi"

    regions = find_regions(peak_bed)
    barcodes = find_barcode(barcode_txt)

    # Initialize matrix based on the number of peaks & barcodes
    matrix = np.zeros((len(barcodes), len(regions)))

    f = open_file(read_bed)
    line = f.readline()
    while line:
        chr,start_index,end_index,barcode = line.strip().split()
        print(chr)
        print(start_index,end_index,barcode)
        line = f.readline()

    while True:
        chunk = next_n_lines(fin, CHUNK_SIZE)
        if(len(chunk) == 0): break
        reads = pybedtools.BedTool(chunk)
        for line in peaks.intersect(reads, wa=True, wb=True):
            elems = str(line).split()
            cur_region = '\t'.join(elems[:3])
            cur_barcode = elems[6]
            if cur_barcode not in barcodes: continue
            if cur_region not in regions: sys.exit("error(main): region not in the list")
            mat[barcodes[cur_barcode], regions[cur_region]] += 1
    fin.close()

    # convert the matrix to a binary matrix
    mat[ np.where( mat > 1 ) ] = 1

    #np.savetxt(output, mat, delimiter='\t', fmt="%d")

if __name__ == '__main__':
    main()
