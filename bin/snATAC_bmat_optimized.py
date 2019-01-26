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
    """
    Put each line of peak information into a list of tuples with the form of
    [('chr','start_index','end_index'), ...] and return the list.
    """
    regions = []
    with open(peak_bed) as fin:
        for line in fin:
            regions.append(line.split()[:3])
    return regions

def build_chr_dic(regions):
    """
    Build a dictionary with chromosomes as its key and a list of tuple of start & end
    indexes of a peak as its item. E.g.
    chr_dic = {'chr1':[(0,80),(99,184), ...]
               'chr2':[...]
               ...}
    """
    chr_dic = {}
    for each in regions:
        if each[0] not in chr_dic:
            chr_dic[each[0]] = []
        chr_dic[each[0]].append((each[1],each[2]))
    return chr_dic

def find_barcode(barcode_txt):
    """
    Put each unique barcode into a list and return it.
    """
    barcodes = []
    with open(barcode_txt) as fin:
        for line in fin:
            barcodes.append(line.strip())
    return barcodes

def find_region_index(start_index,end_index,chr,chr_dic):
    """
    For a read, determine if its start and end index is in the range of a peak.
    If so, return the index of the peak it is in. Else, return -1.
    """
    if chr not in chr_dic:
        return -1
    all_regions = chr[dic]
    for i in range(all_regions):
        region_start = all_regions[i][0]
        region_end = all_regions[i][1]

def next_n_lines(file_opened, N):
    """
    Read N lines at one time.
    """
    return [x for x in islice(file_opened, N)]


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
    chr_dic = build_chr_dic(regions)

    # Initialize matrix based on the number of peaks & barcodes
    matrix = np.zeros((len(barcodes), len(regions)))

    f = open_file(read_bed)
    line = f.readline()
    while line:
        chr,start_index,end_index,barcode = line.strip().split()
        if barcode not in barcodes:
            # The current barcode is not found in the barcode list
            line = f.readline()
            continue

        barcode_index = barcodes.find(barcode)
        region_index = find_region_index(start_index,end_index,chr,chr_dic)
        if region_index != -1:
            matrix[barcode_index,region_index] +=1
        line = f.readline()

    f.close()

    # convert the matrix to a binary matrix
    mat[ np.where( mat > 1 ) ] = 1

    #np.savetxt(output, mat, delimiter='\t', fmt="%d")
    total = 0
    for each in mat:
        total += sum(each)
    print(total) # should be 5952

if __name__ == '__main__':
    main()
