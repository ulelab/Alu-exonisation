#!/usr/bin/env python
'''
Created on April 20 2016

@author: Igor Ruiz de los Mozos

ADD ID to column 4th. Id is formed by the joining of Bed positions

Usage:
        example: $ python get_start_position_from_bed.py input_fname.bed output_fname.bed
'''

import sys

def add_ID(fin_bed, fout_bed):
    fin = open(fin_bed, "rt")
    fout = open(fout_bed, "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        if line.__len__() > 6:
            chr = col[0]
            start = col[1]
            end = col[2]
            strand = col[5]

            fout.write(chr + '\t' + start + '\t' + end + '\t' + col[3] + '\t' + chr + ":" + start + ":" + end + ":" + strand + '\t' + strand + '\n')
        line = fin.readline()



if sys.argv.__len__() == 3:
    fin_bed = sys.argv[1]
    fout_bed = sys.argv[2]
    add_ID(fin_bed, fout_bed)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python get_start_position_from_bed.py input_fname.bed output_fname.bed"
    