#!/usr/bin/env python
'''
Created on 3th February

@author: Igor Ruiz de los Mozos

The script will get the start position from a bed file. Used to get the 3'ss of alu exons.

example:
            $ python get_start_position_from_bed.py input_fname.bed output_fname.bed"

'''

import sys

def start_positions(fin_bed, fout_bed):
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
            if strand == '+':
                shift_start = int(start) -2
                shift_end = int(start) -1 
            else:
                shift_end = int(end) + 1
                shift_start = int(end) 
            if shift_start >= 0:
                fout.write(chr + '\t' + str(shift_start) + '\t' + str(shift_end) + '\t' + col[3] + '\t' + col[4] + '\t' + col[5] + '\n')
        line = fin.readline()



if sys.argv.__len__() == 3:
    fin_bed = sys.argv[1]
    fout_bed = sys.argv[2]
    start_positions(fin_bed, fout_bed)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python get_start_position_from_bed.py input_fname.bed output_fname.bed"
    