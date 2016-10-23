'''
Created on Nov 13, 2014

@author: Igor Ruiz de los Mozos

Script transform bedgraph to bed.

Usage:
        python BEDgraph2BED.py input_file.bed output_file.bed

'''


import sys

def BEDgraph2BED(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    header = fin.readline()
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        count = float(col[3])
        
        if count < 0:
            strand = '-'
        else:
            strand = '+'
            
        fout.write(col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + str(abs(count)) + '\t' + "" + '\t' + strand + '\n')
        line = fin.readline()
                    
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    BEDgraph2BED(fname_in, fname_out)
else:
    print("python BEDgraph2BED.py <input_file> <output_file>")