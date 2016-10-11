'''
Created on Mar 16, 2014

@author: Nejc Haberman


Script will select trancript with the maximum exon number-
'''

import sys

def filter(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    header = fin.readline()
    fout.write(header)
    line = fin.readline()
    exon_num = 0
    max_CDS = 0 #longest CDS length
    max = None
    last_gene = None
    while line:
        col = line.rstrip('\n').rsplit('\t')
        #gene = col[1:5]
        gene = col[0]
        cds_start = int(col[5])
        cds_end = int(col[6])
        cds_len = cds_end - cds_start
        if last_gene == gene or last_gene == None:
            if cds_len > max_CDS:
                max = line
                max_CDS = cds_len
        else:
            fout.write(max)
            max_CDS = cds_len
            max = line
        
        last_gene = gene
        line = fin.readline()
    fout.write(max)
    fout.close()
    fin.close()

'''

filter("/media/skgthab/storage/UCL/2015.04.27@ALUs-Jan/CDSs/hg19-cds.bed","/media/skgthab/storage/UCL/2015.04.27@ALUs-Jan/CDSs/hg19-cds-max_CDS.bed")

'''
if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    filter(fname_in, fname_out)
else:
    print("You need two arguments to run the script")
