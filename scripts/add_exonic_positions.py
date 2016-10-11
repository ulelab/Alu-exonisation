'''
Created on Mar 16, 2014

@author: Nejc Haberman


Script will use intersect data from: bedtools intersect -s -a regulated.exons.bed -b hg19-cds-max_exons2.bed -wb

Results will be CDS genomic positions with infromation of exons on transcript level together with Alu exons 
'''

import sys

def get_positions_plus_strand(cds_start,exon_starts,exon_ends):
    positions = []  #exonic positions inside the transcript positions
    for i in range(0,exon_starts.__len__()):
        start = int(exon_starts[i]) - cds_start + 1 #we need to adjust position to UCSC +1 counting
        end = int(exon_ends[i]) - cds_start 
        if end >= 0:  #we ignore exons that are before coding (CDS) sequence
            position = str(start) + ':' + str(end)
            positions.append(position)
    return positions
            
        
def get_positions_minus_strand(cds_end,exon_starts,exon_ends):
    positions = []  #exonic positions inside the transcript positions
    for i in range(0,exon_starts.__len__()):
        start = (cds_end - int(exon_ends[i])) + 1 #we need to adjust position to UCSC +1 counting 
        end = (cds_end - int(exon_starts[i]))
        if end >= 0:  #we ignore exons that are before coding (CDS) sequence
            position = str(start) + ':' + str(end)
            positions.append(position)
    return positions[::-1]
        
def filter(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        alu_exon_start = int(col[1])
        alu_exon_end = int(col[2])
        alu_exon_id = col[3]
        strand = col[5]
        cds_start = int(col[9])
        cds_end = int(col[10])
        exon_starts = col[11].rstrip(',').rsplit(',')
        exon_ends = col[12].rstrip(',').rsplit(',')

        if strand == '+':   #calculate positions on CDS level
            cds_alu_exon_start = alu_exon_start - cds_start
            cds_alu_exon_end = alu_exon_end - cds_start
            cds_exons = get_positions_plus_strand(cds_start,exon_starts,exon_ends)
        elif strand == '-':
            alu_exon_start -= 1
            alu_exon_end -= 1
            cds_alu_exon_start = cds_end - alu_exon_end
            cds_alu_exon_end = cds_end - alu_exon_start
            cds_exons = get_positions_minus_strand(cds_end,exon_starts,exon_ends)
            
        info = alu_exon_id + '|' + str(cds_alu_exon_start) + ':' + str(cds_alu_exon_end) + '|' + str(cds_exons).replace('[','').replace(']','').replace("'","").replace(" ","") + '|' + strand
        fout.write(chr + '\t' + str(cds_start) + '\t' + str(cds_end) + '\t' + info + '\t' + "" + '\t' + strand + '\n')
            
        line = fin.readline()
    fout.close()
    fin.close()

'''
filter("/media/skgthab/storage/UCL/2015.04.27@ALUs-Jan/CDSs3/hg19-cds-gene_symbols-longest-CDS2-regulated_exons.bed","/media/skgthab/storage/UCL/2015.04.27@ALUs-Jan/CDSs3/temp.bed")

'''
if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    filter(fname_in, fname_out)
else:
    print("You need two arguments to run the script")
