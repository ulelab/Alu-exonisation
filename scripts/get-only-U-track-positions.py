'''
Created on Apr 25, 2014
 
@author: Nejc Haberman
 
 Scrit will receive a fasta format and return positions of U tracks that are longer then 4 nts and smaller then 100 nt
 
input:
- input_fasta
- output_fasta
 
output:
- positions
'''
import sys
 

#function will convert long uridine tracks (more then 4ts) to number 4 in the sequence and it will also return positions (left, middle, right) on the ALU position they are
def get_uridine_track_positions(seq):
    max = 100    #maximum number of T repeats
    min = 4 #minimum number of repeats
    positions = []
    for i in reversed(range(min,max)):
        y_track = "T" * i
        replacement = "4" * i
        position = seq.find(y_track)
        while position != -1:
            positions.append(position)
            seq = seq.replace(y_track, replacement,1)   #replace only the first one
            position = seq.find(y_track)
    return positions
    

def get_U_track_positions(fin_fname, fout_fname):
    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    fout.write("info\tpositions\n")
    line = fin.readline()
    positions = []
    while line:
        if line[0] == '>':
            info = line.rstrip('\n').replace('>','')
        else:
            seq = line.rstrip('\n')
            seq = str(seq).upper()
            positions = get_uridine_track_positions(seq)
            positions.sort()
            positions.reverse()
            fout.write(info + '\t' + str(positions).replace('[','').replace(']','') + '\n')
        line = fin.readline()
    
    fin.close()
    fout.close()
               
'''
fin_fname_fasta = "/media/skgthab/storage/UCL/2015.04.27@ALUs-Jan/HeatMap/HeatMaps8-2/U-tracks/test/test.fasta"
fout_fname_fasta = "/media/skgthab/storage/UCL/2015.04.27@ALUs-Jan/HeatMap/HeatMaps8-2/U-tracks/test/test.csv"
get_U_track_positions(fin_fname_fasta, fout_fname_fasta)
'''
 
if sys.argv.__len__() == 3:
    fin_fname_fasta = sys.argv[1]
    fout_fname_fasta = sys.argv[2]
    get_U_track_positions(fin_fname_fasta, fout_fname_fasta)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python k-mer_coverage.py input_fname.fasta motifs.tab"

