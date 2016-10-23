'''
Created on Aug 20, 2013

@author: Igor Ruiz de los Mozos

The script will get the 3'SS 20 nt inside alu element 
'''
import sys

def get3SS(fin_fname, fout_fname, distance_inside_alu):
    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    dist = distance_inside_alu
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        if line.__len__() > 6:
            chr_alu = col[0]
            start_alu = col[1]
            end_alu = col[2]
            strand_alu = col[5]
            shift_start = shift_end = 0
            if strand_alu == '+':
                shift_start = int(start_alu) + int(dist) - 1
                shift_end = shift_start + 1
                #fout.write(chr_alu + '\t' + str(shift_start) + '\t' + str(shift_end) + '\t' + col[3]  + '\t' + "1" + '\t' + strand_alu + '\t' + col[0] + ':' + col[1] + ':' + col[2] + ':' + strand_alu + '\n')
                fout.write(chr_alu + '\t' + str(shift_start) + '\t' + str(shift_end) + '\t' + col[3]  + '\t' + "1" + '\t' + strand_alu  + '\n')
     
            else:
                shift_end = int(end_alu) - int(dist) 
                shift_start = shift_end -1
                #fout.write(chr_alu + '\t' + str(shift_start) + '\t' + str(shift_end) + '\t' + col[3] + '\t' + "-1" + '\t' + strand_alu + '\t' + col[0] + ':' + col[1] + ':' + col[2] + ':' + strand_alu + '\n')
                fout.write(chr_alu + '\t' + str(shift_start) + '\t' + str(shift_end) + '\t' + col[3] + '\t' + "-1" + '\t' + strand_alu+ '\n')

                     
        line = fin.readline()

if sys.argv.__len__() == 4:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    distance_inside_alu = sys.argv[3]
    get3SS(fin_fname, fout_fname, distance_inside_alu)
else:
    #print str(sys.argv.__len__())
    print "error:\t4 arguments are needed\n" + '\n' +"example:\t $ pythonget3SS_from_random_alu.py input_fname.bed output_fname.bed"