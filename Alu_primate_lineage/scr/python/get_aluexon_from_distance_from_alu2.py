#!/usr/bin/env python

'''
Created on Mar 14, 2016

@author: Igor Ruiz de los Mozos

The script will get the distance between 3ss and Alu start.
In this bedtools file the distance from the intersection to the 3´s comes on the 7th column
The distance will be from the "start" of feature a to "start" feature b taking in acount the strand of them

Feed file comes from bed intersect:
bedtools intersect -wao -a 'All_Aluexons.bed' -b '/media/igor/DATA/UCL/Evolution_Alus/Raw_Alus/rmsk_hg19_full_family_Alu_elements.bed' > './3SS_Alus/All_Aluexons_3SS_temp_distance_Alu.bed'


Returns a bed file with a syntetic Alu element 320 nt long and place the 3´ss on it


Usage:
            python get_aluexon_from_distance_from_alu2.py input_fname.bed output_fname.bed"


'''

import sys





def get_aluexon(fin_bed, fout_bed):
    fin = open(fin_bed, "rt")
    fout = open(fout_bed, "w")
    line = fin.readline()
    count_negatives = 0

    while line:
        col = line.rstrip('\n').rsplit('\t')
        if line.__len__() > 6:
            chr_3ss = col[0]
            start_3ss = col[1]
            end_3ss = col[2]
            alu_type = col[3]
            distance = col[6]
            strand_3ss = col[5]


            if strand_3ss == '+':
                shift_start = int(start_3ss) - int(distance) - 20
                shift_end = int(start_3ss) - int(distance) + 300
                if distance > 0 and shift_start > 0:
                    fout.write(chr_3ss + '\t' + str(shift_start) + '\t' + str(shift_end) + '\t' + col[3] + '\t' + chr_3ss + ":" + str(start_3ss) + ":" + str(end_3ss) + '\t' + strand_3ss + '\n')
                
     
     
            elif strand_3ss == '-':
                shift_start = int(start_3ss) + int(distance) - 300 
                shift_end = int(start_3ss) + int(distance) + 20                               
                if distance > 0 and shift_start > 0:
                    fout.write(chr_3ss + '\t' + str(shift_start) + '\t' + str(shift_end) + '\t' + col[3] + '\t' + chr_3ss + ":" + str(start_3ss) + ":" + str(end_3ss) + '\t' + strand_3ss + '\n')
                
            else:
                print "Some bed lines are not suitable " + chr_3ss + '\t' + str(shift_start) + '\t' + str(shift_end) + '\t' + col[3] + '\t' + chr_3ss + ":" + str(start_3ss) + ":" + str(end_3ss) + ":" + strand_3ss + '\t'  + chr_Alu + ":" + str(start_Alu) + ":" + str(end_Alu) + ":" + strand_Alu +  '\t' + str(distance) + '\t' + str(distance2) + '\t' + str(alu_distance) + '\n'
                                
               
        line = fin.readline()



if sys.argv.__len__() == 3:
    fin_bed = sys.argv[1]
    fout_bed = sys.argv[2]
    get_aluexon(fin_bed, fout_bed)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python get_aluexon_from_distance_from_alu.py input_fname.bed output_fname.bed"
    
     


'''

fin_bed = '/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/All_Aluexons_3SS_AS_distance_to_alu/All_Aluexons_3SS_AS_distance_to_alu_hg19-hg38.bed'
fout_bed = '/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/All_Aluexons_3SS_AS_distance_to_alu/All_Aluexons_3SS_AS_distance_to_alu_hg19-hg38_DELETE.bed'

'''
            
