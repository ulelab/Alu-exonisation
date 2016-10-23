'''
Created on Mar 22, 2016

@author: Igor Ruiz de los Mozos

The script will get the distance between 3ss and Aluexon start (OR end of alu element). Return the 3SS single nt and the distance V$5

Feed file comes from bed intersect:
bedtools intersect -wao -a 'All_Aluexons.bed' -b '/media/igor/DATA/UCL/Evolution_Alus/Raw_Alus/rmsk_hg19_full_family_Alu_elements.bed' > './3SS_Alus/All_Aluexons_3SS_temp_distance_Alu.bed'
 The distance will be from the "start" of feature a to "start" feature b taking in acount the strand of them

Usage:
        $ python get_start_position_from_bed.py input_fname.bed output_fname.bed

'''

import sys



## Get the distance from the 3´ss to the end of Alu element
def get_distance_3SS_alu(fin_bed, fout_bed):
    fin = open(fin_bed, "rt")
    fout = open(fout_bed, "w")
    line = fin.readline()
    count_negatives = 0 
    while line:
        col = line.rstrip('\n').rsplit('\t')
        if line.__len__() > 10:
            chr_3ss = col[0]
            start_3ss = col[1]
            end_3ss = col[2]
            alu_type = col [3]
            strand_3ss = col[5]
            chr_Alu = col[6]
            start_Alu = col[7]
            end_Alu = col[8]
            strand_Alu = col[11]


            ## 3´ss on positive strand and Alu element on negative strand
            if strand_3ss == '+' and strand_Alu == '-':
                distance = int(start_3ss) - int(start_Alu)
                alu_distance =  int(end_Alu) - int(start_Alu)
                fout.write(chr_3ss + '\t' + str(start_3ss) + '\t' + str(end_3ss) + '\t' + col[3] + '\t' + col[4]  + '\t' + strand_3ss + '\t' + str(distance) + '\n')

            ## 3´ss on negative strand and Alu element on positive strand
            elif strand_3ss == '-' and strand_Alu == '+':
                distance = int(end_Alu) - int(end_3ss)
                alu_distance =  int(end_Alu) - int(start_Alu)
                fout.write(chr_3ss + '\t' + str(start_3ss) + '\t' + str(end_3ss) + '\t' + col[3] + '\t' + col[4]  + '\t' + strand_3ss + '\t' + str(distance) + '\n')
            
                
            else:
                print "Some bed lines are not suitable OR not in antisense of Alu " + chr_3ss + '\t' + str(start_3ss) + '\t' + str(end_3ss) + '\t' + col[3] + '\t' + chr_3ss + ":" + str(start_3ss) + ":" + str(end_3ss) + ":" + strand_3ss + '\t'  + chr_Alu + ":" + str(start_Alu) + ":" + str(end_Alu) + ":" + strand_Alu +  '\t' + str(distance) + '\t' + "" + '\t' + str(alu_distance) + '\n'
                                
            if distance < 0:
                print "ERROR NEGATIVE DISTANCE ##############################"
                count_negatives =  count_negatives + 1
                print count_negatives
        line = fin.readline()



## Argument pased to this script
if sys.argv.__len__() == 3:
    fin_bed = sys.argv[1]
    fout_bed = sys.argv[2]
    get_distance_3SS_alu(fin_bed, fout_bed)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python get_start_position_from_bed.py input_fname.bed output_fname.bed"
    
       
           
'''
test

fin_bed = '/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/3SS_Alus/All_Aluexons_3SS_AS_distance.bed'
fout_bed = '/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/3SS_Alus/All_Aluexons_3SS_AS_distance_to_end_temp_DELETE.bed'

'''
