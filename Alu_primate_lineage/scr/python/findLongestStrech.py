#!/usr/bin/env python
'''
Created on 1 March 2016

@author: Igor Ruiz de los Mozos

Find the longest strech of a given letter in a string (case insensitive)

Usage:

        python findLongestStrech.py input_fname.fasta output_fname.fasta LetterToFind
'''
import sys
import itertools
from Bio import SeqIO
from decimal import *


## Function to count the longest streeches of a given nucleotide
def count_streches(sequence, letter_to_find):
    lenght=0
    lenght=max(len(list(y)) for (c,y) in itertools.groupby(sequence) if c==letter_to_find)
    #print sequence
    return(lenght)


def longest_streche(fin, fout, strech_to_find):
    fin = open(fin, "rt")
    fout = open(fout, "w")
    strech_to_find=strech_to_find.upper()  ## 
    count=0
    sum=0

    ## Separate header from sequence
    for seq_record in SeqIO.parse(fin, "fasta"):
        #print(seq_record.id)
        #print(repr(seq_record.seq.upper()))
        #print seq_record.seq.upper()
        header = seq_record.id
        seq = seq_record.seq.upper()
        #print(len(seq_record))
        

        ## Check the are not NNNN
        if "N" in seq_record.seq.upper():
            #print "Found NN in the fasta sequecence \n"
            #print seq_record.seq
            continue
        elif strech_to_find not in seq_record.seq.upper():
            #print "Letter " + strech_to_find + " haven't been found on string " + seq_record.seq.upper()
            continue
        else:

            ## Count longest streech
            lenght=count_streches(seq_record.seq.upper(),strech_to_find)
            sum+=lenght
            count+=1
            ## Write to file
            fout.write(str(header) + '\t' + "" + '\t' + str(count) + '\t' + str(lenght)  + '\n')
            
            #print lenght
            #print sum
            #print count
            #print round(Decimal(sum)/Decimal(count), 3)
            
            if lenght > len(seq_record):
                print "Error lenght of streches cannot be longet that the string"
                quit()
    
    if count > 0:

        ## Print results on screen
        average=round(Decimal(sum)/Decimal(count), 3)
        print "Number of total sequences  " + str(count)            
        print "Average longest " + strech_to_find + "strech  " + str(average)





## Feed script
if sys.argv.__len__() == 4:
    fin_fasta = sys.argv[1]
    fout_fname = sys.argv[2]
    letter = str(sys.argv[3])
    longest_streche(fin_fasta, fout_fname, letter)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python findLongestStrech.py input_fname.bed output_fname.tab LetterToFind"


fin_fasta.close()
fout_fname.close()


'''
fin_fasta='/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/3SS_Alus/All_Aluexons_3SS_20_3_corrected.fasta'
fout_fname='/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/3SS_Alus/All_Aluexons_3SS_20_3_corrected_test_delete.tab'
letter="t"

'''
