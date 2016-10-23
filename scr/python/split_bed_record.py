#!/usr/bin/env python
'''
Created on 5 April 2016

@author: Igor Ruiz de los Mozos

Script feed with a multiple bed file
return and split each line in separate file with file name equal to string repesenting the bed position on the genome

Usage:

        python split_bed_record.py bedIN.bed DirOUT"


'''
import subprocess
import sys

## Run bash comand on python
def runUNIXCommands(commands):
    #print commands
    try:
        retcode = subprocess.check_call(commands, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
        else:
            #print >>sys.stderr, "Child returned", retcode
            pass
    except OSError, e:
        print >>sys.stderr, "Execution failed:", e


## Split each line and save to file with bed name
def split_bed(fin_fname, genome_dir):
    fin = open(fin_fname, "rt")
    line = fin.readline()         # Initialize variables to count each position


    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        start = col[1]
        end = col[2]
        strand = col[5]
        runUNIXCommands("mkdir -p " + genome_dir)
        
        #print "BED line   " + line
        
        ## Asign name equeal to bed positions
        temp_in_v1 = (genome_dir + chr + ":" + start + ":" + end  + ".bed")
        
        out = open(temp_in_v1, "w")
        out.write(line)
        out.close()
        
        line = fin.readline()   
        
        
        
#split_bed(fin_fname, genome_dir)

if sys.argv.__len__() == 3:
    fin_fname = sys.argv[1]
    genome_dir = sys.argv[2]
    split_bed(fin_fname, genome_dir)
else:
    #print str(sys.argv.__len__())
    print "error:\t4 \n Usage: \t python split_bed_record.py bedIN DirOUT"

'''
fin_fname='/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/3SS_Alus/All_Aluexons_3SS_C.bed'
genome_dir='/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/3SS_Alus/test2/'


'''