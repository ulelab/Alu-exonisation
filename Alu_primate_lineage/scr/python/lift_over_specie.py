#!/usr/bin/env python

'''
Created on 2 March 2016

@author: Igor Ruiz de los Mozos

Lift over bed file to a new bed file specifying the specie conversion. Optional flag to get fasta from those liftovers

Specify genome_dir variable to download chains for lift over and genome fasta to get fasta

Usage:  python lift_over_specie.py bedIN lift_from lift_to bedOUT optionalFlag  \n optionalFlag_getFasta=TRUE, FALSE OR empty

Available genomes:
#Human, Chimp, Bonobo, Gibon, Baboon, Rhesus, Marmoset, Tarsier, Lemur, Bushbaby, Tree shrew
#hg38, panTro4, panPan1, nomLeu1, papHam1, rheMac3, calJac3, tarSyr2, micMur1, otoGar1, tupBel1

'''
import sys
import os
import subprocess


## Function allow to run bash programs in python
def runUNIXCommands(commands):
    print commands
    try:
        retcode = subprocess.check_call(commands, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
        else:
            print >>sys.stderr, "Child returned", retcode
    except OSError, e:
        print >>sys.stderr, "Execution failed:", e


## Arguments passed when calling this script
if sys.argv.__len__() >= 5:
    fin_fname = sys.argv[1]
    lift_from = str(sys.argv[2])
    lift_to = str(sys.argv[3])
    fout_name = sys.argv[4]
    optionalFlag = sys.argv[5]
    
else:
    #print str(sys.argv.__len__())
    print "error:\t4 \n Usage: \t python lift_over_specie.py bedIN lift_from lift_to bedOUT optionalFlag  \n optionalFlag_getFasta=TRUE, FALSE OR empty"  


## Asign lift_from and chain_url variable - This depend on starting genome final genome to lift over -



if lift_from=='hg19' and lift_to=='hg38':
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
    genome_url='hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
elif lift_from=='hg38' and lift_to=='panTro4':                                                          ## Chimp
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg38/vsPanTro4/hg38.panTro4.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/panTro4/bigZips/panTro4.fa.gz'
elif lift_from=='hg38' and lift_to=='panPan1':                                                          ## Bonobo
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg38/vsPanPan1/hg38.panPan1.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/panPan1/bigZips/panPan1.fa.gz'
elif lift_from=='hg38' and lift_to=='rheMac3':                                                          ## Rhesus
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg38/vsRheMac3/hg38.rheMac3.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/rheMac3/bigZips/rheMac3.fa.gz'
elif lift_from=='hg38' and lift_to=='tarSyr2':                                                          ## Tarsier
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg38/vsTarSyr2/hg38.tarSyr2.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/tarSyr2/bigZips/tarSyr2.fa.gz'
elif lift_from=='hg19' and lift_to=='otoGar1':                                                          ## hg19 to Bushbaby
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg19/vsOtoGar1/hg19.otoGar1.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/otoGar1/bigZips/otoGar1.fa.gz'
elif lift_from=='hg19' and lift_to=='micMur1':                                                          ## hg19 to Mouse Lemur
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg19/vsMicMur1/hg19.micMur1.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/micMur1/bigZips/micMur1.fa.gz'
elif lift_from=='hg19' and lift_to=='tupBel1':                                                          ## hg19 to Tree shrew
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg19/vsTupBel1/hg19.tupBel1.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/tupBel1/bigZips/tupBel1.fa.gz'
elif lift_from=='hg19' and lift_to=='speTri1':                                                          ## hg19 to Squirrel   NO GENOME FOUND!!!
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg19/vsSpeTri1/hg19.speTri1.all.chain.gz'
    genome_url=''
elif lift_from=='hg19' and lift_to=='nomLeu1':                                                          ## hg19 to Gibbon
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg19/vsNomLeu1/hg19.nomLeu1.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/nomLeu1/bigZips/nomLeu1.fa.gz'
elif lift_from=='hg19' and lift_to=='papHam1':                                                          ## hg19 to Baboon
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg19/vsPapHam1/hg19.papHam1.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/papHam1/bigZips/papHam1.fa.gz'
elif lift_from=='hg19' and lift_to=='calJac3':                                                          ## hg19 to Marmoset
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/hg19/vsCalJac3/hg19.calJac3.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/calJac3/bigZips/calJac3.fa.gz'
elif lift_from=='tarSyr2' and lift_to=='otoGar3':                                                       ## Tarsier to Bushbaby   Give errors
    chain_url='hgdownload.cse.ucsc.edu/goldenPath/tarSyr2/vsOtoGar3/tarSyr2.otoGar3.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/otoGar3/bigZips/otoGar3.fa.gz'
elif lift_from=='hg19' and lift_to=='mm10':                                                             ## Tarsier to Bushbaby
    chain_url='hgdownload.soe.ucsc.edu/goldenPath/hg19/vsMm10/hg19.mm10.all.chain.gz'
    genome_url=''
elif lift_from=='mm10' and lift_to=='otoGar3':                                                          ## Tarsier to Bushbaby
    chain_url='hgdownload.soe.ucsc.edu/goldenPath/mm10/vsOtoGar3/mm10.otoGar3.all.chain.gz'
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/otoGar3/bigZips/otoGar3.fa.gz'     
else:
    print "Unknown species to liftOver \n try again"
    quit()

genome_dir='/media/igor/DATA/UCL/Evolution_Alus/Genomes/'
genome_dir= genome_dir + lift_to + '/'
print chain_url + "\n" + genome_dir + "\n" + genome_url + "\n"



runUNIXCommands("mkdir -p " + genome_dir)


def liftOver_species(bed_in, fromm, to, bed_out, optionalFlag):
    
    ## Download the chain to liftOver
    runUNIXCommands("rsync -avzP rsync://" + chain_url + " " + genome_dir)        ## Download the chain to convert
    
    ## Asign names
    temp_chain_name= genome_dir + chain_url.split("/")[-1].replace(".gz", "") 
    bed_in = bed_in.replace(".bed", "")
    bed_out = bed_out.replace(".bed", "")
    
    ## Decompress
    runUNIXCommands("gunzip " + temp_chain_name + ".gz")
    
    ## LiftOver asigned genomes
    liftoverCommand = "/home/igor/Programs/bin/liftOver " + bed_in + ".bed " + temp_chain_name  + " " + bed_out + ".bed " +  bed_out + "_unmaped.bed"
    runUNIXCommands(liftoverCommand)
    runUNIXCommands("rm " + bed_out+ "_unmaped.bed")
    
    runUNIXCommands("rm " + temp_chain_name)                                           ## Remove or compress chain not both
    #runUNIXCommands("gzip " + temp_chain_name)
    
    ## get fasta from those beds if optionalFlag = TRUE is specified
    
    if optionalFlag=="TRUE":
        print " Do optionalFlag get fasta from liftOvered beds"
        runUNIXCommands("rsync -avzP rsync://" + genome_url + " " + genome_dir)        ## Download the fasta sequence from this specie
        temp_genome_name= genome_dir + genome_url.split("/")[-1].replace(".gz", "")
        runUNIXCommands("gunzip " + temp_genome_name + ".gz")
        runUNIXCommands("bedtools getfasta -s -fi "+ temp_genome_name + " -bed " + bed_out + ".bed" + " -fo " + bed_out + ".fasta")
        #runUNIXCommands("rm " + temp_genome_name)                                           ## Remove or compress genome not both
        runUNIXCommands("gzip " + temp_genome_name)
        
    elif optionalFlag=="FALSE" or optionalFlag=="":
        print " No OptionalFlag get fasta"    
    else:
        print " No OptionalFlag get fasta"
        
liftOver_species(fin_fname, lift_from, lift_to, fout_name, optionalFlag)        



''' Tests
if sys.argv.__len__() == 5:
    fin_fname = sys.argv[1]
    lift_from = str(sys.argv[2])
    lift_to = str(sys.argv[3])
    fout_fname = sys.argv[4]
    liftOver_species(fin_fname, lift_from, lift_to, fout_name)
else:
    #print str(sys.argv.__len__())
    print "error:\t4 \n Usage: \t python lift_over_specie.py bedIN lift_from lift_to bedOUT"       
    
    fin_fname='/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/3SS_Alus/All_Aluexons_3SS_20_3_corrected.bed'
    fout_name='/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/3SS_Alus/All_Aluexons_3SS_20_3_corrected_transformed_delete.bed'
    lift_from='hg19'
    lift_to='otoGar1'
    optionalFlag="TRUE"
'''