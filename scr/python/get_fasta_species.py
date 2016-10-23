#!/usr/bin/env python
'''
Created on 2 March 2016

@author: Igor Ruiz de los Mozos

Get fasta sequence from a specified genome.
Specify genome_dir variable to download chains for lift over and genome fasta to get fasta

Usage:

python get_fasta_species.py bedIN.bed specieVERSION fastaOUT.fasta

Species availables: 

hg38, hg19, panTro4, panPan1, rheMac3, tarSyr2, otoGar1, micMur1, tupBel1, nomLeu1, papHam1, calJac3, otoGar3, mm10



'''
import sys
import os
import subprocess


## Function allow to run system programs on python scripts

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


## Input arguments

if sys.argv.__len__() == 4:
    fin_fname = sys.argv[1]
    fasta_from = str(sys.argv[2])
    fout_name = sys.argv[3]
   
else:
    print str(sys.argv.__len__())
    print "error:\t4 \n Usage: \t python get_fasta_species.py bedIN specieGENOME fastaOUT \n Species availables hg38, hg19, panTro4, panPan1, rheMac3, tarSyr2, otoGar1, micMur1, tupBel1, nomLeu1, papHam1, calJac3, otoGar3, mm10 \n Example:  \n python get_fasta_species.py gene_position_to_get_fasta.bed hg19 output.fasta \n\n"


## Asign URL to download genome in function of genome specie/version

if fasta_from=='hg38':
    genome_url='hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
elif fasta_from=='hg19':                                                             ## hg19
    #genome_url='hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz' 
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz'
elif fasta_from=='panTro4':                                                          ## Chimp
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/panTro4/bigZips/panTro4.fa.gz'
elif fasta_from=='panPan1':                                                          ## Bonobo
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/panPan1/bigZips/panPan1.fa.gz'
elif fasta_from=='rheMac3':                                                          ## Rhesus
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/rheMac3/bigZips/rheMac3.fa.gz'
elif fasta_from=='tarSyr2':                                                          ## Tarsier
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/tarSyr2/bigZips/tarSyr2.fa.gz'
elif fasta_from=='otoGar1':                                                          ## Bushbaby1
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/otoGar1/bigZips/otoGar1.fa.gz'
elif fasta_from=='micMur1':                                                          ## Mouse Lemur
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/micMur1/bigZips/micMur1.fa.gz'
elif fasta_from=='tupBel1':                                                          ## Tree shrew
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/tupBel1/bigZips/tupBel1.fa.gz'
elif fasta_from=='speTri1':                                                          ## Squirrel   NO GENOME FOUND!!!
    genome_url=''
elif fasta_from=='nomLeu1':                                                          ## Gibbon
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/nomLeu1/bigZips/nomLeu1.fa.gz'
elif fasta_from=='papHam1':                                                          ## Baboon
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/papHam1/bigZips/papHam1.fa.gz'
elif fasta_from=='calJac3':                                                          ## Marmoset
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/calJac3/bigZips/calJac3.fa.gz'
elif fasta_from=='otoGar3':                                                          ## Bushbaby
    genome_url='hgdownload.cse.ucsc.edu/goldenPath/otoGar3/bigZips/otoGar3.fa.gz'
elif fasta_from=='mm10':                                                             ## Mouse
    genome_url=''

    
else:
    print "Unknown specie to get fasta \n try again \n\n"
    print "Species availables: \n hg38, hg19, panTro4, panPan1, rheMac3, tarSyr2, otoGar1, micMur1, tupBel1, nomLeu1, papHam1, calJac3, otoGar3, mm10"
    quit()

genome_dir='./Data/Genomes/'
genome_dir= genome_dir + fasta_from + '/'
print "Downloading genome *************** \n" + "genome_dir + "\n" + genome_url + "\n"


## Create folder with genome name
runUNIXCommands("mkdir -p " + genome_dir)


def get_fasta_species(bed_in, fromm, fasta_out):
    
    runUNIXCommands("rsync -avzP rsync://" + genome_url + " " + genome_dir)        ## Download the fasta sequence from this specie
    
    temp_genome_name= genome_dir + genome_url.split("/")[-1].replace(".gz", "")     ## Get genome name
    runUNIXCommands("gunzip " + temp_genome_name + ".gz")

    ## using bedtools to get the fasta seq
    runUNIXCommands("bedtools getfasta -s -fi "+ temp_genome_name + " -bed " + bed_in + " -fo " + fasta_out )
    runUNIXCommands("rm " + temp_genome_name)                                           ## Remove or compress genome not both
    #runUNIXCommands("gzip " + temp_genome_name)
        

        
get_fasta_species(fin_fname, fasta_from, fout_name)        

