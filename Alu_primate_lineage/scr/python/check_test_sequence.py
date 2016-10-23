'''
Created on Mar 9, 2016

@author: Igor

The script will check that the fasta sequence is apropiate for downstream analysis

Check that the sequence only contains ATGC nucleotides
Return Uper case sequences

Usage:

        python check_test_sequence fastaIN.fasta fastaOUT.fasta fastaREJECTED.fasta"

'''


import sys
from Bio import SeqIO

## Function check that the sequence only contains ATGC nucleotides
dna = set("ATGC")
def validate(seq, alphabet=dna):
    "Checks that a sequence only contains values from an alphabet"
    leftover = set(seq.upper()) - alphabet
    
    return not leftover

# using it with other alphabets
#prot = set('ACDEFGHIKLMNPQRSTVWY')
#print validate("mglsdgewql", alphabet=prot)


## Function split sequences in valid and rejected. Also return sequence in CAPITAL
def check_fasta(fin_fname, fout_fname, fout_rejected):

    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    frejected = open(fout_rejected, "w")
    for seq_record in SeqIO.parse(fin, "fasta"):
        header = seq_record.id
        seq = seq_record.seq.upper()
        #print(seq_record.id)
        #print(repr(seq_record.seq.upper()))
        #print seq_record.seq.upper()
        sequence = seq_record.seq.upper()
        #print(len(seq_record))
        if validate(sequence):
            #print sequence

            # Write new fasta + the header
            fout.write(">"+str(header) + '\n')
            fout.write(str(seq) + '\n')
        else:

            ## Reject sequences print o screen and save to  rejected.fata file
            print "NNNNNNNNNNNNNNN         " + sequence
            frejected.write(">"+str(header) + '\n')
            frejected.write(str(seq) + '\n')
               
    


# Argumet passed to this script when called
if sys.argv.__len__() ==4:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    fout_rejected = sys.argv[3]
    check_fasta(fin_fname, fout_fname, fout_rejected)
else:
    #print str(sys.argv.__len__())
    print "error:\t4 \n Usage: \t python check_test_sequence fastaIN fastaOUT fastaREJECTED"  
    


'''


fin_fname = '/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/All_Aluexons_3SS_20_3_corrected/test/All_Aluexons_3SS_20_3_corrected_hg38_to_rheMac3.fasta'
fout_fname = '/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/All_Aluexons_3SS_20_3_corrected/test/All_Aluexons_3SS_20_3_corrected_hg38_to_rheMac3_TEST_OUT_DELETE.fasta'
fout_rejected = '/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/All_Aluexons_3SS_20_3_corrected/test/All_Aluexons_3SS_20_3_corrected_hg38_to_rheMac3_TEST_OUT_REJECTED_DELETE.fasta'


'''