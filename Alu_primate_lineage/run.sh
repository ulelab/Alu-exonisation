#!/bin/bash -l

### This pipeline will produce all the analysis and the plots from this evolutionary analysis on primate lineage
### Each script contain detailed usage information and acctions are commented
### Refer to /src/python or /src/bash /src/R to find all the script used

# Input bed file with Alu exon positions
bed_file=$1

# Create Data folder
mkdir Data

### Feed with the raw table of Allu exon positions
awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $bed_file > './Data/Alu_exons/All_Aluexons.bed'

### Get the 3SS  -2 in positive strand   !!! I detected that file have 2nt mismacth incorrect position on the positive strand. With this is corrected
python './src/python/get_start_position_from_bed_corrected.py' './Data/All_Aluexons.bed' './Data/All_Aluexons_3SS.bed'

########## TEST 3SS localization
### NEGATIVE
cat './Data/All_Aluexons_3SS.bed' | grep '-' > './Data/All_Aluexons_3SS_AS_sort_negative_strand.bed'
### Get the best 3SS fromt 2nt upstream to 2 nt downstream
python './src/python/get_best_3SS.py' './Data/All_Aluexons_3SS_AS_sort_negative_strand.bed' './Data/All_Aluexons_3SS_AS_sort_negative_strand_corrected.bed' './temp/'

# Output obtained
#2 downstream 211
#1 downstream 161
#same 3SS     2610
#1 upstream   134
#2 upstream   242
#Total        3358
#BED REJECTED 0


### POSITIVE
cat './Data/All_Aluexons_3SS.bed' | grep '+' > './Data/All_Aluexons_3SS_AS_sort_positive_strand.bed'
### Get the best 3SS fromt 2nt upstream to 2 nt downstream
python './src/python/get_best_3SS.py' './Data/All_Aluexons_3SS_AS_sort_positive_strand.bed' './Data/All_Aluexons_3SS_AS_sort_positive_strand_corrected.bed' './temp/'

#Outpput obtained on Positive strand
#2 downstream 216
#1 downstream 144
#same 3SS     2778
#1 upstream   154
#2 upstream   209
#Total        3501
#BED REJECTED 3



####################################################
#### Only use the ones that are in antisense
####################################################

### Intersect with Alus to check if they still belong to them -S Require different strandedness 
bedtools intersect -S -wa -a './Data/All_Aluexons_3SS.bed' -b '/media/igor/DATA/UCL/Evolution_Alus/Raw_Alus/rmsk_hg19_full_family_Alu_elements.bed' > './Data/All_Aluexons_3SS_AS.bed'
### Sort and get uniques 
sort -k 1,1 -k2,2n -k3,3n -k6,6 -u './Data/All_Aluexons_3SS_AS.bed' > './Data/All_Aluexons_3SS_AS_sort.bed'


### Add Id to bed file. Consisting of all the bed information join in column 4
python './src/python/add_ID_to_bed.py' './Data/All_Aluexons_3SS_AS_sort.bed' './Data/All_Aluexons_3SS_hg19_ID.bed'


### File to get the distance from the 3'SS to alu element. It intersect bed file with 3´ss position and returns this distance
bedtools intersect -wao -a './Data/All_Aluexons_3SS_hg19_ID.bed' -b '/media/igor/DATA/UCL/Evolution_Alus/Raw_Alus/rmsk_hg19_full_family_Alu_elements.bed' > './Data/All_Aluexons_3SS_hg19_Distance.bed'


### For real 3'SS 
python './src/python/3SS_distance_from_alu.py' './Data/All_Aluexons_3SS_hg19_Distance.bed' './Data/All_Aluexons_3SS_hg19_Distance_to_alu.bed'

############################################################
### 2  Lift over in all the studied genomes and create full table
############################################################

cp './Data/All_Aluexons_3SS_hg19_Distance_to_alu.bed' './Data/All_Aluexons_3SS_hg19_Distance.bed'
bash ./src/bash/lift_and_procces.sh './Data/All_Aluexons_3SS_hg19_Distance.bed'


## From the TE full table get only Alus elements
cat '/media/igor/DATA/UCL/Evolution_Alus/Raw_Alus/rmsk_hg19_full_full.tab' | grep 'Alu' > /media/igor/DATA/UCL/Evolution_Alus/Raw_Alus/Alus.hg19.tab

## Go to Rranges to get the quantiles '/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/import RptMasker output as GRange.R' output > '/media/igor/DATA/UCL/Evolution_Alus/New3SS/3SS_Alus/hg19_allAlus.bed'

##  Intersect with file with quantiles
bedtools intersect -wao -a  './Data/All_Aluexons_3SS_hg19_Distance_to_alu.bed' -b '/media/igor/DATA/UCL/Evolution_Alus/New3SS/3SS_Alus/hg19_allAlus.bed' > './Data/All_Aluexons_3SS_hg19_DQ.bed'



############################################################
### To get negative controls of alus on intronic sequences
############################################################

### Get alus that are in antisense of transcripts (
bedtools intersect -S -wa -a '/Data/intronic.Alu.elements.bed' -b '/Data/hg19-ensemble-whole_gene-v74-longest_trans.bed' > './Data/Random_Aluelements_AS.bed'

### shufle get the first 10000 and sort
shuf './Data/Random_Aluelements_AS.bed' | head -10000 | sort -k 1,1 -k2,2n -u > './Data/Random_alu_exons.bed'


### Sort and get uniques 
sort -k 1,1 -k2,2n -k3,3n -k6,6 -u './Data/Random_alu_exons.bed' > './Data/All_Aluexons_3SS_AS_sort.bed'


cat './Data/Random_alu_exons.bed' | wc -l  ### 9811 Random Alu exons


##Get a imaginary 3SS 20 nt inside the exon
python './src/python/get3SS_from_random_alu.py' './Data/Random_alu_exons.bed' './Data/Random_alu_exons_3SS_AS.bed' 20
### Get the best 3SS fromt 2nt upstream to 2 nt downstream
python './src/python/get_best_3SS.py' './Data/Random_alu_exons_3SS_AS.bed' './Data/Random_alu_exons_3SS_C.bed' './temp/'

#Check output of random Alu 3´ss exons
#2 downstream 2035
#1 downstream 1726
#same 3SS     2293
#1 upstream   2076
#2 upstream   1868
#Total        9998
#BED REJECTED 2


### File to get the distance from the 3'SS to alu element 
bedtools intersect -wao -a './Data/Random_alu_exons_3SS_C.bed' -b '/Data/rmsk_hg19_full_family_Alu_elements.bed' > './Data/Random_alu_exons_3SS_C_distance.bed'


### Get distance from alu exon start OR alu element end CORRECT distance
### For random 3'SS
python './src/python/3SS_distance_from_alu.py' './Data/Random_alu_exons_3SS_C_distance.bed' './Data/Random_alu_exons_3SS_C_distance_to_alu.bed'

## lift over in all the studied genomes and create full table (
cp './Data/Random_alu_exons_3SS_C_distance_to_alu.bed' './Data/Random_alu_exons_3SS_distance.bed'

############################################################
### 2.1  Lift Over and Process full table for Random Alu elements
############################################################

############# Main script that lift over each sequence and get 3´ss strengt, longest U track, Longest U track on right arm, longest U track on left arm
./lift_and_procces.sh './Data/Random_alu_exons_3SS_distance.bed'


############################################################
### 3. Data pre-processing
############################################################
#
# R script will merge all previous data coming from Lift over in a LONG and WIDE format tables.

Rscript /scr/R/Data_Prepocesig.R

