# Alu exon evolution on primate lineage


### 1.  Main Pipeline

This pipeline produces all the analysis and plots from Alu evolution on primate lineage. Each script contains detailed information on usage and each action is commented.

All the programs here described are piped by `bash run.sh inputBedFile.bed` script  [run.sh](run.sh), that will start with some data reordering, filtration and  3´ss testing, antisense Alu selection, sort and get uniques. Those unique positions were intersected with Alu TE from repeat masker. Then we measure the distance from the 3´ss to the Alu TE end.

Same process were taken for random Alus.

Subsequently it calls next sections script until plots the results.


#### Usage:
`bash run.sh bedfile.bed`

#### Scripts called
|Link|Description|
|---|-----------|
|[ get_start_position_from_bed_corrected.py ](scr/python/get_start_position_from_bed_corrected.py )| get the start position from a bed file. Used to get the 3'ss of Alu exonsCheck the best 3SS score in 2nt upstream and 2nt downstream     |
|[ get_best_3SS.py  ](scr/python/get_best_3SS.py )|   Check the best 3SS score in 2nt upstream and 2nt downstream   |
|[ add_ID_to_bed.py ](scr/python/add_ID_to_bed.py )|   ADD ID to column 4th. Id is formed by the joining of Bed positions   |
|[ 3SS_distance_from_alu.py ](scr/python/3SS_distance_from_alu.py )|  The script will get the distance between 3ss and Aluexon start     |
|[ get3SS_from_random_alu.py ]( scr/python/get3SS_from_random_alu.py )|   get the 3'SS 20 nt inside alu element    |


### 2.  Lift Over and process full table

Main script [ lift_and_procces.sh  ](scr/bash/lift_and_procces.sh )

Script is feed by a bed file with Alu exons 3´ss position separated by tab:

    chr      start      end     Alu_exon_ID      Alu_class      strand      3´ss_distance_to_Alu

#
We then used the UCSC Genome Browser LiftOver tool to obtain orthologue genomic loci of the 3’ splice site of human Alu elements in the representative species for new world (marmoset, calJac3), old world (rhesus macaque, rheMac3) and hominidae (gibbon, nomLeu1 and chimpanzee, panTro4) lineages.

This script will do:

1. Lift over the 3´ss on all the genomes of the study
2. Split the bed file on individual files. Each bed line to a different bed file
3. Get MaxEntSplice site score.  Predicted max entropy used
4. Get fasta sequence, check that is correct and measure the longest U streech  - Whole alu
5. Get fasta sequence, check that is correct and measure the longest U streech  - right arm
6. Get fasta sequence, check that is correct and measure the longest U streech  - left arm
7. Return a tabular table with all of those results for each specie:
#

        chr     start   end      aluexon    position    strand      distance_to_alu     X3SSS   LongestUTrack   UTrack_Left     UTrack_right

#
#### Usage:
Script is feed by a bed file with Alu exons 3´ss position separated by tab and with the path of the output directory.
` bash ./src/bash/lift_and_procces.sh Aluexons_3SS_hg19_Distance.bed OutDIR  `

#### Scripts called
|Link|Description|
|---|-----------|
|[ lift_over_specie.py]( scr/python/lift_over_specie.py)|  Lift over bed file to a new bed file specifying the specie conversion. Optional flag to get fasta from those liftovers Seq|
|[ split_bed_record.py ]( scr/python/split_bed_record.py)|   Return and split each line in separate file with file name equal to string representing the bed position on the genome   |
|[ flankBEDpositionsStrandSpecific.py ]( scr/python/flankBEDpositionsStrandSpecific.py )|  The script will flank the region in both directions in a new bed file.   |
|[ get_fasta_species.py ]( scr/python/get_fasta_species.py )|  Get fasta sequence from a specified genome  |
|[ check_test_sequence.py ]( scr/python/check_test_sequence.py )|  Check that the fasta sequence is appropriate for downstream analysis    |
|[ findLongestStrech.py ]( scr/python/findLongestStrech.py)|  Find the longest stretch of a given letter in a string (case insensitive)    |
|[ get_aluexon_from_distance_from_alu2.py ]( scr/python/get_aluexon_from_distance_from_alu2.py)|  The script will get the distance between 3ss and Alu start. Then it will create a synthetic Alu element covering all the predicted Alu seq (at least 320 nt).    |


### 3.  Data pre-processing

Main script [  Data_Prepocesig.R ]( /scr/R/Data_Prepocesig.R )

R script will merge all previous data coming from Lift over in a LONG and WIDE format tables.

- On the LONG table added lots of metadata and produce all the classification used on this study:

        chr     start   end	    aluexon	    position	strand	    distance_to_alu	3SSScore	LongestUTrack	LongestUright	LongestULeft

[ whole_final.LONG.tab ](/Results/whole_final.LONG.tab)    -> All the exons in LONG fortat tab separated

- On the WIDE format table, every row belong to a bed position in human and on the right columns the rest of the species


        chr     start   end	    aluexon	    position	strand	    distance_to_alu	3SSScore	LongestUTrack	LongestUright	LongestULeft

[ whole_final_5sp.WIDE.RICH.tab ](/Results/whole_final_5sp.WIDE.RICH.tab)     -> All the exon in WIDE format plus all the rich metadata of each Alu exon



#### Usage:
` Rscript ./scr/R/Data_Prepocesig.R `

###########################



### 4. Classification of Alu elements by divergence or evolutionary dynamics, and analysis of their 3’ splice sites and U-tracts


Main script [ ** Fig_5.R ]( scr/R/Fig_5.R)


#### Clasification of Alu exon by the evolutive path.


1. First we get the most distant specie in where we could find homologous Alu sequences inserted (FURTHERST)

2. Then we classified the Alu elements based on the evolutionary dynamics of their 3’ss.

   - ‘Emerging’ Alu exons have a 3’ss with a score less than 3 in the most distant species.
   - ‘Stable’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength increased towards human by less than 1.
   - ‘Evolving’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength in human is more than 1+(score in the distant species). For example, if the score in marmoset is 2.5, then the Alu exon is considered as ‘emerging’, if it is 4 in marmoset and 4.5 in human, then it’s considered as ‘stable’, and if it’s 4 in marmoset and 6 in human, then it’s considered as ‘evolving’.

3. Plot 3 splice site strength, Longest U track on Alu, longest U track on left arm, longest U track on right arm



#### Usage:
` Rscipt .src/R/Fig_5.R `


###########################




### 5.  3´ss coupled with U track lenght

Main script [ ** Fig_5.R ]( scr/R/Fig_5.R)


#### Clasification of Alu exon by the evolutive path.


1. First we get the most distant specie in where we could find homologous Alu sequences inserted (FURTHERST)

2. Then we classified the Alu elements based on the evolutionary dynamics of their 3’ss.

   - ‘Emerging’ Alu exons have a 3’ss with a score less than 3 in the most distant species.
   - ‘Stable’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength increased towards human by less than 1.
   - ‘Evolving’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength in human is more than 1+(score in the distant species). For example, if the score in marmoset is 2.5, then the Alu exon is considered as ‘emerging’, if it is 4 in marmoset and 4.5 in human, then it’s considered as ‘stable’, and if it’s 4 in marmoset and 6 in human, then it’s considered as ‘evolving’.

3. Plot 3 splice site strength, Longest U track on Alu, longest U track on left arm, longest U track on right arm



#### Usage:
` Rscipt .src/R/Fig_5.R `


###Source Code Overview
![module diagram](Data/Structure.png "Source Code Overview")
