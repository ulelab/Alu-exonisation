# Alu-exonisation
tools related to identification &amp; characterisation of cryptic Alu-exons

These tools have been used to generate or analyse data tables analysing cryptic Alu-exons. 




U-track length
---------------

get_only_U_track_positions.py
The script returns positions of U tracks that are longer then 4 nts and smaller then 100 nt. The input is fasta format. We used the script to identify U-tracts within Alu elements giving rise to Alu exons. The output will contain any U-tract, and positions within the element can then be compared to the 3' splice site position, to identify U-tracks in front of the 3' SS. 

   

PTC prediction tools
--------------------
(1) get_transcripts_with_the_longest_CDS.py
The script will select the trancript with the longest CDS. It takes as input the transcript table available from UCSC table browser ('knownGene' track). The output generated is hg19-cds-max_exons2.bed

(2) add_exonic_positions.py
The script will add cryptic Alu-exons (or any exon) to the exon/CDS annotation table, based on position. This is identifying if the exon is contributing to the annotated CDS and which gene/transcript it belongs to. We used the output from script (1), hg19-cds-max_exons2.bed, as input. The script uses bedtools to compare annotation and custom exon positions: intersect -s -a regulated.exons.bed -b hg19-cds-max_exons2.bed -wb
The output is a table with CDS genomic positions with information of exons on transcript level together with Alu exons

(3) separate_PTC_exons2.py
The script will receive a fasta file with a header of regulated exon ID, regulated exon position relative to CDS start, all exonic positions relative to CDS start position and a FASTA sequence of CDS. The input is generated from the output of (2), using bedtools getFastafrombed.
The output of the script is a table with annotation of PTCs in all custom exons. 
PTC scoring:
Predict frame of the exon based on CDS annotation of the dominant isoform of the gene and output "PTC in frame 1/2/3" if the sequence contains an UAG,UAA,UGA within frame (i.e. it contains a PTC). If it doesnt have PTC, then check if exon length not dividable by 3. If not, then write YES (== frame-shift and thereby likely to introduce PTC). If it is dividable by 3, then write NO.


