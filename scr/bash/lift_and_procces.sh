#!/bin/bash -l

## Created on 15 February 2016
## @author: Igor Ruiz de los Mozos

#################################################
############## Main Script of this research
#################################################
##
##  Script is feed by a bed file with Alu exons 3´ss position. chr      start   end     Alu_exon_ID    Alu_class    strand  3´ss_distance_to_Alu
##
##  It will do:
##                  1) Lift over the 3´ss to diferent genomes
##                  2) Split the bed file on individual files. Each bed line to a different bed file
##                  3) Get MaxEntSplice site score
##                  4) Get fasta sequence, check that is correct and measure the longest U streech  - Whole alu
##                  5) Get fasta sequence, check that is correct and measure the longest U streech  - right arm
##                  6) Get fasta sequence, check that is correct and measure the longest U streech  - left arm
##                  7) Return a tabular table with all of those results
##
##  Output file will be named as the input file but end on .tab
##
##
## Usage:
##          ./lift_and_procces.sh Aluexons_3SS_hg19_Distance.bed OutDIR




## First argument passed to the script is the bed file containing the distance from Alu on the 7th column
FILES=$1

## Create temporary directory. Remove previous data
OUTDIR=OUTDIR_Random
rm -r ./${OUTDIR}
mkdir -p ./${OUTDIR}


for file in $FILES; do

	filename=$(basename "$file")
	extension="${filename##*.}"
	filename="${filename%.*}"
	echo "$file"
	echo "$filename"
	echo "$extension"


	DIR=lifted
	rm -r ./${DIR}
	mkdir -p ./${DIR}


    ##  1) Lift over the 3´ss to diferent genomes

    ## Prepare file
	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' ${file} > ./${DIR}/${filename}_C.bed

	python './src/python/lift_over_specie.py' ./${DIR}/${filename}_C.bed hg19 hg38 ./${DIR}/${filename}_hg19_hg38.bed FALSE

	python './src/python/lift_over_specie.py' ./${DIR}/${filename}_hg19_hg38.bed hg38 panTro4 ./${DIR}/${filename}_hg38_panTro4.bed FALSE
	python './src/python/lift_over_specie.py' ./${DIR}/${filename}_hg19_hg38.bed hg38 panPan1 ./${DIR}/${filename}_hg38_panPan1.bed FALSE
	python './src/python/lift_over_specie.py' ./${DIR}/${filename}_hg19_hg38.bed hg38 rheMac3 ./"$DIR"/${filename}_hg38_rheMac3.bed FALSE
	python './src/python/lift_over_specie.py' ./${DIR}/${filename}_hg19_hg38.bed hg38 tarSyr2 ./"$DIR"/${filename}_hg38_tarSyr2.bed FALSE
	python './src/python/lift_over_specie.py' ./${DIR}/${filename}_C.bed hg19 otoGar1 ./"$DIR"/${filename}_hg19_otoGar1.bed FALSE
	python './src/python/lift_over_specie.py' ./${DIR}/${filename}_C.bed hg19 micMur1 ./"$DIR"/${filename}_hg19_micMur1.bed FALSE
	python './src/python/lift_over_specie.py' ./${DIR}/${filename}_C.bed hg19 tupBel1 ./"$DIR"/${filename}_hg19_tupBel1.bed FALSE
	python './src/python/lift_over_specie.py' ./${DIR}/${filename}_C.bed hg19 nomLeu1 ./"$DIR"/${filename}_hg19_nomLeu1.bed FALSE
	python './src/python/lift_over_specie.py' ./${DIR}/${filename}_C.bed hg19 papHam1 ./"$DIR"/${filename}_hg19_papHam1.bed FALSE
	python './src/python/lift_over_specie.py' ./${DIR}/${filename}_C.bed hg19 calJac3 ./"$DIR"/${filename}_hg19_calJac3.bed FALSE
	##rm ./${DIR}/${filename}_C.bed
	mv ./${DIR}/${filename}_C.bed ./${DIR}/${filename}_hg19_hg19.bed

	#### longest U and SS scores

	LIFTED_FILES=./${DIR}/*.bed
	cp -r '/home/igor/Programs/MaxEntScan/splicemodels' .    ## this folder is needed to run the program Using predicted max entropy

	for lif_file in $LIFTED_FILES; do


		lif_filename=$(basename "$lif_file")
		lif_extension="${lif_filename##*.}"
		lif_filename="${lif_filename%.*}"
		echo "NEW LIFTED FILE ##########################"
		echo "FILE    " "$lif_file"
		echo 'BaseName  ' "$lif_filename"
		echo "Directory  " "$DIR"

		## get the genome name == var6 witch have been lift over
		IFS=_ read -a arr <<< "$lif_filename"
		echo "Genome to get fasta  var7" ${arr[6]}
		echo 'BaseName  ' "$lif_filename"

		mkdir -p ./"$DIR"/"$lif_filename"


		#mkdir -p ./quantiles/${DIR}
		#DIR=quantiles/${DIR}
		DIR_BED=${DIR}/"$lif_filename"

		mkdir -p ./${DIR_BED}



        ## 2) Split the bed file on individual files. Each bed line to a different bed file


		python './src/python/split_bed_record.py' "$lif_file" ./"${DIR_BED}"/


		BEDS=./$DIR_BED/*.bed

		for bed_position in $BEDS; do

			bed_filename=$(basename "$bed_position")
			bed_extension="${bed_filename##*.}"
			bed_filename="${bed_filename%.*}"
			echo "NEW BED POSITION @@@@@@@@@@@@@@@@@@@@@@@@@@@"
			echo "FILE    " "$bed_position"
			echo 'BaseName  ' "$bed_filename"


            ## Get 19 nt upstream 3´s and 3 downstream on bed file and then get fasta file
			python './src/python/flankBEDpositionsStrandSpecific.py' "$bed_position" ./"$DIR_BED"/"$bed_filename"-19_3.bed 19 3
			python './src/python/get_fasta_species.py' ./"$DIR_BED"/"$bed_filename"-19_3.bed ${arr[6]} ./"$DIR_BED"/"$bed_filename".fasta

			## Validate is DNA without Ns
			python './src/python/check_test_sequence.py' ./"$DIR_BED"/"$bed_filename".fasta ./"$DIR_BED"/"$bed_filename"_valid.fasta ./"$DIR_BED"/"$bed_filename"_REJECTED.fasta

			## Get file with longest U track
			python './src/python/findLongestStrech.py' ./"$DIR_BED"/"$bed_filename"_valid.fasta ./"$DIR_BED"/"$bed_filename"_valid.tab t

            ##  3) Get MaxEntSplice site score
			## Get the MaxEntScan score on each fasta seq (Alu 3´ss)
			perl /home/igor/Programs/MaxEntScan/score3.pl ./"$DIR_BED"/"$bed_filename"_valid.fasta > ./"$DIR_BED"/"$bed_filename"_scores.tab


            ## Get the distance from the 3´ss to the end of Alu retrotransposable element
			python './src/python/get_aluexon_from_distance_from_alu2.py' "$bed_position" ./"$DIR_BED"/"$bed_filename"-320.bed

			## Flank 250 upstream and 70 downstream
			python './src/python/flankBEDpositionsStrandSpecific.py' ./"$DIR_BED"/"$bed_filename"-320.bed ./"$DIR_BED"/"$bed_filename"-70.bed 0 -250
			python './src/python/flankBEDpositionsStrandSpecific.py' ./"$DIR_BED"/"$bed_filename"-320.bed ./"$DIR_BED"/"$bed_filename"-250.bed -70 0


            ##  4) Get fasta sequence, check that is correct and measure the longest U streech  - Whole alu

			## Whole aluexon
			python './src/python/get_fasta_species.py' ./"$DIR_BED"/"$bed_filename"-320.bed ${arr[6]} ./"$DIR_BED"/"$bed_filename"-320.fasta
			python './src/python/check_test_sequence.py' ./"$DIR_BED"/"$bed_filename"-320.fasta ./"$DIR_BED"/"$bed_filename"_valid-320.fasta ./"$DIR_BED"/"$bed_filename"_REJECTED-320.fasta
			python './src/python/findLongestStrech.py' ./"$DIR_BED"/"$bed_filename"_valid-320.fasta ./"$DIR_BED"/"$bed_filename"_valid-320.tab t

            ## 5) Get fasta sequence, check that is correct and measure the longest U streech  - right arm

			## U1 region is the alu right arm
			python './src/python/get_fasta_species.py' ./"$DIR_BED"/"$bed_filename"-70.bed ${arr[6]} ./"$DIR_BED"/"$bed_filename"-70.fasta
			python './src/python/check_test_sequence.py' ./"$DIR_BED"/"$bed_filename"-70.fasta ./"$DIR_BED"/"$bed_filename"_valid-70.fasta ./"$DIR_BED"/"$bed_filename"_REJECTED-70.fasta
			python './src/python/findLongestStrech.py' ./"$DIR_BED"/"$bed_filename"_valid-70.fasta ./"$DIR_BED"/"$bed_filename"_valid-70.tab t

            ## 5) Get fasta sequence, check that is correct and measure the longest U streech  - left arm

			## U2 region is the alu left arm
			python './src/python/get_fasta_species.py' ./"$DIR_BED"/"$bed_filename"-250.bed ${arr[6]} ./"$DIR_BED"/"$bed_filename"-250.fasta
			python './src/python/check_test_sequence.py' ./"$DIR_BED"/"$bed_filename"-250.fasta ./"$DIR_BED"/"$bed_filename"_valid-250.fasta ./"$DIR_BED"/"$bed_filename"_REJECTED-250.fasta
			python './src/python/findLongestStrech.py' ./"$DIR_BED"/"$bed_filename"_valid-250.fasta ./"$DIR_BED"/"$bed_filename"_valid-250.tab t




            ################################################
            ###### Create tab file with all the results ########
            ################################################

            ##  6) Return a tabular table with all of previous results

            ## Create file
			cat "$bed_position" > ./${OUTDIR}/temp_file.tab

            ## Grep 3´ss score from intermediate file on the last column
			score=0
			score=$(cat ./"$DIR_BED"/"$bed_filename"_scores.tab | awk '{print $2}')
			## Insert 3´ss score on the last column
			sed -i "s/$/\t$score/" ./${OUTDIR}/temp_file.tab

            ## Grep Longest U track from intermediate file (whole Alu exon) on the last column
			polyU_whole=0
			polyU_whole=$(cat ./"$DIR_BED"/"$bed_filename"_valid-320.tab | awk '{print $3}')

			## If cannot find longest track asign "NA" to output
			if [[ -z $polyU_whole ]]; then
				polyU_whole="NA"
			fi
			## Insert on the last columnn longest U track
			sed -i "s/$/\t$polyU_whole/" ./${OUTDIR}/temp_file.tab

            ## Same but for right arm
			U1=0
			U1=$(cat ./"$DIR_BED"/"$bed_filename"_valid-70.tab | awk '{print $3}')
			if [[ -z $U1 ]]; then
				U1="NA"
			fi
			sed -i "s/$/\t$U1/" ./${OUTDIR}/temp_file.tab

            ## Same but for left arm
			U2=0
			U2=$(cat ./"$DIR_BED"/"$bed_filename"_valid-250.tab | awk '{print $3}')		## get the value from colum 3
			if [[ -z $U2 ]]; then									                    ## Set variable to 'NA' if its empty
				U2="NA"
			fi
			sed -i "s/$/\t$U2/" ./${OUTDIR}/temp_file.tab						        ## Add to the end


			cat ./${OUTDIR}/temp_file.tab >> ./${OUTDIR}/"$lif_filename".tab			## Apend to output file
			rm -r ./${OUTDIR}/temp_file.tab                                             ## Remove intermediate file


		#rm -r ./${DIR_BED}
		done

	sed -i '1ichr\tstart\tend\taluexon\tposition\tstrand\tdistance_to_alu\tX3SSS\tWU\tU1\tU2' ./${OUTDIR}/"$lif_filename".tab   		## Add headers to the file separated by tab
	done ## Nested loop
	rm -r './splicemodels'

done




