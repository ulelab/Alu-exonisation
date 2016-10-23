#!/usr/bin/env Rscript

## Created on 8th of may  2016
## @author: Igor Ruiz de los Mozos


#### R script will merge all previous data comming from Lift over in a LONG and WIDE table
#     On the wide table added lots of metadata 
#
# and produce all the clasification used on this study:
#
#   whole_final.LONG.tab    -> All the exons in LONG fortat tab separated 
#   
#       chr  start  end	aluexon	position	strand	distance_to_alu	3SSScore	LongestU	LongestUright	LongestULeft
#
#
#   whole_final_5sp.WIDE.RICH.tab     -> All the exon in WIDE format plus all the rich metada of each Alu exon
##


library(ggplot2)
require(gplots)
require(reshape)
library(system)
library(parallel)

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("GMD")) {
  install.packages("GMD", dependencies = TRUE)
  library(GMD)
}


## get individual data comming from each Alu exon on different genomes version
## Read a tab files that have on each line:
## chr  start	end	aluexon	position	strand	distance_to_alu	3SSScore	LongestU	LongestUright	LongestULeft

setwd('./Data')  ## Tables storage folder

hg19 <- read.table(("Alu_exons/All_Aluexons_3SS_hg19_Distance_hg19_hg19.tab"), sep="\t",  header = TRUE)
hg19$region <- "hg19"

panTro4<- read.table("Alu_exons/All_Aluexons_3SS_hg19_Distance_hg38_panTro4.tab", sep="\t", header = TRUE)
panTro4$region <- "panTro4"

rheMac3<- read.table("Alu_exons/All_Aluexons_3SS_hg19_Distance_hg38_rheMac3.tab", sep="\t",  header = TRUE)
rheMac3$region <- "rheMac3"

calJac3<- read.table("Alu_exons/All_Aluexons_3SS_hg19_Distance_hg19_calJac3.tab", sep="\t",  header = TRUE)
calJac3$region <- "calJac3"

nomLeu1<- read.table("Alu_exons/All_Aluexons_3SS_hg19_Distance_hg19_nomLeu1.tab", sep="\t",  header = TRUE)
nomLeu1$region <- "nomLeu1"


## Merge all 3´ss characteristis in a LONG table (one after the other)
whole <- rbind(hg19, panTro4, nomLeu1, rheMac3, calJac3)
#whole <- rbind(hg19, panTro4, panPan1, nomLeu1, papHam1, rheMac3, calJac3) 

colMax <- function(whole) sapply(whole, max, na.rm = TRUE) ## get biger value in a column
colMin <- function(whole) sapply(whole, min, na.rm = TRUE) ## get min value in a column
colMax(whole$X3SSS)

## Write LONG table
write.table(whole_final, file = "Results/whole_final.LONG.tab", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)


## Same than previous table reading but this time in WIDE format
 ## Asign column names
names(hg19)<- c("chr_hg19",  "start_hg19",  "end_hg19",  "aluexon",	"position_hg19",	"strand_hg19",	"distance_hg19",	"3SSS_hg19",	"WU_hg19",	"U1_hg19",	"U2_hg19", "region_hg19")

names(panTro4)<- c("chr_panTro4",	"start_panTro4",	"end_panTro4",	"aluexon",	"position_panTro4",	"strand_panTro4",	"distance_panTro4",	"3SSS_panTro4",	"WU_panTro4",	"U1_panTro4",	"U2_panTro4", "region_panTro4")

names(rheMac3)<- c("chr_rheMac3",	"start_rheMac3",	"end_rheMac3",	"aluexon",	"position_rheMac3",	"strand_rheMac3",	"distance_rheMac3",	"3SSS_rheMac3",	"WU_rheMac3",	"U1_rheMac3",	"U2_rheMac3", "region_rheMac3")

names(calJac3)<- c("chr_calJac3",	"start_calJac3",	"end_calJac3",	"aluexon",	"position_calJac3",	"strand_calJac3",	"distance_calJac3",	"3SSS_calJac3",	"WU_calJac3",	"U1_calJac3",	"U2_calJac3", "region_calJac3")

names(nomLeu1)<- c("chr_nomLeu1",	"start_nomLeu1",	"end_nomLeu1",	"aluexon",	"position_nomLeu1",	"strand_nomLeu1",	"distance_nomLeu1",	"3SSS_nomLeu1",	"WU_nomLeu1",	"U1_nomLeu1",	"U2_nomLeu1", "region_nomLeu1")

genome_list <-list(hg19, panTro4, nomLeu1, rheMac3, calJac3)

## Open rich table on Alu exon metadata to complement the 3´ss score and the U track lenght 
total <- read.table("Alu_exons/All_Aluexons_3SS_hg19_DQ.bed", sep="\t")
colnames(total)<- c("chr_hg19_t",  "start_hg19_t",	"end_hg19_t",	"aluexon",	"position_hg19_t",	"strand_hg19_t",	"distance_hg19_t",	"chr_alu_hg19_t",	"start_alu_hg19_t", "end_alu_hg19_t", "alu_element_hg19_t", "alu_substitutions_hg19_t",	"alu_strand_hg19_t",	
                    "number")
total$number<-NULL

## Merge rich metadata table with 3'ss score and U track of each Alu exon lifter Over fom human to all the ederly primates
library(plyr)
for (genome in genome_list) {
  
  total <- merge.data.frame(total, genome, by="aluexon", all=TRUE) 
  
  
}



#################################
## Get the furthest specie in were an Alu exon have been lift over
#################################


# temporary file
rnames <- total[,1]                                            # Get line names
temp <- data.frame(total[,seq(24, ncol(total), 11)])           # Get column with 3´SSScore
rownames(temp) <- rnames

## Get the furtherst specie to fifted over

temp1 <- data.frame(total[,1])      # Get line names


## Go over each line and add 3´ss score on a list
for (i in 1:nrow(temp)) {
  lista <- list()  
  lista <- unname(unlist(temp[i,]))
  lista_n <- lista[!is.na(lista)]
  temp1$furthest_all[i] <- as.character(list(lista))
  temp1$furthest[i] <- as.character(tail(lista_n, n=1))       # Last element on the list is the FURTHEST
  
}
head(temp1)
class(temp1)

## Add FURTHEST info on table
names(temp1) <- c("aluexon", "furthest_all", "furthest")
total_furthest <- merge.data.frame(total, temp1, by="aluexon", all=TRUE) 
total_furthest <- cbind(total, temp1)

## Write WIDE table
write.table(total_furthest, file = "Results/whole_final_5sp.WIDE.tab.txt", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)


#################################
## Import UCSC Alu exon type clasification
#################################

## Import definitions of splice type exons (Constitutive, Alternative .... Cryptic...and NAs!!)
load("./Data/all.known.Aluexons.withID_andsplicetype.RData")
#unique(Igorsexons.gr$UCSCtype)
## Merege data frames to get splice type exons in total_furthest 
temp <- data.frame(Igorsexons.gr$AluexonID)
temp <- cbind(temp, as.character(Igorsexons.gr$UCSCtype))
temp <- temp[,c(4,6)]
class(temp)
typeof(temp)
showMethods(temp)
as.character(data.frame(c(Igorsexons.gr$AluexonID, Igorsexons.gr$UCSCtype)))


head(temp)
names(temp) <- c("aluexon", "UCSCtype")
colnames(temp)
ncol(total_furthest)
total_furthest <- total_furthest[,-c(69)]
total_furthest <- merge.data.frame(total_furthest, temp, by="aluexon", all=TRUE) 


#################################
## Position of 3´ss Righr or left Alu arm
#################################

total_furthest$dist <- NULL
total_furthest$dist[total_furthest$distance_hg19_t>70] <- "Left arm"
total_furthest$dist[total_furthest$distance_hg19_t<=70] <- "Right arm"
total_furthest$dist <- as.factor(total_furthest$dist) 


final_data_frame_exon_clasification <- subset(final_data_frame, final_data_frame$region == "hg19")
final_data_frame_exon_clasification <- final_data_frame_exon_clasification[,c(1,13)]
names(final_data_frame_exon_clasification) <- c("aluexon", "Exon_Originated")

temp <- merge.data.frame(total_furthest, final_data_frame_exon_clasification, by="aluexon", all=TRUE)

total_furthest <- temp



#################################
## Clasification of Alu exon by the evolutive path. 
#################################

# We classified the Alu elements based on the evolutionary dynamics of their 3’ SS.
# ‘Emerging’ Alu exons have a 3’ss with a score less than 3 in the most distant species.
# ‘Stable’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength increased towards human by less than 1.
# ‘Evolving’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength in human is more than 1+(score in the distant species). For example, if the score in marmoset is 2.5, then the Alu exon is considered as ‘emerging’, if it is 4 in marmoset and 4.5 in human, then it’s considered as ‘stable’, and if it’s 4 in marmoset and 6 in human, then it’s considered as ‘evolving’.



######## Alu exons present in marmoset

## CalJack 3ss higher than 3
cal_high_3ss <- subset(queries, queries$"X3SSS_calJac3" >=3) # 1406
## Create new column 
cal_high_3ss$human_minus_marmoset <- cal_high_3ss$X3SSS_hg19 - cal_high_3ss$X3SSS_calJac3
## 
higher_marmoset <- subset(cal_high_3ss, cal_high_3ss$human_minus_marmoset <= 1 )
higher_human <- subset(cal_high_3ss, cal_high_3ss$human_minus_marmoset > 1 )
drops <- c("human_minus_marmoset")
higher_marmoset <- higher_marmoset[, !(names(higher_marmoset) %in% drops)]
higher_human <- higher_human[, !(names(higher_human) %in% drops)]


cal_lower_3ss <- subset(queries, queries$"X3SSS_calJac3" <3)
## higher_human <- rbind(cal_lower_3ss, higher_human)  ### Jernej says to join this two toghether but I did not

higher_marmoset_l <- wide_to_long(higher_marmoset, "Constant Exon")
higher_human_l <- wide_to_long(higher_human, "Evolving Exon")
cal_lower_3ss_l <- wide_to_long(cal_lower_3ss, "Emerging Exon")

final_data_frame <- rbind(cal_lower_3ss_l, higher_human_l,  higher_marmoset_l)

exons_originated_marmoset <- final_data_frame

write.table(exons_originated_marmoset, file = "Results/exons_originated_marmoset.tab", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)


######## Exons originated in Rhesus or Baboon papHam1"="Baboon", "rheMac3

rest1 <- queries[!(rownames(queries) %in% rownames(cal_high_3ss)) & !(rownames(queries) %in% rownames(cal_lower_3ss)),] # & !(rownames(queries) %in% rownames(mac_high_3ss_no_oldMonkey)) & !(rownames(queries) %in% rownames(mac_low_3ss_no_oldMonkey)),] #1278

## Check that the number of row is the same in all of them
nrow(queries) == nrow(cal_high_3ss) + nrow(cal_lower_3ss) + nrow(rest1) #+ nrow(mac_low_3ss_no_oldMonkey) + nrow(newers_primates) + nrow(old_monkeys)

## Rhesus 3ss higher than 3
rhes_high_3ss <- subset(rest1, rest1$X3SSS_rheMac3 >= 3 ) 

## Create new column 3ss human - (3ss rhesus OR 3ss baboon) < 1

rhes_high_3ss$human_minus_rhes_papH  <- (rhes_high_3ss$X3SSS_hg19  - rhes_high_3ss$X3SSS_rheMac3)


higher_rhes <- subset(rhes_high_3ss, rhes_high_3ss$human_minus_rhes_papH < 1)
higher_human <- subset(rhes_high_3ss, rhes_high_3ss$human_minus_rhes_papH >= 1)
## Remove new column
drops <- c("human_minus_rhes_papH")
higher_human_rhes <- higher_human[, ! (names(higher_human) %in% drops)]
drops <- c("human_minus_rhes_papH")
higher_rhes <- higher_rhes[, !(names(higher_rhes) %in% drops)]

rhes_papH_lower_3ss <- subset(rest1, rest1$X3SSS_rheMac3 < 3)

### Two groups to plot rhes_papH_lower_3ss and higher_human_rhes_papH

higher_rhes_l <- wide_to_long(higher_rhes, "Constant Exon")
higher_human_rhes_l <- wide_to_long(higher_human_rhes, "Evolving Exon")
rhes_papH_lower_3ss_l <- wide_to_long(rhes_papH_lower_3ss, "Emerging Exon")

final_data_frame <- rbind(final_data_frame, rhes_papH_lower_3ss_l, higher_human_rhes_l, higher_rhes_l )  ## Join with previous clasificatin in Marmoset

exons_originated_rhesus <- final_data_frame

write.table(exons_originated_rhesus, file = "Results/exons_originated_rhesus.tab", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)



######## Exons originated in Chimp or Gibon ## Exons originated in Rhesus or Baboon papHam1"="Baboon", "rheMac3

rest2 <- queries[!(rownames(queries) %in% rownames(cal_high_3ss)) & !(rownames(queries) %in% rownames(cal_lower_3ss)) & !(rownames(queries) %in% rownames(higher_rhes)) & !(rownames(queries) %in% rownames(higher_human_rhes)) & !(rownames(queries) %in% rownames(rhes_papH_lower_3ss)),] # & !(rownames(queries) %in% rownames(mac_high_3ss_no_oldMonkey)) & !(rownames(queries) %in% rownames(mac_low_3ss_no_oldMonkey)),] #1278

## Check that the number of row is the same in all of them
nrow(queries) == nrow(cal_high_3ss) + nrow(cal_lower_3ss) + nrow(rhes_high_3ss) + nrow(rhes_papH_lower_3ss) + nrow(rest2) #+ nrow(mac_low_3ss_no_oldMonkey) + nrow(newers_primates) + nrow(old_monkeys)

## Chimp Bonobo or Gibon 3ss higher than 3
chimp_gibo_high_3ss <- subset(rest2, rest2$X3SSS_nomLeu1 >= 3 | rest2$X3SSS_panTro4 >= 3) 

rnames <- chimp_gibo_high_3ss[,1]
temp <- data.frame(chimp_gibo_high_3ss[,seq(35, 46, 11)])
rownames(temp) <- rnames

temp[1,]

## Get the furtherst specie to fifted over
temp1 <- data.frame(chimp_gibo_high_3ss[,1])
for (i in 1:nrow(temp)) {
  lista <- list()  
  lista <- unname(unlist(temp[i,]))
  lista_n <- lista[!is.na(lista)]
  temp1$furthest_all_rhes_papH[i] <- as.character(list(lista))
  temp1$furthest_rhes_papH[i] <- as.character(tail(lista_n, n=1))
  
}
head(temp1)



names(temp1) <- c("aluexon", "furthest_all_chimp_gibo", "furthest_chimp_gibo")
total_furthest_chimp_gibo <- merge.data.frame(chimp_gibo_high_3ss, temp1, by="aluexon", all=TRUE) 


## Create new column 3ss human - (3ss the oldest of chimp OR bono OR gibo) < 1
for (i in 1:nrow(total_furthest_chimp_gibo)){
  if (total_furthest_chimp_gibo$furthest_chimp_gibo[i]  == "nomLeu1"){
    
    total_furthest_chimp_gibo$human_minus_chimp_gibo[i]  <- (total_furthest_chimp_gibo$X3SSS_hg19[i]  - total_furthest_chimp_gibo$X3SSS_nomLeu1[i] )
    
  }else if (total_furthest_chimp_gibo$furthest_chimp_gibo[i]  == "panTro4"){
    total_furthest_chimp_gibo$human_minus_chimp_gibo[i]  <- (total_furthest_chimp_gibo$X3SSS_hg19[i]  - total_furthest_chimp_gibo$X3SSS_panTro4[i] )
    
  }else{
    total_furthest_chimp_gibo$human_minus_chimp_gibo[i]  <- 'NA'
  }
}


higher_chimp_gibo <- subset(total_furthest_chimp_gibo, total_furthest_chimp_gibo$human_minus_chimp_gibo < 1)
higher_human_chimp_gibo <- subset(total_furthest_chimp_gibo, total_furthest_chimp_gibo$human_minus_chimp_gibo >= 1)
## Remove new column
drops <- c("human_minus_chimp_gibo", "aluexon.1", "aluexon.2" , "furthest_all_chimp_gibo", "furthest_chimp_gibo")
higher_chimp_gibo <- higher_chimp_gibo[, !(names(higher_chimp_gibo) %in% drops)]

higher_human_chimp_gibo <- higher_human_chimp_gibo[, !(names(higher_human_chimp_gibo) %in% drops)]

rest <- rest2[!(rownames(rest2) %in% rownames(chimp_gibo_high_3ss)),]
chimp_gibo_lower_3ss <- subset(rest, rest$X3SSS_nomLeu1 < 3 | rest$X3SSS_panTro4 < 3 )

higher_chimp_gibo_l <- wide_to_long(higher_chimp_gibo, "Constant Exon")
higher_human_chimp_gibo_l <- wide_to_long(higher_human_chimp_gibo, "Evolving Exon")
chimp_gibo_lower_3ss_l <- wide_to_long(chimp_gibo_lower_3ss, "Emerging Exon")

final_data_frame <- rbind(final_data_frame, chimp_gibo_lower_3ss_l, higher_human_chimp_gibo_l, higher_chimp_gibo_l)

exons_originated_chimp_gibbon <- final_data_frame

write.table(exons_originated_chimp_gibbon, file = "Results/exons_originated_chimp_gibbon.tab", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)



#################################
## SAVE final table
#################################

write.table(total_furthest, file = "Results/whole_final_5sp.WIDE.RICH.tab", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)


