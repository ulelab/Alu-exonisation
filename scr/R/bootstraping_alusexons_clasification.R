### Bootstrapping

library(ggplot2)
require(ggplot2)
require(reshape)
library(system)
library(boot)
library(parallel)
#library(doParallel)
#registerDoParallel(cores=detectCores(all.tests=TRUE))

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}


setwd('/Users/Igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data')
setwd('/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data')  ## Curro


hg19 <- read.table(("All_Aluexons_3SS_hg19_Distance_hg19_hg19.tab"), sep="\t",  header = TRUE)
hg19$region <- "hg19"

panTro4<- read.table("All_Aluexons_3SS_hg19_Distance_hg38_panTro4.tab", sep="\t", header = TRUE)
panTro4$region <- "panTro4"

rheMac3<- read.table("All_Aluexons_3SS_hg19_Distance_hg38_rheMac3.tab", sep="\t",  header = TRUE)
rheMac3$region <- "rheMac3"

calJac3<- read.table("All_Aluexons_3SS_hg19_Distance_hg19_calJac3.tab", sep="\t",  header = TRUE)
calJac3$region <- "calJac3"

nomLeu1<- read.table("All_Aluexons_3SS_hg19_Distance_hg19_nomLeu1.tab", sep="\t",  header = TRUE)
nomLeu1$region <- "nomLeu1"


#genome_list <-list(hg19, hg38, panTro4, panPan1, nomLeu1, papHam1, rheMac3, calJac3, tarSyr2, micMur1, otoGar1, tupBel1) 
whole <- rbind(hg19, panTro4, nomLeu1, rheMac3, calJac3)
#whole <- rbind(hg19, panTro4, panPan1, nomLeu1, papHam1, rheMac3, calJac3) 

colMax <- function(whole) sapply(whole, max, na.rm = TRUE) ## get biger value in a column
colMin <- function(whole) sapply(whole, min, na.rm = TRUE) ## get min value in a column
colMax(whole$X3SSS)

delete <- subset(whole, select=c(X3SSS, WU, U1, U2))
colMax(delete)    ## max 3sss 14.69
colMin(delete)

all_table <- data.frame(table(whole$region))
all_table 


#write.table(whole_final, file = "whole_final.bed", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)


names(hg19)<- c("chr_hg19",  "start_hg19",  "end_hg19",  "aluexon",	"position_hg19",	"strand_hg19",	"distance_hg19",	"3SSS_hg19",	"WU_hg19",	"U1_hg19",	"U2_hg19", "region_hg19")

names(panTro4)<- c("chr_panTro4",	"start_panTro4",	"end_panTro4",	"aluexon",	"position_panTro4",	"strand_panTro4",	"distance_panTro4",	"3SSS_panTro4",	"WU_panTro4",	"U1_panTro4",	"U2_panTro4", "region_panTro4")

names(rheMac3)<- c("chr_rheMac3",	"start_rheMac3",	"end_rheMac3",	"aluexon",	"position_rheMac3",	"strand_rheMac3",	"distance_rheMac3",	"3SSS_rheMac3",	"WU_rheMac3",	"U1_rheMac3",	"U2_rheMac3", "region_rheMac3")

names(calJac3)<- c("chr_calJac3",	"start_calJac3",	"end_calJac3",	"aluexon",	"position_calJac3",	"strand_calJac3",	"distance_calJac3",	"3SSS_calJac3",	"WU_calJac3",	"U1_calJac3",	"U2_calJac3", "region_calJac3")

names(nomLeu1)<- c("chr_nomLeu1",	"start_nomLeu1",	"end_nomLeu1",	"aluexon",	"position_nomLeu1",	"strand_nomLeu1",	"distance_nomLeu1",	"3SSS_nomLeu1",	"WU_nomLeu1",	"U1_nomLeu1",	"U2_nomLeu1", "region_nomLeu1")


#("hg19"="1", "hg38"="2", "panTro4"="3", "panPan1"="4", "nomLeu1"="5", "papHam1"="6", "rheMac3"="7", "calJac3"="8", "tarSyr2"="9", "micMur1"="10", "otoGar1"="11", "tupBel1"="12")  

#genome_list <-list(hg19, panTro4, panPan1, nomLeu1, papHam1, rheMac3, calJac3) 
genome_list <-list(hg19, panTro4, nomLeu1, rheMac3, calJac3)

total <- read.table("All_Aluexons_3SS_hg19_DQ.bed", sep="\t")


colnames(total)<- c("chr_hg19_t",  "start_hg19_t",	"end_hg19_t",	"aluexon",	"position_hg19_t",	"strand_hg19_t",	"distance_hg19_t",	"chr_alu_hg19_t",	"start_alu_hg19_t", "end_alu_hg19_t", "alu_element_hg19_t", "alu_substitutions_hg19_t",	"alu_strand_hg19_t",	
                    "number")
total$number<-NULL


library(plyr)
for (genome in genome_list) {
  
  total <- merge.data.frame(total, genome, by="aluexon", all=TRUE) 
  
  
}

###########
rnames <- total[,1]
temp <- data.frame(total[,seq(24, ncol(total), 11)])
rownames(temp) <- rnames

## Get the furtherst specie to fifted over
#temp <- temp[1:10,]
temp1 <- data.frame(total[,1])
#temp1 <- temp1[1:10,]
for (i in 1:nrow(temp)) {
  lista <- list()  
  lista <- unname(unlist(temp[i,]))
  lista_n <- lista[!is.na(lista)]
  temp1$furthest_all[i] <- as.character(list(lista))
  temp1$furthest[i] <- as.character(tail(lista_n, n=1))
  
}
head(temp1)
class(temp1)

names(temp1) <- c("aluexon", "furthest_all", "furthest")
total_furthest <- merge.data.frame(total, temp1, by="aluexon", all=TRUE) 
total_furthest <- cbind(total, temp1)




#setwd('/media/igor/DATA/UCL/Evolution_Alus/New3SS/Output_tables')
#write.table(total_furthest, file = "total_furthest.tab", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)


queries <- total_furthest
## Change colnames to introduce X in front of 3SS score...
names(queries)[names(queries) == '3SSS_hg19'] <- 'X3SSS_hg19'
names(queries)[names(queries) == '3SSS_panTro4'] <- 'X3SSS_panTro4'
#names(queries)[names(queries) == '3SSS_panPan1'] <- 'X3SSS_panPan1'
names(queries)[names(queries) == '3SSS_nomLeu1'] <- 'X3SSS_nomLeu1'
#names(queries)[names(queries) == '3SSS_papHam1'] <- 'X3SSS_papHam1'
names(queries)[names(queries) == '3SSS_rheMac3'] <- 'X3SSS_rheMac3'
names(queries)[names(queries) == '3SSS_calJac3'] <- 'X3SSS_calJac3'

colnames(queries)


setwd('/media/igor/DATA/UCL/Evolution_Alus/New3SS/Output_R')

wide_to_long <- function(data_frame, new_column){
  
  # Function transform the wide table of lift ovet to long table in where each row has a factor that represent its procedence and used for later plotting
  #new_column <- "cal_high_3ss"
  #data_frame<-queries
  
  colnames(data_frame) <- NULL
  hg19_q <- data.frame(data_frame[,1], data_frame[,14:24])
  pantro4_q <- data.frame(data_frame[,1], data_frame[,25:35])
  nomLeu1_q<-  data.frame(data_frame[,1], data_frame[,36:46])
  rheMac3_q <- data.frame(data_frame[,1], data_frame[,47:57])
  calJac3_q <- data.frame(data_frame[,1], data_frame[,58:68])
  
  # <- data.frame(data_frame[,1], data_frame[,69:79])
  # <- data.frame(data_frame[,1], data_frame[,80:90])
  #rest <- data.frame(data_frame[,1], data_frame[,91:101])
  
  
  
  out_data_frame <- rbind(hg19_q, pantro4_q, nomLeu1_q, rheMac3_q, calJac3_q)
  
  colnames(out_data_frame) <- c("aluexon", "chr", "start", "end", "position", "strand", "distance_to_alu", "X3SSS", "WU", "U1", "U2", "region")# , "group")  
  
  out_data_frame$group <- as.factor(new_column)
  out_data_frame <- out_data_frame[complete.cases(out_data_frame),]
  
  return(out_data_frame)
  
}

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


#cal_lower_3ss <- subset(queries, queries$"X3SSS_calJac3" <3) # 1748
#cal_lower_3ss <- rbind(cal_lower_3ss, higher_marmoset)

cal_lower_3ss <- subset(queries, queries$"X3SSS_calJac3" <3)
## higher_human <- rbind(cal_lower_3ss, higher_human)  ### Jernej says to join this two toghether but I did not

### Two groups to plot higher_marmoset and higher_human
# Maybe it will be nice to get the rest... on somethins similat to== rest <- queries[!(rownames(queries) %in% rownames(cal_high_3ss)) & !(rownames(queries) %in% rownames(cal_low_3ss)) & !(rownames(queries) %in% rownames(mac_high_3ss_no_oldMonkey)) & !(rownames(queries) %in% rownames(mac_low_3ss_no_oldMonkey)),] #1901

higher_marmoset_l <- wide_to_long(higher_marmoset, "Constant Exon")
higher_human_l <- wide_to_long(higher_human, "Evolving Exon")
cal_lower_3ss_l <- wide_to_long(cal_lower_3ss, "Emerging Exon")

final_data_frame <- rbind(cal_lower_3ss_l, higher_human_l,  higher_marmoset_l)

#final_data_frame[] <- lapply(final_data_frame, as.character)



#################################
## Exons originated in Rhesus or Baboon papHam1"="Baboon", "rheMac3
#################################
##!!!!!! Load again the previous tables before function wide to long
rest1 <- queries[!(rownames(queries) %in% rownames(cal_high_3ss)) & !(rownames(queries) %in% rownames(cal_lower_3ss)),] # & !(rownames(queries) %in% rownames(mac_high_3ss_no_oldMonkey)) & !(rownames(queries) %in% rownames(mac_low_3ss_no_oldMonkey)),] #1278

## Check that the number of row is the same in all of them
nrow(queries) == nrow(cal_high_3ss) + nrow(cal_lower_3ss) + nrow(rest1) #+ nrow(mac_low_3ss_no_oldMonkey) + nrow(newers_primates) + nrow(old_monkeys)


## Rhesus 3ss higher than 3
rhes_high_3ss <- subset(rest1, rest1$X3SSS_rheMac3 >= 3 ) 


## Create new column 3ss human - (3ss rhesus OR 3ss baboon) < 1

rhes_high_3ss$human_minus_rhes_papH  <- (rhes_high_3ss$X3SSS_hg19  - rhes_high_3ss$X3SSS_rheMac3)



# display.brewer.pal(9, "Greens")
# col <- brewer.pal(9, "Greens")[c(2:6)]

higher_rhes <- subset(rhes_high_3ss, rhes_high_3ss$human_minus_rhes_papH < 1)
higher_human <- subset(rhes_high_3ss, rhes_high_3ss$human_minus_rhes_papH >= 1)
## Remove new column
drops <- c("human_minus_rhes_papH")
higher_human_rhes <- higher_human[, ! (names(higher_human) %in% drops)]
drops <- c("human_minus_rhes_papH")
higher_rhes <- higher_rhes[, !(names(higher_rhes) %in% drops)]



rhes_papH_lower_3ss <- subset(rest1, rest1$X3SSS_rheMac3 < 3)
#higher_human_rhes_papH <- rbind(higher_human_rhes_papH, rhes_papH_lower_3ss)   ### Jernej says to join this two toghether but I did not


### Two groups to plot rhes_papH_lower_3ss and higher_human_rhes_papH
# Maybe it will be nice to get the rest... on somethins similat to== rest <- queries[!(rownames(queries) %in% rownames(cal_high_3ss)) & !(rownames(queries) %in% rownames(cal_low_3ss)) & !(rownames(queries) %in% rownames(mac_high_3ss_no_oldMonkey)) & !(rownames(queries) %in% rownames(mac_low_3ss_no_oldMonkey)),] #1901

higher_rhes_l <- wide_to_long(higher_rhes, "Constant Exon")
higher_human_rhes_l <- wide_to_long(higher_human_rhes, "Evolving Exon")
rhes_papH_lower_3ss_l <- wide_to_long(rhes_papH_lower_3ss, "Emerging Exon")


final_data_frame <- rbind(final_data_frame, rhes_papH_lower_3ss_l, higher_human_rhes_l, higher_rhes_l )  ## Join with previous clasificatin in Marmoset


#################################
## Exons originated in Chimp or Gibon
#################################



##!!!!!! Load again the previous tables before function wide to long
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
#total_furthest_chimp_gibo <- cbind(chimp_gibo_high_3ss, temp1)


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

## Rest of queries$X3SSS_nomLeu1 >= 3 | queries$X3SSS_panPan1 >= 3 | queries$X3SSS_panTro4 >= 3  ## Since its OR (could be < than 3 in other columns.. we have to deal with the rest
rest <- rest2[!(rownames(rest2) %in% rownames(chimp_gibo_high_3ss)),]
chimp_gibo_lower_3ss <- subset(rest, rest$X3SSS_nomLeu1 < 3 | rest$X3SSS_panTro4 < 3 )

higher_chimp_gibo_l <- wide_to_long(higher_chimp_gibo, "Constant Exon")
higher_human_chimp_gibo_l <- wide_to_long(higher_human_chimp_gibo, "Evolving Exon")
chimp_gibo_lower_3ss_l <- wide_to_long(chimp_gibo_lower_3ss, "Emerging Exon")


final_data_frame <- rbind(final_data_frame, chimp_gibo_lower_3ss_l, higher_human_chimp_gibo_l, higher_chimp_gibo_l)

final_data_frame <- rbind(cal_lower_3ss_l, higher_human_l,  higher_marmoset_l, chimp_gibo_lower_3ss_l, higher_human_chimp_gibo_l, higher_chimp_gibo_l, rhes_papH_lower_3ss_l, higher_human_rhes_l, higher_rhes_l )                   ## All joined


## Bootstraping

tot <- data.frame(stringsAsFactors=FALSE)
#names(tot) <- c("Median", "Normal_low", "Normal_hight", "Basic_low", "Basic_hight", "Region", "Group")

tot <- data.frame(Median=as.numeric(), 
                  Normal_low=as.numeric(),
                  Normal_hight=as.numeric(),
                  Basic_low=as.numeric(),
                  Basic_hight=as.numeric(),
                  Region=as.character(),
                  Group=as.character(),
                  stringsAsFactors=FALSE)

lista <- c("hg19", "hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3")

#lista <- c("hg19", "hg19", "panTro4") #, "nomLeu1", "rheMac3", "calJac3")

for (genome in lista){
  
  print(genome)
  toboot <- final_data_frame [ final_data_frame$region == genome,]
  
  #toboot <- final_data_frame [ final_data_frame$region == "panTro4",]
  #toboot <- toplot.panelD [ as.numeric(toplot.panelD$variable) == 2,]
  
  aggregate(toboot$WU, by=list(var=toboot$group), median)
  boot.statistics = function(data, indices){
    d <- data[ indices, ]
    d <- aggregate(d$WU, by=list(var=d$group), median)
    return(d$x)
  }
  
  #toboot$group
  
  
  bootobject1 <- boot(data=toboot, statistic = boot.statistics, R=20000,  parallel = "multicore", ncpus=4)
  
  plot(bootobject1, index=3)
  print(bootobject1)
  plot(bootobject1)
  
  emerging <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=1) ## Emerging Exon
  evolving <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=2) ## Evolving Exon
  constant <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=3) ## Constant Exon
  
  line_emerging <- c(bootobject1$t0[1], emerging$normal[2], emerging$normal[3], emerging$basic[4], emerging$basic[5], genome, "Emerging")
  line_evolving <- c(bootobject1$t0[2], evolving$normal[2], evolving$normal[3], evolving$basic[4], evolving$basic[5], genome, "Evolving")
  line_constant<- c(bootobject1$t0[3], constant$normal[2], constant$normal[3], constant$basic[4], constant$basic[5], genome, "Constant")
  

  
  tot[,1:5] <- sapply(tot[,1:5], as.numeric)
  tot[,6:7] <- sapply(tot[,6:7], as.character)
  tot <- rbind(tot, line_emerging, line_evolving, line_constant)
  names(tot) <- c("Median", "Normal_low", "Normal_hight", "Basic_low", "Basic_hight", "Region", "Group")
  
}

tot <- tot[-c(1:3),]


pdf_name <- paste0("3SS_" , experiment_name, ".pdf")
#pdf(pdf_name , width=20, height=13)
ggtitle_name <- paste0("3SS ", experiment_name)

### All i one
tot <- data.frame()
names(tot) <- c("Median", "Normal_low", "Normal_hight", "Basic_low", "Basic_hight", "Region", "Group")
lista <- c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3")



bootstrap_to_table <- function(final_data_frame, organism_subset) {
  
  organism_subset <- paste0(organism_subset)
  print(organism_subset)
  toboot <- final_data_frame [ final_data_frame$region == "hg19",]
  
}

bootstrap_to_table(final_data_frame, "hg19")

  # here you can subset by ==1, ==2, ==3, etc. if  dataframe$variable is a factor (is.factor(dataframe$variable))
  
  #toboot <- toplot.panelD [ as.numeric(toplot.panelD$variable) == 2,]
  
  aggregate(toboot$WU, by=list(var=toboot$group), median)
  boot.statistics = function(data, indices){
    d <- data[ indices, ]
    d <- aggregate(d$WU, by=list(var=d$group), median)
    return(d$x)
  }
  
  toboot$group
  
  
  bootobject1 <- boot(data=toboot, statistic = boot.statistics, R=20000,  parallel = "multicore", ncpus=4)
  plot(bootobject1, index=3)
  
  print(bootobject1)
  plot(bootobject1)
  
  emerging_hg19 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=1) ## Emerging Exon
  evolving_hg19 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=2) ## Evolving Exon
  constant_hg19 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=3) ## Constant Exon
  
  line_emerging <- c(bootobject1$t0[1], emerging_hg19$normal[2], emerging_hg19$normal[3], emerging_hg19$basic[4], emerging_hg19$basic[5], genome, "Emerging")
  line_evolving <- c(bootobject1$t0[2], evolving_hg19$normal[2], evolving_hg19$normal[3], evolving_hg19$basic[4], evolving_hg19$basic[5], genome, "Evolving")
  line_constant<- c(bootobject1$t0[3], constant_hg19$normal[2], constant_hg19$normal[3], constant_hg19$basic[4], constant_hg19$basic[5], genome, "Constant")
  

  
}



names(bootobject1)
bootobject1$statistic()
emerging_hg19[2]
emerging_hg19[6[[2]]]










### hg19



source("http://www.anthro.utah.edu/~rogers/datanal/R/scatboot.r")



toboot <- final_data_frame [ as.factor(final_data_frame$region) == "hg19",]
# here you can subset by ==1, ==2, ==3, etc. if  dataframe$variable is a factor (is.factor(dataframe$variable))

#toboot <- toplot.panelD [ as.numeric(toplot.panelD$variable) == 2,]

aggregate(toboot$WU, by=list(var=toboot$group), median)
boot.statistics = function(data, indices){
  d <- data[ indices, ]
  d <- aggregate(d$WU, by=list(var=d$group), median)
  return(d$x)
}

toboot$group


bootobject1 <- boot(data=toboot, statistic = boot.statistics, R=20000,  parallel = "multicore", ncpus=4)
plot(bootobject1, index=3)

print(bootobject1)
plot(bootobject1)

emerging_hg19 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=1) ## Emerging Exon
evolving_hg19 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=2) ## Evolving Exon
constant_hg19 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=3) ## Constant Exon

### panTro4

toboot <- final_data_frame [ as.factor(final_data_frame$region) == "panTro4",]
# here you can subset by ==1, ==2, ==3, etc. if  dataframe$variable is a factor (is.factor(dataframe$variable))

#toboot <- toplot.panelD [ as.numeric(toplot.panelD$variable) == 2,]

aggregate(toboot$WU, by=list(var=toboot$group), median)
boot.statistics = function(data, indices){
  d <- data[ indices, ]
  d <- aggregate(d$WU, by=list(var=d$group), median)
  return(d$x)
}

toboot$group


bootobject1 <- boot(data=toboot, 
                    statistic = boot.statistics, R=20000,  parallel = "multicore", ncpus=4)
plot(bootobject1, index=3)

print(bootobject1)
plot(bootobject1)

emerging_panTro4 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=1) ## Emerging Exon
evolving_panTro4 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=2) ## Evolving Exon
constant_panTro4 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=3) ## Constant Exon


### nomLeu1

toboot <- final_data_frame [ as.factor(final_data_frame$region) == "nomLeu1",]
# here you can subset by ==1, ==2, ==3, etc. if  dataframe$variable is a factor (is.factor(dataframe$variable))

#toboot <- toplot.panelD [ as.numeric(toplot.panelD$variable) == 2,]

aggregate(toboot$WU, by=list(var=toboot$group), median)
boot.statistics = function(data, indices){
  d <- data[ indices, ]
  d <- aggregate(d$WU, by=list(var=d$group), median)
  return(d$x)
}

toboot$group


bootobject1 <- boot(data=toboot, 
                    statistic = boot.statistics, R=20000,  parallel = "multicore", ncpus=4)
plot(bootobject1, index=3)

print(bootobject1)
plot(bootobject1)

emerging_nomLeu1 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=1) ## Emerging Exon
evolving_nomLeu1 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=2) ## Evolving Exon
constant_nomLeu1 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=3) ## Constant Exon



### rheMac3

toboot <- final_data_frame [ as.factor(final_data_frame$region) == "rheMac3",]
# here you can subset by ==1, ==2, ==3, etc. if  dataframe$variable is a factor (is.factor(dataframe$variable))

#toboot <- toplot.panelD [ as.numeric(toplot.panelD$variable) == 2,]

aggregate(toboot$WU, by=list(var=toboot$group), median)
boot.statistics = function(data, indices){
  d <- data[ indices, ]
  d <- aggregate(d$WU, by=list(var=d$group), median)
  return(d$x)
}

toboot$group


bootobject1 <- boot(data=toboot, 
                    statistic = boot.statistics, R=20000,  parallel = "multicore", ncpus=4)
plot(bootobject1, index=3)

print(bootobject1)
plot(bootobject1)

emerging_rheMac3<- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=1) ## Emerging Exon
evolving_rheMac3 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=2) ## Evolving Exon
constant_rheMac3 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=3) ## Constant Exon


### calJac3

toboot <- final_data_frame [ as.factor(final_data_frame$region) == "calJac3",]
# here you can subset by ==1, ==2, ==3, etc. if  dataframe$variable is a factor (is.factor(dataframe$variable))

#toboot <- toplot.panelD [ as.numeric(toplot.panelD$variable) == 2,]

aggregate(toboot$WU, by=list(var=toboot$group), median)
boot.statistics = function(data, indices){
  d <- data[ indices, ]
  d <- aggregate(d$WU, by=list(var=d$group), median)
  return(d$x)
}

toboot$group


bootobject1 <- boot(data=toboot, 
                    statistic = boot.statistics, R=20000,  parallel = "multicore", ncpus=4)
plot(bootobject1, index=3)

print(bootobject1)
plot(bootobject1)

emerging_calJac3 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=1) ## Emerging Exon
evolving_calJac3 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=2) ## Evolving Exon
constant_calJac3 <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=3) ## Constant Exon

names(bootobject1)
bootobject1$t0[1]
#### extract boot object
library("parallel") 
tot <- data.frame()
names(tot) <- c("Median", "Normal_low", "Normal_hight", "Basic_low", "Basic_hight", "Region", "Group")

names(bootobject1)
bootobject1$statistic()
emerging_hg19[2]
emerging_hg19[6[[2]]]

lista <- c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3")

line_emerging <- c(bootobject1$t0[1], emerging_hg19$normal[2], emerging_hg19$normal[3], emerging_hg19$basic[4], emerging_hg19$basic[5], genome, "Emerging")
line_evolving <- c(bootobject1$t0[2], evolving_hg19$normal[2], evolving_hg19$normal[3], evolving_hg19$basic[4], evolving_hg19$basic[5], genome, "Evolving")
line_constant<- c(bootobject1$t0[3], constant_hg19$normal[2], constant_hg19$normal[3], constant_hg19$basic[4], constant_hg19$basic[5], genome, "Constant")



tot <- rbind(tot, line_emerging, line_evolving, line_constant)

names(emerging_hg19)
emerging_hg19$basic

class(emerging_hg19)
plot(emerging_hg19)

boot <- cbind(emerging_hg19, evolving_hg19, constant_hg19)

as.character(boot[,c(4,5,6)])
as.character(boot)

toboot <- dataframe [ as.numeric(dataframe$variable) == 2,]
# here you can subset by ==1, ==2, ==3, etc. if  dataframe$variable is a factor (is.factor(dataframe$variable))

head(toboot)

#set median as bootstatistic and draw 2000 times
trap.median <- function(data, indices){
  d <- data[ indices]
  return(median(d))
}
set.seed = 1939
bootobject.test <- boot(data=toboot$value, 
                        statistic = trap.median, R=1000,  parallel = "multicore", ncpus=4)
print(bootobject.test)
boot.ci(bootobject.test, conf=0.95, type = c("norm", "basic", "BCa"))





