#!/usr/bin/env Rscript

## Created on 10th of may  2016
## @author: Igor Ruiz de los Mozos


#################################
## Classification of Alu elements by divergence or evolutionary dynamics, and analysis of their 3’ splice sites and U-tracts
#################################
#
# 1. First we get the most distant specie in wehere we could find homologous Alu sequences inserted.
#
# 2. Then we classified the Alu elements based on the evolutionary dynamics of their 3’ss.
#   ‘Emerging’ Alu exons have a 3’ss with a score less than 3 in the most distant species.
#   ‘Stable’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength increased towards human by less than 1.
#   ‘Evolving’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength in human is more than 1+(score in the distant species). For example, if the score in marmoset is 2.5, then the Alu exon is considered as ‘emerging’, if it is 4 in marmoset and 4.5 in human, then it’s considered as ‘stable’, and if it’s 4 in marmoset and 6 in human, then it’s considered as ‘evolving’.
#
# 3. Plot 3 splice site strenght, Longest U track on Alu, longest U track on left arm, longest U track on right arm
#
#


## Load data table from data preprocesing step
read.table(total, file = "/Results/whole_final_5sp.WIDE.RICH.tab.txt", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)


#################################
## Get the most distant species to which we could lift over its genome coordinates
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
  lista_n <- lista[!is.na(lista)]                             # Remove any NA (if NA is present is because Alu exon could not been lift over on that specie)
  temp1$furthest_all[i] <- as.character(list(lista))
  temp1$furthest[i] <- as.character(tail(lista_n, n=1))       # Last element on the list is the FURTHEST specie to been lift over

}
head(temp1)
class(temp1)

## Add FURTHEST coloumn info on table
names(temp1) <- c("aluexon", "furthest_all", "furthest")
total_furthest <- merge.data.frame(total, temp1, by="aluexon", all=TRUE)
total_furthest <- cbind(total, temp1)


#################################
## Classification Alu elements based on the evolutionary dynamics of their 3’ss.
#################################
#
# 2. Then we classified the Alu elements based on the evolutionary dynamics of their 3’ss.
#   ‘Emerging’ Alu exons have a 3’ss with a score less than 3 in the most distant species.
#   ‘Stable’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength increased towards human by less than 1.
#   ‘Evolving’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength in human is more than 1+(score in the distant species). For example, if the score in marmoset is 2.5, then the Alu exon is considered as ‘emerging’, if it is 4 in marmoset and 4.5 in human, then it’s considered as ‘stable’, and if it’s 4 in marmoset and 6 in human, then it’s considered as ‘evolving’.
#


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

### Two groups to plot higher_marmoset and higher_human

higher_marmoset_l <- wide_to_long(higher_marmoset, "Constant Exon")
higher_human_l <- wide_to_long(higher_human, "Evolving Exon")
cal_lower_3ss_l <- wide_to_long(cal_lower_3ss, "Emerging Exon")

final_data_frame <- rbind(cal_lower_3ss_l, higher_human_l,  higher_marmoset_l)

exons_originated_marmoset <- final_data_frame

write.table(exons_originated_marmoset, file = "exons_originated_marmoset.tab", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)





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

higher_rhes_l <- wide_to_long(higher_rhes, "Constant Exon")
higher_human_rhes_l <- wide_to_long(higher_human_rhes, "Evolving Exon")
rhes_papH_lower_3ss_l <- wide_to_long(rhes_papH_lower_3ss, "Emerging Exon")


final_data_frame <- rbind(final_data_frame, rhes_papH_lower_3ss_l, higher_human_rhes_l, higher_rhes_l )  ## Join with previous clasificatin in Marmoset


exons_originated_rhesus <- final_data_frame

write.table(exons_originated_rhesus, file = "exons_originated_rhesus.tab", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)





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

## Rest of queries$X3SSS_nomLeu1 >= 3 | queries$X3SSS_panPan1 >= 3 | queries$X3SSS_panTro4 >= 3  ## Since its OR (could be < than 3 in other columns.. we have to deal with the rest
rest <- rest2[!(rownames(rest2) %in% rownames(chimp_gibo_high_3ss)),]
chimp_gibo_lower_3ss <- subset(rest, rest$X3SSS_nomLeu1 < 3 | rest$X3SSS_panTro4 < 3 )

higher_chimp_gibo_l <- wide_to_long(higher_chimp_gibo, "Constant Exon")
higher_human_chimp_gibo_l <- wide_to_long(higher_human_chimp_gibo, "Evolving Exon")
chimp_gibo_lower_3ss_l <- wide_to_long(chimp_gibo_lower_3ss, "Emerging Exon")


final_data_frame <- rbind(final_data_frame, chimp_gibo_lower_3ss_l, higher_human_chimp_gibo_l, higher_chimp_gibo_l)



exons_originated_chimp_gibbon <- final_data_frame

write.table(exons_originated_chimp_gibbon, file = "exons_originated_chimp_gibbon.tab", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)



######################################################
######## Plot Fig 5 A B & Fig5S B C D
######################################################
#
# 3. Plot 3 splice site strenght, Longest U track on Alu, longest U track on left arm, longest U track on right arm
#

###############

setwd('./Results' )

## Set name of the experiment to "Exons Marmoset Originated" Alu exons that have been originated in
experiment_name <- "Exons Marmoset Originated"

## Call Funcion passing the data_frame to plot
plot_custom_bar_plots(exons_originated_marmoset, experiment_name)

## Now with Alu exons orginated in rhesus
experiment_name <- "Exons rhesus Originated"
plot_custom_bar_plots(exons_originated_rhesus, experiment_name)

## Now with Alu exons orginated in chimp or gibbon
experiment_name <- "Exons chimp_gibbon Originated"
plot_custom_bar_plots(exons_originated_chimp_gibbon, experiment_name)



#####################################   Functions   #############################################################

## In this scipt I have used two funcions
#  1. WIDE to LONG allow to extract any column from WIDE table, factorise the column/elemnt studied and return a LONG table ready to plot
#  2. Function that is feed with dataframe in LONG format of AluID 3SSS UTRack_whole UtrackLeft UtraclRight
# Call Funcion passing the data_frame to plot
#


############################################
## Function creates a LONG table (need to plot) from the WIDE
############################################
#
# WIDE to LONG allow to extract any column from WIDE table, factorise the column/elemnt studied and return a LONG table ready to plot
#

wide_to_long <- function(data_frame, new_column){
  
  # Function transform the wide table of lift ovet to long table in where each row has a factor that represent its procedence and used for later plotting
  
  colnames(data_frame) <- NULL
  hg19_q <- data.frame(data_frame[,1], data_frame[,14:24])
  pantro4_q <- data.frame(data_frame[,1], data_frame[,25:35])
  nomLeu1_q<-  data.frame(data_frame[,1], data_frame[,36:46])
  rheMac3_q <- data.frame(data_frame[,1], data_frame[,47:57])
  calJac3_q <- data.frame(data_frame[,1], data_frame[,58:68])

  out_data_frame <- rbind(hg19_q, pantro4_q, nomLeu1_q, rheMac3_q, calJac3_q)
  
  colnames(out_data_frame) <- c("aluexon", "chr", "start", "end", "position", "strand", "distance_to_alu", "X3SSS", "WU", "U1", "U2", "region")# , "group")  
  
  out_data_frame$group <- as.factor(new_column)
  out_data_frame <- out_data_frame[complete.cases(out_data_frame),]
  
  return(out_data_frame)
  
}


############################################
## Plots for X3SSS WU U1 and U2
############################################

# Function that is feed with dataframe in LONG format of AluID 3SSS UTRack_whole UtrackLeft UtraclRight
# Usage:
#         plot_custom_bar_plots(final_data_frame, experiment_name)
#
# Example:
#         experiment_name <- "Exons Marmoset Originated"
#         plot_custom_bar_plots(exons_originated_marmoset, experiment_name)

## Function to plot 3SSS WU U1 and U2 at the same time ## Now ussing median  
plot_custom_bar_plots <- function(final_data_frame, experiment_name) {
  
  #### 3SSS## Plots
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  #Mean
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = mean(X3SSS, na.rm = TRUE))
  #Mode
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = getmode(X3SSS, na.rm = TRUE))
  #Median
  p_meds <- ddply(final_data_frame, .(region, group), summarise, med = median(X3SSS, na.rm = TRUE))
  p_meds
  
  a<-  ggplot(final_data_frame, aes(factor(region), X3SSS), alpha = 1, colour = "black") + 
    geom_boxplot() + 
    #scale_y_log10() + 
    theme_bw() +
    ggtitle("median quartile 25%, 50%, 75%") + 
    xlab("") + 
    ylab("length (nt)") +
    theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"))
  
  a 
  
  count_table <- data.frame(table(final_data_frame$region, final_data_frame$group))
  colnames(count_table) <- c("region", "group", "freq")
  count_table
  
  

  ### quartile <- median.quartile(clusters$length)     Nejc what its clusters
  
  dodge <- position_dodge(width = 0.9)
  
  pdf_name <- paste0("3SS_" , experiment_name, ".pdf")
  #pdf(pdf_name , width=20, height=13)
  ggtitle_name <- paste0("3SS ", experiment_name)
  
  plot3SS<- ggplot(final_data_frame, aes(factor(region), X3SSS, fill = group), alpha = 1, width = 0.5) + 
    geom_violin(position = dodge, alpha = 0.6) + 
    geom_boxplot(width=0.3, position = dodge, alpha = 1, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
    theme_bw() +
    #theme(panel.background = element_rect(fill='white', colour='black')) +
    scale_y_continuous(limits = c(-10 , 16)) +
    #scale_y_log10() +
    ggtitle(ggtitle_name) + 
    xlab("") + 
    ylab("3SS score") +
    #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
    scale_fill_manual(values=c("#F42837","#36B1BF","#2C4259")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
    geom_text(data = p_meds, aes(x = region, y = med, label = format((med+1), digits=2)), size = 3, vjust = -1.5, position = dodge, color="white") +
    geom_text(data = count_table, aes(x = region, y = 15 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"), labels=c("hg19"="Human", "panTro4"="Chimp", "nomLeu1"="Gibbon", "rheMac3"="Rhesus", "calJac3"="Marmoset"))
  
  
  plot3SS
  
  ggsave(pdf_name, width=20, height=13)
  
  #### U track Whole aluexon WU ##
  
  #Mean
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = mean(WU, na.rm = TRUE))
  #Mode
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = getmode(WU, na.rm = TRUE))
  #Median
  p_meds <- ddply(final_data_frame, .(region, group), summarise, med = median(WU, na.rm = TRUE))
  p_meds
  
  a<-  ggplot(final_data_frame, aes(factor(region), WU), alpha = 1, colour = "black") + 
    geom_boxplot() + 
    #scale_y_log10() + 
    theme_bw() +
    ggtitle("median quartile 25%, 50%, 75%") + 
    xlab("") + 
    ylab("length (nt)") +
    theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"))
  
  a 
  
  count_table <- data.frame(table(final_data_frame$region, final_data_frame$group))
  colnames(count_table) <- c("region", "group", "freq")
  count_table
  
  
  library(plyr)
  ### quartile <- median.quartile(clusters$length)     Nejc what its clusters
  
  dodge <- position_dodge(width = 0.9)
  
  pdf_name <- paste0("WU_" , experiment_name, ".pdf")
  ggtitle_name <- paste0("WU ", experiment_name)
  
  plotWU<- ggplot(final_data_frame, aes(factor(region), WU, fill = group), alpha = 1, width = 0.5) + 
    geom_violin(position = dodge, alpha = 0.6) + 
    geom_boxplot(width=0.3, position = dodge, alpha = 0.6, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
    theme_bw() +
    #theme(panel.background = element_rect(fill='white', colour='black')) +
    scale_y_continuous(limits = c(0 , 30)) +
    #scale_y_log10() +
    ggtitle(ggtitle_name) + 
    xlab("") + 
    ylab("U track lenght") +
    #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
    scale_fill_manual(values=c("#F42837","#36B1BF","#2C4259")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
    geom_text(data = p_meds, aes(x = region, y = med, label = format((med+1), digits=2)), size = 3, vjust = -1.5, position = dodge, color="white") +
    geom_text(data = count_table, aes(x = region, y = 0 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"), labels=c("hg19"="Human", "panTro4"="Chimp", "nomLeu1"="Gibbon", "rheMac3"="Rhesus", "calJac3"="Marmoset"))
  
  
  plotWU
  
  ggsave(pdf_name, width=20, height=13)
  
  #### U1 track aluexon  ##
  
  #Mean
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = mean(U1, na.rm = TRUE))
  #Mode
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = getmode(U1, na.rm = TRUE))
  #Median
  p_meds <- ddply(final_data_frame, .(region, group), summarise, med = median(U1, na.rm = TRUE))
  p_meds
  
  a<-  ggplot(final_data_frame, aes(factor(region), U1), alpha = 1, colour = "black") + 
    geom_boxplot() + 
    #scale_y_log10() + 
    theme_bw() +
    ggtitle("median quartile 25%, 50%, 75%") + 
    xlab("") + 
    ylab("length (nt)") +
    theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"))
  
  a 
  
  count_table <- data.frame(table(final_data_frame$region, final_data_frame$group))
  colnames(count_table) <- c("region", "group", "freq")
  count_table
  
  
  library(plyr)
  ### quartile <- median.quartile(clusters$length)     Nejc what its clusters
  
  dodge <- position_dodge(width = 0.9)
  
  pdf_name <- paste0("U1_" , experiment_name, ".pdf")
  ggtitle_name <- paste0("U1 ", experiment_name)
  
  plotU1<- ggplot(final_data_frame, aes(factor(region), U1, fill = group), alpha = 1, width = 0.5) + 
    geom_violin(position = dodge, alpha = 0.6) + 
    geom_boxplot(width=0.3, position = dodge, alpha = 1, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
    theme_bw() +
    #theme(panel.background = element_rect(fill='white', colour='black')) +
    scale_y_continuous(limits = c(0 , 30)) +
    #scale_y_log10() +
    ggtitle(ggtitle_name) + 
    xlab("") + 
    ylab("U track lenght") +
    #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
    scale_fill_manual(values=c("#F42837","#36B1BF","#2C4259")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
    geom_text(data = p_meds, aes(x = region, y = med, label = format((med+1), digits=2)), size = 2, vjust = -1.5, position = dodge, color="white") +
    geom_text(data = count_table, aes(x = region, y = 0 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"), labels=c("hg19"="Human", "panTro4"="Chimp", "nomLeu1"="Gibbon", "rheMac3"="Rhesus", "calJac3"="Marmoset"))
  
  
  plotU1
  
  ggsave(pdf_name, width=20, height=13)
  
  
  #### U2 track aluexon  ##
  
  #Mean
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = mean(U2, na.rm = TRUE))
  #Mode
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = getmode(U2, na.rm = TRUE))
  #Median
  p_meds <- ddply(final_data_frame, .(region, group), summarise, med = median(U2, na.rm = TRUE))
  p_meds
  
  
  a<-  ggplot(final_data_frame, aes(factor(region), U2), alpha = 1, colour = "black") + 
    geom_boxplot() + 
    #scale_y_log10() + 
    theme_bw() +
    ggtitle("median quartile 25%, 50%, 75%") + 
    xlab("") + 
    ylab("length (nt)") +
    theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"))
  
  a 
  
  count_table <- data.frame(table(final_data_frame$region, final_data_frame$group))
  colnames(count_table) <- c("region", "group", "freq")
  count_table
  
  
  library(plyr)
  ### quartile <- median.quartile(clusters$length)     Nejc what its clusters
  
  dodge <- position_dodge(width = 0.9)
  
  pdf_name <- paste0("U2_" , experiment_name, ".pdf")
  ggtitle_name <- paste0("U2 ", experiment_name)
  
  plotU2<- ggplot(final_data_frame, aes(factor(region), U2, fill = group), alpha = 1, width = 0.5) + 
    geom_violin(position = dodge, alpha = 1, alpha = 0.6) + 
    geom_boxplot(width=0.3, position = dodge, alpha = 1, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
    theme_bw() +
    #theme(panel.background = element_rect(fill='white', colour='black')) +
    scale_y_continuous(limits = c(0 , 30)) +
    #scale_y_log10() +
    ggtitle(ggtitle_name) + 
    xlab("") + 
    ylab("U track lenght") +
    #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
    scale_fill_manual(values=c("#F42837","#36B1BF","#2C4259")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
    geom_text(data = p_meds, aes(x = region, y = med, label = format((med+1), digits=2)), size = 2, vjust = -1.5, position = dodge, color="white") +
    geom_text(data = count_table, aes(x = region, y = 0 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"), labels=c("hg19"="Human", "panTro4"="Chimp", "nomLeu1"="Gibbon", "rheMac3"="Rhesus", "calJac3"="Marmoset"))
  
  
  plotU2
  
  ggsave(pdf_name, width=20, height=13)
  
  
  
  dev.off()
} ## End function

