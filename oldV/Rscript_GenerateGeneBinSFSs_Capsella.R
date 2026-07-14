library(ggplot2)
library("gridExtra")
library(tidyr)
library(dplyr)
library(lme4)
library(MASS)


###############################################################################################
#############                        Capsella data version                      ###############
###############################################################################################



### Removing multiple transcripts.
### Rules for picking transcript copy- of all the numbers keep the longest, if no difference, choose randomly- then use the Gene name column (which lacks the '.N' number),
### and keep the transcript number column 
### geneinfofile_ordered will then be used in combination with all downstream analysis- TranscriptID retains the '.N' number for correct matching.
geneinfofile <- read.table('/Volumes/MyPassport/Capsella_AlexM/JGIfilesCaprub1_0/Josephs_etal._Capsella_GeneData/evx068_Supp/josephsetal_genedata.txt', header=TRUE)
#geneinfofile$ProteinLength = geneinfofile$sumConnec/nrow(geneinfofile)

genenames <- read.table('/Volumes/MyPassport/Capsella_AlexM/Capsella_GeneName.tsv', header = TRUE)

geneinfofile_names <- merge(genenames, geneinfofile, by.x = 'Pacid', by.y = 'PAC')

#### Some testing for sanity checks
cor.test(geneinfofile_names$ProteinLength, geneinfofile_names$length)
plot(geneinfofile_names$ProteinLength, geneinfofile_names$length/3)

nonsynfile <- read.csv('/Volumes/MyPassport/Capsella_AlexM/NCBIfilesCaprub1_0/corientalis_grandiflora.SNPS.exonic.0fold.SFS.csv', header=TRUE)
synfile <- read.csv('/Volumes/MyPassport/Capsella_AlexM/NCBIfilesCaprub1_0/corientalis_grandiflora.SNPS.exonic.4fold.SFS.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

synfile$Pol <- rowSums(synfile[2:50])

polfile <- merge(geneinfofile_names, synfile, by.x = 'ProteinName', by.y = 'Gene')
polfile <- merge(polfile, nonsynfile, by.x = 'ProteinName', by.y = 'Gene')


###bin by splitting genes into equal size groups- start with 20 groups
### did for norm.exp, norm.exp, norm.exp
n <- 20

# rank by factor
polexp <- polfile
exp_split <- polexp[order(polexp$dnds),]
# identify groups by sumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polfile$Pol) / n)) + 1
### 1395.8 was the number of polymorphisms in A thal- in C grandiflora it's 5743.25 
sum(polfile$Pol) / n

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()

for (group in exp_split){
  print(sum(group$Pol))
  print(mean(group$dnds))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$dnds), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$dnds), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Volumes/MyPassport/Capsella_AlexM/NCBIfilesCaprub1_0/Grouped20_dnds_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Volumes/MyPassport/Capsella_AlexM/NCBIfilesCaprub1_0/Grouped20_dnds_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)






###############################################################################################
#############                         SIFT data inclusion                       ###############
###############################################################################################


### Removing multiple transcripts.
### Rules for picking transcript copy- of all the numbers keep the longest, if no difference, choose randomly- then use the Gene name column (which lacks the '.N' number),
### and keep the transcript number column 
### geneinfofile_ordered will then be used in combination with all downstream analysis- TranscriptID retains the '.N' number for correct matching.
geneinfofile <- read.table('/Volumes/MyPassport/Capsella_AlexM/JGIfilesCaprub1_0/Josephs_etal._Capsella_GeneData/evx068_Supp/josephsetal_genedata.txt', header=TRUE)
#geneinfofile$ProteinLength = geneinfofile$sumConnec/nrow(geneinfofile)

genenames <- read.table('/Volumes/MyPassport/Capsella_AlexM/Capsella_GeneName.tsv', header = TRUE)

geneinfofile_names <- merge(genenames, geneinfofile, by.x = 'Pacid', by.y = 'PAC')

nonsynfile <- read.csv('/Volumes/MyPassport/Capsella_AlexM/NCBIfilesCaprub1_0/corientalis_grandiflora.SNPS.exonic.0fold.SFS.csv', header=TRUE)
synfile <- read.csv('/Volumes/MyPassport/Capsella_AlexM/NCBIfilesCaprub1_0/corientalis_grandiflora.SNPS.exonic.4fold.SFS.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

synfile$Pol <- rowSums(synfile[2:50])

polfile <- merge(geneinfofile_names, synfile, by.x = 'ProteinName', by.y = 'Gene')
polfile <- merge(polfile, nonsynfile, by.x = 'ProteinName', by.y = 'Gene')

SIFTfile <- read.csv('/Volumes/MyPassport/Capsella_AlexM/SIFTScore_PerGene_MeanStdOnly.txt', header = TRUE)

polfile <- merge(polfile, SIFTfile, by.x = 'CarubvID', by.y = 'Gene')


###bin by splitting genes into equal size groups- start with 20 groups
### did for norm.exp, norm.exp, norm.exp
n <- 15

# rank by factor
polexp <- polfile
exp_split <- polexp[order(polexp$Mean_sift),]
# identify groups by sumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polfile$Pol) / n)) + 1
### 1395.8 was the number of polymorphisms in A thal- in C grandiflora SIFT score version it's 4180.2
sum(polfile$Pol) / n

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()

for (group in exp_split){
  print(sum(group$Pol))
  print(mean(group$Mean_sift))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$Mean_sift), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$Mean_sift), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Volumes/MyPassport/Capsella_AlexM/NCBIfilesCaprub1_0/Grouped15_MeanSIFT_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Volumes/MyPassport/Capsella_AlexM/NCBIfilesCaprub1_0/Grouped15_MeanSIFT_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)







