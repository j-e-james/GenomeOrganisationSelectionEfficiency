library(ggplot2)
library("gridExtra")
library(tidyr)
library(dplyr)
library(lme4)
library(MASS)


### Contains a lot of very repetitive code, which is unneccessary


### Removing multiple transcripts.
### Rules for picking transcript copy- of all the numbers keep the longest, if no difference, choose randomly- then use the Gene name column (which lacks the '.N' number),
### and keep the transcript number column 
### geneinfofile_ordered will then be used in combination with all downstream analysis- TranscriptID retains the '.N' number for correct matching.
geneinfofile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Protein_table_GeneRefs_and_Metrics_Athaliana_UID138.csv', header=TRUE)
nonsynfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.0fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
synfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.4fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

synfile$Pol <- rowSums(synfile[2:50])


### for synonymous and nonsynonymous
### need to do concordantly so that the same genes end up in the synonymous and nonsynonymous bins.
polfile <- merge(geneinfofile, synfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- merge(polfile, nonsynfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- polfile[order(polfile$ProteinLength, decreasing = TRUE), ] 
polfile_ordered <- polfile[!duplicated(polfile$GeneID), ]

### some sanity checks, matching gene and protein lengths:
cor.test(polfile_ordered$Ln.x+polfile_ordered$Ln.y, polfile_ordered$ProteinLength*3)
cor.test(polfile_ordered$Ln.x, polfile_ordered$ProteinLength*3)
cor.test(polfile_ordered$Ln.y, polfile_ordered$ProteinLength*3)

### This one could be argued to relate to level of pleiotropy
expressionfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/GeneExpressionGSE43858_RAW/Results_GSE43858_RAW.txt', header=TRUE)
expressionfile <- expressionfile[c(1,2,3)]
polexp <- merge(polfile_ordered, expressionfile, by.x = "GeneID", by.y = "Gene")
hist(log10(expressionfile$Mean_fpkm))


###bin by splitting genes into groups of even numbers of mutations, ranked by factor
n <- 10
polexp <- polexp[which(polexp$Mean_fpkm != 0),]
# rank by factor
exp_split <- polexp[order(polexp$Mean_fpkm),]
# identify groups by sumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polexp$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()

for (group in exp_split){
  print(mean(group$Mean_fpkm))
  print(sum(group$Pol))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$Mean_fpkm), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$Mean_fpkm), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_fpkmpolcounts_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_fpkmpolcounts_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)



#######Repeat process for length
exp_split <- polexp[order(polexp$ProteinLength),]
# identify groups by sumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polexp$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()

for (group in exp_split){
  print(mean(group$ProteinLength))
  print(sum(group$Pol))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$ProteinLength), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$ProteinLength), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_lengthpolcounts_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_lengthpolcounts_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)



#######Repeat process for SIFT
geneinfofile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Protein_table_GeneRefs_and_Metrics_Athaliana_UID138.csv', header=TRUE)
nonsynfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.0fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
synfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.4fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

synfile$Pol <- rowSums(synfile[2:50])


###first for synonymous
### no- need to do concordantly so that the same genes end up in the synonymous and nonsynonymous bins.
polfile <- merge(geneinfofile, synfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- merge(polfile, nonsynfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- polfile[order(polfile$ProteinLength, decreasing = TRUE), ] 
polfile_ordered <- polfile[!duplicated(polfile$GeneID), ]

### SIFT file
recombfile <- read.csv("/Users/jennyjames/Dropbox/1001GenomesProjectData/SIFTScore/SIFTScore_PerGene_MeanStdOnly.txt", header = TRUE)
polrecomb <- merge(polfile_ordered, recombfile, by.x = "TranscriptID", by.y = "Gene")
hist(recombfile$Mean_sift)


polexp <- polrecomb[which(polrecomb$Mean_sift != 0),]
# rank by factor
exp_split <- polexp[order(polexp$Mean_sift),]
exp_split$SumPol <- cumsum(exp_split$Pol)
# identify groups by sumulative total of Pol. This includes an n + 1 to smooth out the group size
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polexp$Pol) / n + 1)) + 1



###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)




fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()

for (group in exp_split){
  print(mean(group$Mean_sift))
  print(sum(group$Pol))
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
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_SIFTpolcounts_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_SIFTpolcounts_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)






#######Repeat process for recombination
geneinfofile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Protein_table_GeneRefs_and_Metrics_Athaliana_UID138.csv', header=TRUE)
nonsynfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.0fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
synfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.4fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

synfile$Pol <- rowSums(synfile[2:50])


###first for synonymous
### no- need to do concordantly so that the same genes end up in the synonymous and nonsynonymous bins.
polfile <- merge(geneinfofile, synfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- merge(polfile, nonsynfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- polfile[order(polfile$ProteinLength, decreasing = TRUE), ] 
polfile_ordered <- polfile[!duplicated(polfile$GeneID), ]

### This one could be argued to relate to level of pleiotropy
recombfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Recomination_doi_10.5061_dryad.v655ns36__v1/Recomb10MbResults.txt', header = TRUE)
polrecomb <- merge(polfile_ordered, recombfile, by.x = "TranscriptID", by.y = "Gene")
hist(log10(recombfile$MeanRecomb))

###bin by splitting genes into equal size groups- start with 20 groups

# rank by factor
polexp <- polrecomb 

exp_split <- polexp[order(polexp$MeanRecomb),]
# identify groups by sumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polrecomb$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)


fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()

for (group in exp_split){
  print(mean(group$MeanRecomb))
  print(sum(group$Pol))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$MeanRecomb), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$MeanRecomb), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_Recombpolcounts_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_Recombpolcounts_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)





#######Repeat process for density
geneinfofile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Protein_table_GeneRefs_and_Metrics_Athaliana_UID138.csv', header=TRUE)
nonsynfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.0fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
synfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.4fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

synfile$Pol <- rowSums(synfile[2:50])


###first for synonymous
### no- need to do concordantly so that the same genes end up in the synonymous and nonsynonymous bins.
polfile <- merge(geneinfofile, synfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- merge(polfile, nonsynfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- polfile[order(polfile$ProteinLength, decreasing = TRUE), ] 
polfile_ordered <- polfile[!duplicated(polfile$GeneID), ]

### density file
genedensityfile <- read.table('/Users/jennyjames/Dropbox/1001GenomesProjectData/Tair10Ref/Density.txt', sep = '\t', header=FALSE)
colnames(genedensityfile) <- c('Gene', 'DensityPer10kb', 'DistanceTo10Genes')
genedensityfile['DensityPerGene'] <- 10/genedensityfile$DistanceTo10Genes
recombfile <- genedensityfile
polrecomb <- merge(polfile_ordered, recombfile, by.x = "GeneID", by.y = "Gene")
hist(recombfile$DensityPer10kb)
hist(recombfile$DensityPerGene)

###bin by splitting genes into equal size groups- start with 20 groups
exp_split <- polrecomb[order(polrecomb$DensityPerGene),]

# identify groups by sumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polrecomb$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)


fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()

for (group in exp_split){
  print(mean(group$DensityPerGene))
  print(sum(group$Pol))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$DensityPerGene), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$DensityPerGene), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_DensityPerGenepolcounts_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_DensityPerGenepolcounts_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)





##### My weird factors, which I struggle to break into groups evenly by polymorphism due to their fairly discrete nature

#######Repeat process for gene age

### Removing multiple transcripts.
### Rules for picking transcript copy- of all the numbers keep the longest, if no difference, choose randomly- then use the Gene name column (which lacks the '.N' number),
### and keep the transcript number column 
### geneinfofile_ordered will then be used in combination with all downstream analysis- TranscriptID retains the '.N' number for correct matching.
geneinfofile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Protein_table_GeneRefs_and_Metrics_Athaliana_UID138.csv', header=TRUE)
nonsynfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.0fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
synfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.4fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

synfile$Pol <- rowSums(synfile[2:50])

###first for synonymous
### no- need to do concordantly so that the same genes end up in the synonymous and nonsynonymous bins.
polfile <- merge(geneinfofile, synfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- merge(polfile, nonsynfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- polfile[order(polfile$ProteinLength, decreasing = TRUE), ] 
polfile_ordered <- polfile[!duplicated(polfile$GeneID), ]

### Gene age is contained in the geneinfofile
polexp<- polfile_ordered
hist(log10(polexp$AgeOfOldestPfam_LUCA_LECA_Updated))

###with age, is not possible to bin by splitting genes into equal size groups- ages are so unevenly distributed.
### discrete groups, of all unique values
exp_split <- polexp[order(polexp$AgeOfOldestPfam_LUCA_LECA_Updated),]

exp_split<- split(exp_split, exp_split$AgeOfOldestPfam_LUCA_LECA_Updated)


fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()

for (group in exp_split){
  print(mean(group$AgeOfOldestPfam_LUCA_LECA_Updated))
  print(sum(group$Pol))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$AgeOfOldestPfam_LUCA_LECA_Updated), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$AgeOfOldestPfam_LUCA_LECA_Updated), SFS_nonsyn ) )
}
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_Age_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_Age_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)


###we will downsample groups to a set number of genes (rather than polymorphisms), to ensure that our results hold.
fpkm_syn_downsample <- data.frame()
fpkm_nonsyn_downsample <- data.frame()

for (fullgroup in exp_split){
  print(nrow(fullgroup))
  if (nrow(fullgroup) > 200){ 
    print(mean(fullgroup$AgeOfOldestPfam_LUCA_LECA_Updated))
    print(sum(fullgroup$Pol))  
    ###working on a version that also takes a sample of the largest groups.
    group <- fullgroup[sample(nrow(fullgroup), 200),]
  

    colnames(group)
    SFS_syn <- dplyr::select(group, ends_with('.x'))
    SFS_syn <- colSums(SFS_syn)
    SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
    SFS_nonsyn <- colSums(SFS_nonsyn)
    fpkm_syn_downsample <- rbind(fpkm_syn_downsample, c( mean(group$AgeOfOldestPfam_LUCA_LECA_Updated), SFS_syn ) )
    fpkm_nonsyn_downsample <- rbind(fpkm_nonsyn_downsample, c( mean(group$AgeOfOldestPfam_LUCA_LECA_Updated), SFS_nonsyn ) )
  }
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn_downsample, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_Agepoldown_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn_downsample, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_Agepoldown_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)



###we will restrict to the largest groups, and limit to 3000 polymorphisms
fpkm_syn_downsample <- data.frame()
fpkm_nonsyn_downsample <- data.frame()

for (fullgroup in exp_split){
  if (sum(fullgroup$Pol) > 2500){ 
    ###take a subset of data, based on a cumulative sum of polymorphisms
    group <- subset(fullgroup, cumsum(fullgroup$Pol) <= 3300)
    print(mean(group$AgeOfOldestPfam_LUCA_LECA_Updated))
    print(sum(group$Pol))  

    colnames(group)
    SFS_syn <- dplyr::select(group, ends_with('.x'))
    SFS_syn <- colSums(SFS_syn)
    SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
    SFS_nonsyn <- colSums(SFS_nonsyn)
    fpkm_syn_downsample <- rbind(fpkm_syn_downsample, c( mean(group$AgeOfOldestPfam_LUCA_LECA_Updated), SFS_syn ) )
    fpkm_nonsyn_downsample <- rbind(fpkm_nonsyn_downsample, c( mean(group$AgeOfOldestPfam_LUCA_LECA_Updated), SFS_nonsyn ) )
  }
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn_downsample, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_Agepolcounts_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn_downsample, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_Agepolcounts_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)






#######Repeat process for gene network connectivity groups
geneinfofile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Protein_table_GeneRefs_and_Metrics_Athaliana_UID138.csv', header=TRUE)
nonsynfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.0fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
synfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.4fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

synfile$Pol <- rowSums(synfile[2:50])

###first for synonymous
### no- need to do concordantly so that the same genes end up in the synonymous and nonsynonymous bins.
polfile <- merge(geneinfofile, synfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- merge(polfile, nonsynfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- polfile[order(polfile$ProteinLength, decreasing = TRUE), ] 
polfile_ordered <- polfile[!duplicated(polfile$GeneID), ]


### Read in file that specifies the gene connectivity groups
networkfile <- read.table('/Users/jennyjames/Dropbox/1001GenomesProjectData/GeneNetwork/genes_module_54.txt', header = TRUE, sep = '\t')
polrecomb <- merge(networkfile, polfile_ordered, by.y = "GeneID", by.x = "genes")
#### exclude genes not in a network group- group 0
polrecomb <- polrecomb[which(polrecomb$module != "0"),]
### a lot more genes are present in the lower numbered groups.
### My dataset is only the A thaliana nuclear genome- plastid genomes etc are removed, hence the loss of 3 distinct modules
hist(polrecomb$module)

###bin by module
exp_split <- polrecomb[order(polrecomb$module),]
exp_split<- split(exp_split, exp_split$module)


fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()

for (group in exp_split){
  print(paste(unique(group$module), length(group$module), sep = ' '))
  print(sum(group$Pol))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( unique(group$module)[[1]], SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( unique(group$module)[[1]], SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_GeneNetwork_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_GeneNetwork_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)



###we will downsample groups to a set number of genes (rather than polymorphisms), to ensure that our results hold.
fpkm_syn_downsample <- data.frame()
fpkm_nonsyn_downsample <- data.frame()

for (fullgroup in exp_split){
  if (nrow(fullgroup) > 200){ 
    group <- fullgroup[sample(nrow(fullgroup), 200),]
 #   print(nrow(group))
  #  print(paste(unique(group$module), length(group$module), sep = ' '))
  #  print(paste(unique(group$module)))
    print(sum(group$Pol))  
    colnames(group)
    SFS_syn <- dplyr::select(group, ends_with('.x'))
    SFS_syn <- colSums(SFS_syn)
    SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
    SFS_nonsyn <- colSums(SFS_nonsyn)
    fpkm_syn_downsample <- rbind(fpkm_syn_downsample, c( unique(group$module)[[1]], SFS_syn ) )
    fpkm_nonsyn_downsample <- rbind(fpkm_nonsyn_downsample, c( unique(group$module)[[1]], SFS_nonsyn ) )
  }
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn_downsample, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_GeneNetworkpoldown_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn_downsample, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped10_GeneNetworkpoldown_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)








geneinfofile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Protein_table_GeneRefs_and_Metrics_Athaliana_UID138.csv', header=TRUE)
nonsynfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.0fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
synfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.4fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

synfile$Pol <- rowSums(synfile[2:50])


### for synonymous and nonsynonymous
### need to do concordantly so that the same genes end up in the synonymous and nonsynonymous bins.
polfile <- merge(geneinfofile, synfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- merge(polfile, nonsynfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- polfile[order(polfile$ProteinLength, decreasing = TRUE), ] 
polfile_ordered <- polfile[!duplicated(polfile$GeneID), ]




#######Slightly different process for GO groups
GO_terms <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/GO_over400.csv')

### using polfile_ordered


colnames(GO_terms)
df_list <- list()

for (i in 1:length(colnames(GO_terms))){
  print(i)
  print(colnames(GO_terms)[i])
  gene_names <- GO_terms[,colnames(GO_terms)[i]]
  
  new_name <- paste('polfile_ordered_',colnames(GO_terms)[i], sep = '') 
  polfile_ordered_i <- polfile_ordered[polfile_ordered$GeneID %in% gene_names, ]
  if (nrow(polfile_ordered_i) > 200){
    df_list <- append(df_list, new_name)
    assign(new_name, polfile_ordered_i)
  }
  
}

### other groups (based on genomic factors) have approx 3386.95 synonymous pols.
### For GO groups: the lowest number of syn pols is 1303: downsample all the dataframes 
### to approx this number, just over a third of the other groups. we can but hope!

df_pol_restrict_list <- list()

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()

for (i in df_list){
  print(i)
  new_name <- paste('restricted_',i, sep = '')
  df <- get(i)
  print(sum(df$Pol)) ### the lowest number of syn pols is 1303: downsample all the dataframes to approx this number.
  df$cumsum <- cumsum(df$Pol)
  df_restrict <- df[which(df$cumsum <= 1350),]
  df_pol_restrict_list <- append(df_pol_restrict_list, new_name)
  assign(new_name, df_restrict)
  
  print(sum(df_restrict$Pol))
  print(mean(df_restrict$Ln.x))
  print(mean(df_restrict$Ln.y))
  GO_term <- read.table(text = i, sep = "_", as.is = TRUE)[[3]]
  
  colnames(df_restrict)
  SFS_syn <- dplyr::select(df_restrict, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(df_restrict, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c(GO_term, SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c(GO_term, SFS_nonsyn ) )  
  
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped1350_GOterms_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped1350_GOterms_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE, quote=FALSE)







