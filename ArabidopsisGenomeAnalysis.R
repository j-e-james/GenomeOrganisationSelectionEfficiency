library(ggplot2)
library("gridExtra")
library(tidyr)
library(dplyr)
library(lme4)
library(MASS)
library(corrplot)
library(ggfortify)
library(factoextra)
library(forecast)
library(corrplot)
library(lmerTest) # gets p-values in a straightforward way
library(gridExtra)
library(ggcorrplot)
library(ggfortify)
library(MuMIn) # gets Pseudo-R-squared values from GLMM 
library(broom)
library(car)
library(pls)
library(cowplot)
library(viridis)
library(stringr)
library(gridExtra)
library(ape)
library(patchwork)

########################################## A. thaliana ##########################################

###Useful functions for plotting

###function for printing out correlation results
cor_vals <- function(df, x, y){
  pear <- cor.test(x, y)
  spear <- cor.test(x, y, method = "spearman")
  spear_pval <- spear[3]
  pear_pval <- pear[3]
  spear_rho <- format(spear[4], digits = 2)
  pear_R <- format(pear[4], digits = 2)
  if (spear_pval < 0.05){
    spear_sig <- '*' } 
  else {
    spear_sig <- ''}
  
  if (pear_pval < 0.05){
    pear_sig <- '*' } 
  else {
    pear_sig <- ''}
  
  pear_pval <- format(pear[3], digits = 2)
  spear_pval <- format(spear[3], digits = 2)
  pear_pval <- paste(pear_pval, pear_sig, sep = '')
  spear_pval <- paste(spear_pval, spear_sig, sep = '')
  spear_string <- paste(spear_rho,spear_pval,sep = ", ")
  pear_string <- paste(pear_R,pear_pval,sep = ", ")
  return(paste(pear_string,spear_string, sep = '\n'))
}


#### function for producing linear model printed string
lm_vals <- function(df, x, y){
  m <- lm(x ~ y, df);
  
  pval <- summary(m)$coefficients[2,4]
  
  if (pval < 0.05){
    sig <- '*'
  } else {
    sig <- ''
  }
  pval <- format(summary(m)$coefficients[2,4], digits = 3)
  
  pval <- paste(pval, sig, sep = '')
  R2 <- format(summary(m)$r.squared, digits = 3)
  slope <- format(summary(m)$coefficients[2,1], scientific = TRUE, digits = 3)
  return(paste(slope, R2, pval, sep = ", "))
}




########################################## Generate polfile_ordered, without any genome biology traits ##########################################

### Removing multiple transcripts.
### Rules for picking transcript copy- of all the numbers keep the longest, if no difference, choose randomly- then use the Gene name column (which lacks the '.N' number),
### and keep the transcript number column 
### geneinfofile_ordered will then be used in combination with all downstream analysis- TranscriptID retains the '.N' number for correct matching.
geneinfofile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Protein_table_GeneRefs_and_Metrics_Athaliana_UID138.csv', header=TRUE)
nonsynfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.0fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
synfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/1001genomes_snp-short-indel_with_tair10_only_ACGTN.4fold.Poly.Bi.alignmentpolarized.50.SFS.haploid.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

### Generate a Pol column, which is Ps, number of neutral variants per gene
synfile$Pol <- rowSums(synfile[2:50])

### for synonymous and nonsynonymous: need to do concordantly so that the same genes end up in the synonymous and nonsynonymous bins.
polfile <- merge(geneinfofile, synfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- merge(polfile, nonsynfile, by.x = 'TranscriptID', by.y = 'Gene')
polfile <- polfile[order(polfile$ProteinLength, decreasing = TRUE), ] 
polfile_ordered <- polfile[!duplicated(polfile$GeneID), ]

plot(polfile_ordered$ProteinLength*3, polfile_ordered$Ln.x+polfile_ordered$Ln.y, ylab = '0 fold + 4 fold sites', xlab = 'Protein length (aa) X 3', pch = 21, bg = "gold")

pi_calculate <- function(sample, xory){
  sample_size <- 50
  L_name = paste("Ln.", xory, sep = '')
  Ln <- sample[, which(colnames(sample)==L_name)]
  SFS <- sample[, which(colnames(sample)== paste("1.", xory, sep = '')):which(colnames(sample)== paste("49.", xory, sep = ''))]
  # compute allele frequency by dividing column entries in bin by sample size
  SFS_freq <- colnames(SFS)
  SFS <- data.frame(t(SFS))
  colnames(SFS) <- c('count')
  SFS_freq <- as.numeric(sub( paste('.', xory, sep = ''), '', SFS_freq))
  # compute allele frequency by dividing column entries in bin by sample size
  SFS_freq <- SFS_freq/sample_size
  SFS['SFS_freq'] <- SFS_freq
  # compute pi per site: 2 * p * (1 - p) * bessel correction * count-per-allele-freq
  SFS['pi_site'] <- 2 * SFS$SFS_freq * (1 - SFS$SFS_freq) *
    (sample_size) / (sample_size - 1) * SFS$count
  pi <- sum(SFS$pi_site) / (sum(SFS$count) + Ln)
  return(pi)
}

sample_row <- polfile_ordered[1,] 
pi_s <- pi_calculate(sample_row, 'x')

pi_s <- by(polfile_ordered, seq_len(nrow(polfile_ordered)), pi_calculate, 'x')
pi_n <- by(polfile_ordered, seq_len(nrow(polfile_ordered)), pi_calculate, 'y')

polfile_ordered$'pi_S' <- pi_s
polfile_ordered$'pi_N' <- pi_n



### expression file (Mean_fpkm)
expressionfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/GeneExpressionGSE43858_RAW/Results_GSE43858_RAW.txt', header=TRUE)
expressionfile <- expressionfile[c(1,2,3)]
Pi_Biology_byGene <- merge(polfile_ordered, expressionfile, by.x = "GeneID", by.y = "Gene", all.x=F)

### network file (connections - the sum total of connections)
networkfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/AraNetConnectivity/AraNet_connections.txt', header = TRUE)
Pi_Biology_byGene <- merge(Pi_Biology_byGene, networkfile, by.x = "GeneID", by.y = "gene", all.x=F)
Pi_Biology_byGene$avgConnec <- Pi_Biology_byGene$connections/nrow(Pi_Biology_byGene)
 
### recombination file (rec.rate)
recombfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Recombination_Brazier2022/Recomb100KbBrazierResults.txt', header = TRUE)
Pi_Biology_byGene <- merge(Pi_Biology_byGene, recombfile, by.x = "TranscriptID", by.y = "Gene", all.x=F)

### density file (NumberOfGenesPerWindow)
Densityfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Tair10Ref/Density100KbResults.txt', header = TRUE)
Pi_Biology_byGene <- merge(Pi_Biology_byGene, Densityfile, by.x = "GeneID", by.y = "Gene", all.x=F)

GOfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GOTermCounts.csv')
Pi_Biology_byGene <- merge(Pi_Biology_byGene, GOfile, by.x = "GeneID", by.y = "Gene", all.x=F)

### ProteinLength

attach(Pi_Biology_byGene)

log_reduced_df <- as.data.frame(cbind(log10(Pi_Biology_byGene$pi_S), log10(Pi_Biology_byGene$pi_N/Pi_Biology_byGene$pi_S),
                                  log10(Pi_Biology_byGene$rec.rate), log10(Pi_Biology_byGene$Mean_fpkm),
                                  log10(Pi_Biology_byGene$NumberOfGenesPerWindow), log10(Pi_Biology_byGene$ProteinLength), log10(Pi_Biology_byGene$connections)))


colnames(log_reduced_df) <- c('pi4', 'pi0/pi4', 'RecRate', 'Meanfpkm', 'GeneDensityperWindow', 'Length', 'Connectivity')

log_reduced_df$RecRate[is.infinite(log_reduced_df$RecRate)] <- NA
log_reduced_df <- na.omit(log_reduced_df)


p.mat <- model.matrix(~0+., data= log_reduced_df) %>% 
  cor_pmat(use="pairwise.complete.obs")

cor_plot <- 
  model.matrix(~0+., data= log_reduced_df) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag=FALSE, type="upper",
             hc.order = TRUE,
             lab=TRUE, lab_size=4, 
             p.mat = p.mat,
             insig = "pch",
             colors = c("#6D9EC1", "white", "#E46726"),
             ggtheme = ggplot2::theme_bw())

attach(log_reduced_df)

summary(lm(Meanfpkm ~ Connections + Length + GeneDensityperWindow))
summary(lm(Meanfpkm ~ Connections + GeneDensityperWindow))
summary(lm(Meanfpkm ~ Connections * GeneDensityperWindow))


summary(lm(Connections ~ Meanfpkm + Length + GeneDensityperWindow))
summary(lm(Connections ~ Meanfpkm + Length))
summary(lm(Connections ~ Meanfpkm * Length))




### expression file (Mean_fpkm)
expressionfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/GeneExpressionGSE43858_RAW/Results_GSE43858_RAW.txt', header=TRUE)
expressionfile <- expressionfile[c(1,2,3)]
Pi_Biology_byGene <- merge(polfile_ordered, expressionfile, by.x = "GeneID", by.y = "Gene", all.x=F)

### network file (connections - the sum total of connections)
networkfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/AraNetConnectivity/AraNet_connections.txt', header = TRUE)
Pi_Biology_byGene <- merge(Pi_Biology_byGene, networkfile, by.x = "GeneID", by.y = "gene", all.x=F)
Pi_Biology_byGene$avgConnec <- Pi_Biology_byGene$connections/nrow(Pi_Biology_byGene)

### density file (NumberOfGenesPerWindow)
Densityfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Tair10Ref/Density100KbResults.txt', header = TRUE)
Pi_Biology_byGene <- merge(Pi_Biology_byGene, Densityfile, by.x = "GeneID", by.y = "Gene", all.x=F)


attach(Pi_Biology_byGene)

plot(connections, Mean_fpkm, log = "xy")
plot(avgConnec, Mean_fpkm, log = "xy")

logmod_connections_expression <- lm(log10(connections) ~ log10(Mean_fpkm) * log10(ProteinLength))
logmod_expression_connections <- lm(log10(Mean_fpkm) ~ log10(connections) + log10(NumberOfGenesPerWindow))

expression_connecindependent <- resid(logmod_expression_connections)
connections_expressindependent <- resid(logmod_connections_expression)
connections_expression <- lm(connections ~ Mean_fpkm)
length(resid(logmod_connections_expression))


Pi_Biology_byGene$expressions_corrected <- expression_connecindependent
Pi_Biology_byGene$connections_corrected <- connections_expressindependent


n <- 20

###expression
exp_split <- Pi_Biology_byGene[order(Pi_Biology_byGene$expressions_corrected),]
# identify groups by cumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(Pi_Biology_byGene$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in exp_split){
  print(median(group$expressions_corrected))
  print(sum(group$Pol))
  print(length(group$expressions_corrected))
  mean_vals <- c(mean_vals, mean(group$expressions_corrected))
  sd_vals <- c(sd_vals, sd(group$expressions_corrected))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$expressions_corrected), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$expressions_corrected), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/Grouping/Grouped_express_correc_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/Grouping/Grouped_express_correc_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)



###connectivity
con_split <- Pi_Biology_byGene[order(Pi_Biology_byGene$connections_corrected),]
# identify groups by cumulative total of Pol
con_split$group <- cumsum(con_split$Pol) %/% (ceiling(sum(Pi_Biology_byGene$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
con_split <- split(con_split, con_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in con_split){
  print(mean(group$connections_corrected))
  print(sum(group$Pol))
  print(length(group$connections_corrected))
  mean_vals <- c(mean_vals, mean(group$connections_corrected))
  sd_vals <- c(sd_vals, sd(group$connections_corrected))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$connections_corrected), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$connections_corrected), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/Grouping/Grouped_connec_correc_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/Grouping/Grouped_connec_correc_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)




###GEAR expression
Gear <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/MBQN_20250715.csv')

Gear_simplified <- Gear[4:length(Gear),]

Gear_numeric <- Gear_simplified %>% mutate(across(-c(GEAR.ID), as.numeric) )

Gear_numeric <- data.frame(ID=Gear_numeric[,1], GearMeans=rowMeans(Gear_numeric[,-1], na.rm = TRUE))

Pi_Biology_byGene <- merge(Pi_Biology_byGene, Gear_numeric, by.x = "GeneID", by.y = "ID", all.x=T)


Gear_split <- Pi_Biology_byGene[order(Pi_Biology_byGene$GearMeans),]
# identify groups by cumulative total of Pol
Gear_split$group <- cumsum(Gear_split$Pol) %/% (ceiling(sum(Pi_Biology_byGene$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
Gear_split <- split(Gear_split, Gear_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in Gear_split){
  print(mean(group$GearMeans))
  print(sum(group$Pol))
  print(length(group$GearMeans))
  mean_vals <- c(mean_vals, mean(group$GearMeans))
  sd_vals <- c(sd_vals, sd(group$GearMeans))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$GearMeans), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$GearMeans), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/Grouping/Grouped_GearMeans_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/Grouping/Grouped_GearMeans_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)



###GO counts
con_split <- Pi_Biology_byGene[order(Pi_Biology_byGene$GOTermCount),]
# identify groups by cumulative total of Pol
con_split$group <- cumsum(con_split$Pol) %/% (ceiling(sum(Pi_Biology_byGene$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
con_split <- split(con_split, con_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in con_split){
  print(mean(group$GOTermCount))
  print(sum(group$Pol))
  print(length(group$GOTermCount))
  mean_vals <- c(mean_vals, mean(group$GOTermCount))
  sd_vals <- c(sd_vals, sd(group$GOTermCount))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$GOTermCount), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$GOTermCount), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/Grouping/Grouped5_GOTermCount_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/Grouping/Grouped5_GOTermCount_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)




###Expression coeff var
Pi_Biology_byGene <- Pi_Biology_byGene[which(Pi_Biology_byGene$Mean_fpkm != 0),]
Pi_Biology_byGene <- Pi_Biology_byGene[which(Pi_Biology_byGene$Std_fpkm != 0),]

con_split <- Pi_Biology_byGene[order(Pi_Biology_byGene$Std_fpkm/Pi_Biology_byGene$Mean_fpkm),]
# identify groups by cumulative total of Pol
con_split$group <- cumsum(con_split$Pol) %/% (ceiling(sum(Pi_Biology_byGene$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
con_split <- split(con_split, con_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in con_split){
  print(mean(group$Std_fpkm/group$Mean_fpkm))
  print(sum(group$Pol))
  print(length(group$Std_fpkm/group$Mean_fpkm))
  mean_vals <- c(mean_vals, mean(group$Std_fpkm/group$Mean_fpkm))
  sd_vals <- c(sd_vals, sd(group$Std_fpkm/group$Mean_fpkm))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$Std_fpkm/group$Mean_fpkm), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$Std_fpkm/group$Mean_fpkm), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/Grouping/Grouped20_CoeffVarfpkm_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/Grouping/Grouped20_CoeffVarfpkm_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)




