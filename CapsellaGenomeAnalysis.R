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

########################################## C. grandiflora ##########################################

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


geneinfofile <- read.table('/Users/jennyjames/Desktop/Capsella_offHarddrive/evx068_Supp/josephsetal_genedata.txt', header=TRUE)
geneinfofile$avgConnec <- geneinfofile$sumConnec/nrow(geneinfofile)
genenames <- read.table('/Users/jennyjames/Desktop/Capsella_offHarddrive/Capsella_GeneName.tsv', header = TRUE)
geneinfofile_names <- merge(genenames, geneinfofile, by.x = 'Pacid', by.y = 'PAC')


#### adding in Arabidopsis homolog information, for matching with aranet
Cap_Athal_homologs <- read.csv('/Users/jennyjames/Desktop/Capsella_offHarddrive/Crubellav1.0/download.20260318.153414/Phytozome/PhytozomeV12/Crubella/annotation/Crubella_183_v1.0.annotation_info.txt', sep = '\t', header = TRUE)
Athal_homologs <- Cap_Athal_homologs[c(1,10,11)]
geneinfofile_names <- merge(geneinfofile_names, Athal_homologs, by.x = 'Pacid', by.y = 'X.pacId')
### this is from aranet- network file (connections - the sum total of connections)
networkfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/AraNetConnectivity/AraNet_connections.txt', header = TRUE)
geneinfofile_names$Best.hit.arabi.name <- str_sub(geneinfofile_names$Best.hit.arabi.name, end = -3)
geneinfofile_names <- merge(geneinfofile_names, networkfile, by.x = "Best.hit.arabi.name", by.y = "gene", all.x=F)
### this is for A thlaiana GO group analysis.
GOfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GOTermCounts.csv')
geneinfofile_names <- merge(geneinfofile_names, GOfile, by.x = "Best.hit.arabi.name", by.y = "Gene", all.x=F)


#### Some testing for sanity checks
cor.test(geneinfofile_names$ProteinLength, geneinfofile_names$length)
plot(geneinfofile_names$ProteinLength, geneinfofile_names$length/3)


nonsynfile <- read.csv('/Users/jennyjames/Desktop/Capsella_offHarddrive/corientalis_grandiflora.SNPS.exonic.0fold.SFS.csv', header=TRUE)
synfile <- read.csv('/Users/jennyjames/Desktop/Capsella_offHarddrive/corientalis_grandiflora.SNPS.exonic.4fold.SFS.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

synfile$Pol <- rowSums(synfile[2:50])

polfile <- merge(geneinfofile_names, synfile, by.x = 'ProteinName', by.y = 'Gene')
polfile <- merge(polfile, nonsynfile, by.x = 'ProteinName', by.y = 'Gene')


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

sample_row <- polfile[1,] 
pi_s <- pi_calculate(sample_row, 'x')

pi_s <- by(polfile, seq_len(nrow(polfile)), pi_calculate, 'x')
pi_n <- by(polfile, seq_len(nrow(polfile)), pi_calculate, 'y')

polfile$'pi_S' <- pi_s
polfile$'pi_N' <- pi_n

### expression file (Mean_fpkm)
polfile$exp.means
polfile$norm.exp
### connectivity (Mean_fpkm)
polfile$sumConnec
polfile$avgConnec
### connectivityAranet
polfile$connections

### Ranking for joint expression and connectivity
polfile$ExpressionRank <- NA
polfile$ConnectivityRank <- NA
polfile$ExpressionRank[order(polfile$norm.exp)]  <- 1:nrow(polfile)
polfile$ConnectivityRank[order(polfile$connections)]  <- 1:nrow(polfile)
polfile$RankingMean <- (polfile$ConnectivityRank + polfile$ExpressionRank)/2


ggplot(polfile, aes(log10(pi_N/pi_S),log10(dnds))) +
  geom_point(aes(colour = norm.exp)) +
  theme_classic() 

ggplot(polfile, aes(log10(pi_N/pi_S),log10(dnds))) +
  geom_point(aes(colour = log10(connections))) +
  theme_classic() 




attach(polfile)

log_reduced_df <- as.data.frame(cbind(log10(polfile$pi_S), log10(polfile$pi_N/polfile$pi_S), log10(polfile$connections), polfile$norm.exp, log10(polfile$length)))
colnames(log_reduced_df) <- c('pi4', 'pi0/pi4', 'Connectivity', 'Expression', 'Length')

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



plot(norm.exp, avgConnec)
plot(norm.exp, connections, log = "y")

### 
normmod_connections_expression <- lm(avgConnec ~ norm.exp + length)
normmod_expression_connections <- lm(norm.exp ~ avgConnec + length)

normod_araconnections_expression <- lm(log10(connections) ~ norm.exp + log10(length))

summary(lm(connections ~ norm.exp + length))

expression_connecindependent <- resid(normmod_expression_connections)
connections_expressindependent <- resid(normmod_connections_expression)
aranetconnec_expressindependent <- resid(normod_araconnections_expression)

polfile$expressions_corrected <- expression_connecindependent
polfile$connections_corrected <- connections_expressindependent
polfile$aranetconnec_corrected <- aranetconnec_expressindependent

n <- 10


###expression
exp_split <- polfile[order(polfile$expressions_corrected),]
# identify groups by cumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polfile$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in exp_split){
  print(mean(group$expressions_corrected))
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
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped10_express_correc_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped10_express_correc_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)


###connectivity
exp_split <- polfile[order(polfile$connections_corrected),]
# identify groups by cumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polfile$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in exp_split){
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
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped10_connec_correc_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped10_connec_correc_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)


###Aranet connectivity
exp_split <- polfile[order(polfile$connections),]
# identify groups by cumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polfile$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in exp_split){
  print(mean(group$connections))
  print(sum(group$Pol))
  print(length(group$connections))
  mean_vals <- c(mean_vals, mean(group$connections))
  sd_vals <- c(sd_vals, sd(group$connections))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$connections), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$connections), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped15_aranetconnec_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped15_aranetconnec_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)



### GO group analysis, based on Arabi orthologs
exp_split <- polfile[order(polfile$GOTermCount),]
# identify groups by cumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polfile$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in exp_split){
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
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped10_GOTermCount_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped10_GOTermCount_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)




exp_split <- polfile[order(polfile$connections_corrected),]
# identify groups by cumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polfile$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in exp_split){
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
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped10_connec_correc_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped10_connec_correc_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)


###aranet connectivity, corrected
exp_split <- polfile[order(polfile$aranetconnec_corrected),]
# identify groups by cumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polfile$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in exp_split){
  print(mean(group$aranetconnec_corrected))
  print(sum(group$Pol))
  print(length(group$aranetconnec_corrected))
  mean_vals <- c(mean_vals, mean(group$aranetconnec_corrected))
  sd_vals <- c(sd_vals, sd(group$aranetconnec_corrected))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$aranetconnec_corrected), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$aranetconnec_corrected), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped20_aranetconnec_correc_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped20_aranetconnec_correc_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)


###ranking
exp_split <- polfile[order(polfile$RankingMean),]
# identify groups by cumulative total of Pol
exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polfile$Pol) / n)) + 1

###have generated 20 dataframes, sorted by expression
exp_split <- split(exp_split, exp_split$group)

fpkm_syn <- data.frame()
fpkm_nonsyn <- data.frame()
mean_vals <- list()
sd_vals <- list()

for (group in exp_split){
  print(mean(group$RankingMean))
  print(sum(group$Pol))
  print(length(group$RankingMean))
  mean_vals <- c(mean_vals, mean(group$RankingMean))
  sd_vals <- c(sd_vals, sd(group$RankingMean))
  print(mean(group$Ln.x))
  print(mean(group$Ln.y))
  colnames(group)
  SFS_syn <- dplyr::select(group, ends_with('.x'))
  SFS_syn <- colSums(SFS_syn)
  SFS_nonsyn <- dplyr::select(group, ends_with('.y'))
  SFS_nonsyn <- colSums(SFS_nonsyn)
  fpkm_syn <- rbind(fpkm_syn, c( mean(group$RankingMean), SFS_syn ) )
  fpkm_nonsyn <- rbind(fpkm_nonsyn, c( mean(group$RankingMean), SFS_nonsyn ) )
}

###we now need to print out the SFSs to a file, and run dadi to project the sample size down.
write.table(fpkm_syn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped20_ranking_synonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(fpkm_nonsyn, file = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Manuscript/ManuscriptRevisions/GroupingCapsella/Grouped20_ranking_nonsynonymous.csv', append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)




