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

########################################## Supplementary figure 1 ##########################################


p4_plot <- ggplot(polfile_ordered, aes(pi_S)) + 
  geom_histogram(colour="black", fill="orange", bins=20) + 
  theme_bw(base_size=16) + labs(y = expression(italic('Arabidopsis thaliana'))) +  
  scale_x_log10(expression(pi[4])) +
  annotate("text", label = paste('Median = ', signif(median(polfile_ordered$pi_S), digits = 3)), x= I(0.2), y = I(0.8))+
  coord_cartesian(xlim = c(1e-05, 0.2))

p_plot <- ggplot(polfile_ordered, aes(pi_N/pi_S)) + 
  geom_histogram(colour="black", fill="midnightblue", bins=20) + 
  theme_bw(base_size=16,) + labs(y = '') + 
  scale_x_log10(expression(pi[0]/pi[4])) +
  annotate("text", label = paste('Median = ', signif(median(polfile_ordered$pi_N/polfile_ordered$pi_S), digits = 3)), x= I(0.2), y = I(0.8))+
  coord_cartesian(xlim = c(1e-04, 1e03))

plot_grid(p4_plot, p_plot)


########################################## Create Pi_Biology_byGene, merging polfile_ordered with genome biology traits ##########################################

### dnds file
### data from https://academic.oup.com/genetics/article/224/2/iyad074/7140272#490983955
dndsfile <- read.csv('/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/MilesJosephs2023_dnds.csv', header = TRUE)
dndsfile <- dndsfile[c(1,20,21)]
dndsfile$dnds <- dndsfile$dN/dndsfile$dS
Pi_Biology_byGene <- merge(polfile_ordered, dndsfile, by.x = "GeneID", by.y = "name", all.x=T)


### expression file
expressionfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/GeneExpressionGSE43858_RAW/Results_GSE43858_RAW.txt', header=TRUE)
expressionfile <- expressionfile[c(1,2,3)]
Pi_Biology_byGene <- merge(Pi_Biology_byGene, expressionfile, by.x = "GeneID", by.y = "Gene", all.x=T)

### SIFT file
Siftfile <- read.csv("/Users/jennyjames/Dropbox/1001GenomesProjectData/SIFTScore/SIFTScore_PerGene_MeanStdOnly.txt", header = TRUE)
Pi_Biology_byGene <- merge(Pi_Biology_byGene, Siftfile, by.x = "TranscriptID", by.y = "Gene", all.x=T)

### length
### just polfile_ordered$ProteinLength)

Densityfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Tair10Ref/Density100KbResults.txt', header = TRUE)
### collumn headers are 'PercentCDS' and 'NumberOfGenesPerWindow'
Pi_Biology_byGene <- merge(Pi_Biology_byGene, Densityfile, by.x = "GeneID", by.y = "Gene", all.x=T)

networkfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/AraNetConnectivity/AraNet_connections.txt', header = TRUE)
Pi_Biology_byGene <- merge(Pi_Biology_byGene, networkfile, by.x = "GeneID", by.y = "gene", all.x=T)

### recombination file
recombfile <- read.csv('/Users/jennyjames/Dropbox/1001GenomesProjectData/Recombination_Brazier2022/Recomb100KbBrazierResults.txt', header = TRUE)
Pi_Biology_byGene <- merge(Pi_Biology_byGene, recombfile, by.x = "TranscriptID", by.y = "Gene", all.x=T)


reduced_df <- as.data.frame(cbind(log10(Pi_Biology_byGene$pi_S), log10(Pi_Biology_byGene$pi_N/Pi_Biology_byGene$pi_S),
                                  Pi_Biology_byGene$rec.rate, Pi_Biology_byGene$Mean_fpkm, Pi_Biology_byGene$Mean_sift,
                                  Pi_Biology_byGene$NumberOfGenesPerWindow, Pi_Biology_byGene$ProteinLength, Pi_Biology_byGene$connections, Pi_Biology_byGene$dS))

colnames(reduced_df) <- c('pi4', 'pi0pi4', 'RecRate', 'Meanfpkm', 'Mean_sift', 'GeneDensityperWindow', 'Length', 'Connections', 'dS')



p.mat <- model.matrix(~0+., data= reduced_df) %>% 
  cor_pmat(use="pairwise.complete.obs")



########################################## Supplementary figure 2 ##########################################

cor_plot <- 
  model.matrix(~0+., data= reduced_df) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag=FALSE, type="upper",
             hc.order = TRUE,
             lab=TRUE, lab_size=4, 
             p.mat = p.mat,
             insig = "pch",
             colors = c("#6D9EC1", "white", "#E46726"),
             ggtheme = ggplot2::theme_bw())




### just including the predictor variables in a PCA
Property_pca <- prcomp(na.omit(reduced_df)[,3:9], scale = TRUE)

print(Property_pca)
summary(Property_pca)
get_eigenvalue(Property_pca)
fviz_eig(Property_pca, col.var="blue")

fviz_plot <- fviz_pca_var(Property_pca, 
                          col.var = "cos2",
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                          repel = TRUE # Avoid text overlapping
)


plot_grid(cor_plot, fviz_plot, labels = c("A", "B"))

########################################## Principle component regression ##########################################

### performing a principle component regression
set.seed(10)
model <- pcr(pi4 ~ RecRate + Meanfpkm+ Mean_sift+ GeneDensityperWindow+ Length+ Connections+ dS, data = reduced_df, scale=TRUE, validation="CV")
summary(model)
### suggests that 3 principle components is the correct number to use.
validationplot(model)
loadings(model)


### performing a principle component regression
set.seed(10)
model <- pcr(pi0pi4 ~ RecRate + Meanfpkm+ Mean_sift+ GeneDensityperWindow+ Length+ Connections+ dS, data = reduced_df, scale=TRUE, validation="CV")
summary(model)
### suggests that 2 principle components is teh correct number to use
validationplot(model)
loadings(model)



########################################## Load DFE file data - full variables, 20 groups. Generates pop_statistics_rename ##########################################


source('/Applications/polyDFE-master/postprocessing.R')

#### generating a data frame for polyDFE results for every genome feature considered
pop_statistics <- data.frame()

features_list <- c('DensityGeneWindowpolcounts', 'lengthpolcounts', 'RecombBrazierpolcounts', 'SIFTpolcounts', 
                   'NetworkConnectivity20polcounts', 'fpkmpolcounts', 'GOterms')


for (feature in features_list){
  path = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/PolyDFE/'
  ###ideally we would sort these files
  files <- list.files(path = path, pattern = paste(feature, '.+\\.output', sep = ''))
  print(files)
  pattern_df <- as.data.frame(do.call(rbind, str_split(files, pattern = '_')))
  files[with(pattern_df, order(V2, V3))]
  
  for (file in files){
    print(file)
    pattern <- str_split(file, pattern = '_')[[1]][2]
    pattern_value <- as.numeric(pattern)
    
    est <-parseOutput(paste(path, file, sep = ""))
    dfe_for_plot <- getDiscretizedDFE(est[[1]], c(-100, -10, -1, 0, 1))
    #    print(barplot(dfe_for_plot, ylim = c(0, 0.9), col = viridis(7)[4] ))
    text(5, 0.8, labels = pattern_value)
    
    dfe <- data.frame(getDiscretizedDFE(est[[1]], c(-100, -10, -1, 0, 1)))[1,]
    values <- est[[1]]$values[[1]]
    b_val <- values["b"]
    S_d_val <- values["S_d"]
    p_b_val <- values["p_b"]
    S_b_val <- values["S_b"]
    criteria_val <- est[[1]]$criteria
    eps_an_val <- values["eps_an"]
    
    alpha_val<- est[[1]]$alpha[[1]]    
    
    prime_values <- data.frame(feature, file, pattern_value, dfe, b_val, S_d_val, p_b_val, S_b_val, alpha_val, eps_an_val, criteria_val)
    colnames(prime_values) <- c('GroupingCategory', 'Filename', 'GroupValue', "<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1", 'b', 'S_d', 'p_b', 'S_b', 'alpha', 'eps_an', 'criteria')
    pop_statistics <- rbind(prime_values, pop_statistics)
    
  }
}


#### some hard-coded elements here for GO term analyses. These tend to be run separately from the continuous features
#### N.B. check the labelling error for GO terms due to rounding.
pop_statistics_rename <- pop_statistics
GO_names <- pop_statistics_rename$Filename[pop_statistics_rename$GroupingCategory == 'GOterms']
GO_names <-sub("GOterms[0-9]+_GO.", "", GO_names)
GO_names <-sub("_PolyDFE.output", "", GO_names)
pop_statistics_GO <- pop_statistics_rename[which(pop_statistics_rename$GroupingCategory == 'GOterms'),]
pop_statistics_GO$GroupValue <- as.numeric(GO_names)
pop_statistics_GO$GroupValue <- c(40387,37947,1184030,101062,179616,1181057,63685,52937,101896,22191,2813,1119,56094,30278,135084,121123)
pop_statistics_rest <- pop_statistics_rename[which(pop_statistics_rename$GroupingCategory != 'GOterms'),]
pop_statistics_rename <- rbind(pop_statistics_rest, pop_statistics_GO)


#### For figures, labels are taken from these 'pretty names'
names <- pop_statistics_rename$GroupingCategory
pretty_names <- list()

for (name in names){
  if (name == 'DensityGeneWindowpolcounts'){
    pretty_names <- append(pretty_names, 'Gene density per 100kb')}
  if (name == 'DensityPercentCDSpolcounts'){
    pretty_names <- append(pretty_names, '%CDS per 100kb')}  
  if (name == 'fpkmpolcounts'){
    pretty_names <- append(pretty_names, 'Expression level')} 
  if (name == 'RecombBrazierpolcounts'){
    pretty_names <- append(pretty_names, 'Recombination rate')} 
  if (name == 'SIFTpolcounts'){
    pretty_names <- append(pretty_names, 'SIFT score')}
  if (name == 'GOterms'){
    pretty_names <- append(pretty_names, 'GO group size')}
  if (name == 'NetworkConnectivity20polcounts'){ 
    pretty_names <- append(pretty_names, 'Connectivity')}
  if (name == 'lengthpolcounts'){ 
    pretty_names <- append(pretty_names, 'Length')}
}

pop_statistics_rename$pretty_names <- pretty_names
pop_statistics_rename$Factor <- sapply(str_split(pop_statistics_rename$Filename, '_'), `[`, 1)



########################################## Perform transform and update values in pop_statistics_rename ##########################################

#### Do not do this- perform transforms later

#### write in log transform function- for features
#pop_statistics_rename[pop_statistics_rename$GroupingCategory == 'fpkmpolcounts', 'GroupValue'] <- log10(pop_statistics_rename[pop_statistics_rename$GroupingCategory == 'fpkmpolcounts', 'GroupValue'])
#pop_statistics_rename[pop_statistics_rename$GroupingCategory == 'NetworkConnectivity20polcounts', 'GroupValue'] <- log10(pop_statistics_rename[pop_statistics_rename$GroupingCategory == 'NetworkConnectivity20polcounts', 'GroupValue'])
#pop_statistics_rename[pop_statistics_rename$GroupingCategory == 'lengthpolcounts', 'GroupValue'] <- log10(pop_statistics_rename[pop_statistics_rename$GroupingCategory == 'lengthpolcounts', 'GroupValue'])
#pop_statistics_rename[pop_statistics_rename$GroupingCategory == 'GOterms', 'GroupValue'] <- log10(pop_statistics_rename[pop_statistics_rename$GroupingCategory == 'GOterms', 'GroupValue'])





########################################## Calculate pi over PolyDFE input files, and merge down into pop_statistics_rename ##########################################

pi_calculate <- function(SFS, L, sample_size){
  colnames(SFS) <- as.numeric(str_replace(colnames(SFSsel), "V", ""))
  SFS_freq <- as.numeric(str_replace(colnames(SFSsel), "V", ""))
  
  # compute allele frequency by dividing column entries in bin by sample size
  SFS_freq <- SFS_freq/sample_size
  SFS <- data.frame(t(SFS))
  colnames(SFS) <- c('count')
  SFS['SFS_freq'] <- SFS_freq
  # compute pi per site: 2 * p * (1 - p) * bessel correction * count-per-allele-freq
  SFS['pi_site'] <- 2 * SFS$SFS_freq * (1 - SFS$SFS_freq) *
    (sample_size) / (sample_size - 1) * SFS$count
  #  print(SFS)
  pi <- sum(SFS$pi_site) / (sum(SFS$count) + L)
  return(pi)  
}

variable_df <- data.frame()
for (variable in features_list){
  ###input file rows are ordered neut, then sel 
  path = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/PolyDFE/'
  files <- list.files(path = path, pattern = paste(variable, '.+\\.input', sep = ''))
  pattern_df <- as.data.frame(do.call(rbind, str_split(files, pattern = '_')))
  pi0 <- list()
  pi4 <- list()
  pi04 <- list() 
  for (file in files){
    SFS <- read.table(paste(path,file, sep = ""), header = FALSE, skip = 1)
    SFSneut <- SFS[1,][-length(SFS)]
    Lneut <- SFS[1,length(SFS)]
    SFSsel <-  SFS[2,][-length(SFS)]
    Lsel <- SFS[2,length(SFS)] 
    sample <- length(SFS)
    
    pi4_item <- pi_calculate(SFSneut, Lneut, sample)
    pi0_item <- pi_calculate(SFSsel, Lsel, sample)
    pi4 <- append(pi4, pi4_item) 
    pi0 <- append(pi0, pi0_item) 
    pi04 <- append(pi04, pi0_item/pi4_item)
  }
  pattern_df <- as.data.frame(pattern_df)
  colnames(pattern_df) <- c('Factor', 'Variable', 'source')
  pattern_df$pi4 <- as.numeric(pi4)
  pattern_df$pi0 <- as.numeric(pi0)
  pattern_df$pi04 <- as.numeric(pi04)
  pattern_df$Variable <- as.numeric(pattern_df$Variable)
  pattern_df$GroupingVariable <- variable
  
  variable_df <- rbind(variable_df, pattern_df)  
}


variable_df_rename <- variable_df
variable_df_GO <- variable_df_rename[which(variable_df_rename$GroupingVariable == 'GOterms'),]
variable_df_GO$Variable <- c(40387,37947,1184030,101062,179616,1181057,63685,52937,101896,22191,2813,1119,56094,30278,135084,121123)
variable_df_rest <- variable_df_rename[which(variable_df_rename$GroupingVariable != 'GOterms'),]
variable_df_rename <- rbind(variable_df_rest, variable_df_GO)

pop_statistics_rename <- merge(pop_statistics_rename, variable_df_rename)



########################################## Generate mutation_groups- mutation rates per group  ##########################################

mutation_groups <- data.frame()

#### continuous variable options
variable_list <-  c('rec.rate', 'Mean_fpkm', 'Mean_sift', 'NumberOfGenesPerWindow', 'ProteinLength', 'connections', 'dS')

#### match up names in the continuous variables used with factors in pop_statistics_rename
variable_list <-  c('rec.rate', 'Mean_fpkm', 'Mean_sift', 'NumberOfGenesPerWindow', 'ProteinLength', 'connections')
pop_statistics_list <- c('RecombBrazierpolcounts', 'fpkmpolcounts', 'SIFTpolcounts', 'DensityGeneWindowpolcounts', 'lengthpolcounts', 'NetworkConnectivity20polcounts')


for (var_index in 1:length(variable_list)) {
  i <- variable_list[var_index]
  print(i)
  print(summary(Pi_Biology_byGene[[i]]))
  variable <- pop_statistics_list[var_index]
  
  n <- 20
  
  polexp <- Pi_Biology_byGene[which(Pi_Biology_byGene[[i]] != 'NA'),]
  polexp <- polexp[which(polexp[[i]] > 0),]
  
  # rank by factor
  exp_split <- polexp[order(polexp[[i]]),]
  # identify groups by cumulative total of Pol
  exp_split$group <- cumsum(exp_split$Pol) %/% (ceiling(sum(polexp$Pol) / n)) + 1
  
  ###have generated 20 dataframes, sorted by variable
  exp_split <- split(exp_split, exp_split$group)
  
  count = 1
  for (df in exp_split){
    print(paste(variable, count, sep = ''))
    print(summary(df[[i]]))
    SD <- sd(df[[i]])
    len_var <- (length(df$GeneID))
    subset_list <- as.list(df$GeneID)
    mutation_groups <- rbind(mutation_groups, c(paste(variable, count, sep = ''), SD, mean(df[[i]]), len_var, mean(df$dS, na.rm = TRUE), mean(df$dN, na.rm = TRUE)))
    count <- count + 1
    
  }
}

colnames(mutation_groups) <- c('Factor', 'SD', 'FactorMean', 'NumGenes', 'dS', 'dN')


#### discrete variables - GO terms
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
### to approx this number, just over a third of the other groups.
GO_term_names <- list()
df_pol_restrict_list <- list()
count = 1
for (i in df_list){
  print(i)
  new_name <- paste('GOterms',count, sep = '')
  df <- get(i)
  print(sum(df$Pol)) ### the lowest number of syn pols is 1303: downsample all the dataframes to approx this number.
  df$cumsum <- cumsum(df$Pol)
  df_restrict <- df[which(df$cumsum <= 1350),]
  
  df_pol_restrict_list <- append(df_pol_restrict_list, new_name)
  
  subset_list <- as.list(df_restrict$GeneID)
  
  muts_subset <- subset(Pi_Biology_byGene, GeneID %in% subset_list)
  print(mean(muts_subset$dS,na.rm = TRUE))
  
  len_var <- length(muts_subset$dS)
  print(new_name)
  
  
  mutation_groups <- rbind(mutation_groups, c(new_name, 0, 0, len_var, mean(muts_subset$dS, na.rm = TRUE), mean(muts_subset$dN, na.rm = TRUE)))
  
  
  GO_term <- read.table(text = i, sep = "_", as.is = TRUE)[[3]]
  GO_term_names <- append(GO_term_names, c(new_name, GO_term))
  
  count = count + 1
}



########################################## Generate df_overall- combine DFE data (pop_statistics_rename) with mutation rate (mutation_groups) estimated per gene group ########################################## 

df_overall <- merge(pop_statistics_rename, mutation_groups, by = "Factor")

df_overall$dS <- as.numeric(df_overall$dS)
df_overall$dN <- as.numeric(df_overall$dN)
df_overall$SD <- as.numeric(df_overall$SD)
df_overall$NumGenes <- as.numeric(df_overall$NumGenes)

### Mutation rate here is approximated as ds, but we need to account for the number of generations since the separation from A thaliana, which was how these 
### estimates were calculated, to get the per generation rate- taken from Timetree.org, as 6.2MYA
df_overall$MutationRate <- df_overall$dS/(6.2*10^6)

df_overall$'Ne' <- df_overall$pi4 / (4 * df_overall$MutationRate)

df_to_write <- apply(df_overall,2,as.character)
write.csv(df_to_write, "/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/df_overall_DFEGenomeBiodSMutRates.csv", row.names = FALSE)

########################################## Create overall figures ########################################## 

df_overall <- read.csv("/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/df_overall_DFEGenomeBiodSMutRates.csv", header = TRUE, check.names=FALSE)

### splitting by GroupingCategory
df_overall_split <- split(df_overall, df_overall$GroupingCategory)

print(range(df_overall$pi4))
print(range(df_overall$pi04))
print(range(df_overall$S_d*-1/df_overall$Ne))
print(range(df_overall$b))


count = 0
for (factor in df_overall_split){
  count = count + 1
  print(count)
  factor_name <- names(df_overall_split)[count]
  print(factor_name)
  factor_name <-sub("pol.*", "", factor_name)
  pretty_name <- unique(factor$pretty_name)
  
  print(summary(lm(factor$pi04 ~ factor$GroupValue)) )
  print(summary(lm(factor$pi4 ~ factor$GroupValue)) )
  print(summary(lm(factor$Ne ~ factor$GroupValue)) )
  
  disc_dfe <- factor[c(4:10)]
  
  df1<- cbind('GroupValue' = disc_dfe$GroupValue, '<-100' = disc_dfe$`<-100`, 'categories' = "<-100") 
  df2<- cbind('GroupValue' = disc_dfe$GroupValue, '-100 to -10' = disc_dfe$`-100 to -10`, 'categories' = "-100 to -10") 
  df3<- cbind('GroupValue' = disc_dfe$GroupValue, '-10 to -1' = disc_dfe$`-10 to -1`, 'categories' = "-10 to -1") 
  df4<- cbind('GroupValue' = disc_dfe$GroupValue, '-1 to 0' = disc_dfe$`-1 to 0`, 'categories' = "-1 to 0") 
  df5<- cbind('GroupValue' = disc_dfe$GroupValue, '0 to 1' = disc_dfe$`0 to 1`, 'categories' = "0 to 1") 
  df6<- cbind('GroupValue' = disc_dfe$GroupValue, '>1' = disc_dfe$`>1`, 'categories' = ">1") 
  
  disc_dfe <- rbind(df1, df2, df3, df4, df5, df6)
  
  colnames(disc_dfe) <- c('GroupValue', 'fraction', 'categories')
  disc_dfe <- data.frame(disc_dfe)
  disc_dfe[2] <- sapply(disc_dfe[,2], as.numeric) 
  disc_dfe <- as.data.frame(lapply(disc_dfe, unlist))
  
  disc_dfe$categories <- factor(disc_dfe$categories, levels = c("<-100", "-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  disc_dfe$GroupValue <- disc_dfe$GroupValue
  
  disc_dfe$GroupValue <- signif(as.numeric(disc_dfe$GroupValue), digits = 3)
  disc_dfe$GroupValue <- factor(disc_dfe$GroupValue, sort(as.numeric(unique(disc_dfe$GroupValue))) )
  
  
  disc_dfe <- ggplot(disc_dfe, aes(fill = factor(GroupValue), x=categories, y=fraction)) +
    scale_fill_manual(values = viridis(length(unique(disc_dfe$GroupValue)))) +
    geom_bar(position="dodge", stat="identity") +
    theme_classic() + labs(y = "Fraction of mutations", x = expression(paste("Fitness effects (", italic("N"['e']), italic("s"), ")" ) )  ) +
    coord_cartesian(ylim = c(0,1)) + guides(fill=guide_legend(title=pretty_name))  + theme(legend.position = c(0.85, 0.55))  
  
  
  p_plot <- ggplot(factor, aes(factor$GroupValue, factor$pi04)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#BFFFAD", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, factor$pi04,log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(pi[0]/pi[4]))+ scale_x_log10() + ylim(0, 0.44)
  
  p4_plot <- ggplot(factor, aes(factor$GroupValue, factor$pi4)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#79E378", size = 2)+
    annotate("text",  x = I(0), y = I(0), label = lm_vals(factor, factor$pi4 ,log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(pi[4])) + scale_x_log10() + ylim(0.0025, 0.01)
  
  S_plot <- ggplot(factor, aes(factor$GroupValue, factor$S_d*-1/factor$Ne)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#45425A", size = 2)+
    annotate("text",  x = I(0), y = I(0), label = lm_vals(factor, factor$S_d*-1/factor$Ne,log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(italic(s))) + scale_x_log10() + ylim(-0.01, 0.025)
  
  b_plot <- ggplot(factor, aes(factor$GroupValue, factor$b)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#464A45", size = 2)+
    theme_classic() + labs(x = pretty_name, y = expression(italic(beta))) + 
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor,factor$b,log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    scale_x_log10() + ylim(0.15, 0.42)
  
  
  se_plot <- ggplot(factor, aes(factor$GroupValue, factor$S_d*-1)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#BFFFAD", size = 2)+
    geom_errorbar(aes(xmin = factor$GroupValue-as.numeric(factor$SD), xmax = factor$GroupValue+as.numeric(factor$SD))) +
    annotate("text",  x = I(0), y = I(0), label = lm_vals(factor, factor$S_d,log10(factor$GroupValue)) , vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = 'standard deviation') + scale_x_log10()
  
  Sd_plot <- ggplot(factor, aes(factor$GroupValue, factor$S_d*-1)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#BFFFAD", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, factor$S_d,log10(factor$GroupValue)) , vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(italic(Sd))) + scale_x_log10() 
  
  Ne_plot <- ggplot(factor, aes(factor$GroupValue, factor$Ne)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#6C7D47", size = 2)+
    annotate("text",  x = I(0), y = I(0), label = lm_vals(factor, factor$Ne,log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(italic(Ne)))  + scale_x_log10()
  
  alpha_plot <- ggplot(factor, aes(factor$GroupValue, factor$alpha)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#464A45", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, factor$alpha,log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(italic(alpha)))+ scale_x_log10()
 
  dnds_plot <- ggplot(factor, aes(factor$GroupValue, factor$dN/factor$dS)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#464A45", size = 2)+
    annotate("text", x = -Inf, y = -Inf, label = lm_vals(factor, factor$dN/factor$dS, log10(factor$GroupValue)), vjust = -1, hjust = -0.1) + 
    theme_classic() + labs(x = pretty_name, y = expression(italic(d[N])/italic(d[S])) ) + scale_x_log10()
   
  alpha_ds <- ggplot(factor, aes(factor$GroupValue, factor$alpha*factor$dN/factor$dS)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#464A45", size = 2)+
    annotate("text", x = -Inf, y = -Inf, label = lm_vals(factor, factor$alpha*factor$dN/factor$dS, log10(factor$GroupValue)), vjust = -1, hjust = -0.1) + 
    theme_classic() + labs(x = pretty_name, y = expression(italic(omega[alpha]))) + scale_x_log10()

  
  print(plot_grid(alpha_plot, dnds_plot, alpha_ds, labels = c("A", "B", "C"), nrow = 1, ncol = 3))  
  #print(grid.arrange(se_plot, p4_plot, Ne_plot, Sd_plot, S_plot, p_plot, b_plot, alpha_plot, ncol = 2, nrow = 4))
  print(disc_dfe)
  print(Ne_plot)
  
  print(p4_plot / 
          p_plot / 
          S_plot / 
          b_plot)
  
  
}





