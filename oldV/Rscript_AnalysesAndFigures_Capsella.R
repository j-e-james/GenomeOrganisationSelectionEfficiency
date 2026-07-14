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



########################################## Generate polfile_ordered, without any genome biology traits ##########################################
nonsynfile <- read.csv('/Volumes/MyPassport/Capsella_AlexM/NCBIfilesCaprub1_0/corientalis_grandiflora.SNPS.exonic.0fold.SFS.csv', header=TRUE)
synfile <- read.csv('/Volumes/MyPassport/Capsella_AlexM/NCBIfilesCaprub1_0/corientalis_grandiflora.SNPS.exonic.4fold.SFS.csv', header=TRUE)
colnames(synfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')
colnames(nonsynfile) <- c('Gene', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,'Ln')

synfile$Pol <- rowSums(synfile[2:50])

### for synonymous and nonsynonymous
### need to do concordantly so that the same genes end up in the synonymous and nonsynonymous bins.
polfile <- merge(synfile, nonsynfile, by = 'Gene')
#polfile <- polfile[order(polfile$ProteinLength, decreasing = TRUE), ] 
polfile_ordered <- polfile[!duplicated(polfile$Gene), ]



#plot(polfile_ordered$ProteinLength*3, polfile_ordered$Ln.x+polfile_ordered$Ln.y, ylab = '0 fold + 4 fold sites', xlab = 'Protein length (aa) X 3', pch = 21, bg = "gold")

pi_calculate <- function(sample, xory){
  sample_size <- 50
  L_name = paste("Ln.", xory, sep = '')
  Ln <- sample[, which(colnames(sample)==L_name)]
  SFS <- sample[, which(colnames(sample)== paste("1.", xory, sep = '')):which(colnames(sample)== paste("49.", xory, sep = ''))]
  # compute allele frequency by deviding column entries in bin by sample size
  SFS_freq <- colnames(SFS)
  SFS <- data.frame(t(SFS))
  colnames(SFS) <- c('count')
  SFS_freq <- as.numeric(sub( paste('.', xory, sep = ''), '', SFS_freq))
  # compute allele frequency by deviding column entries in bin by sample size
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


########################################## Supplementary figure 3 ##########################################

p4_plot <- ggplot(polfile_ordered, aes(pi_S)) + 
  geom_histogram(colour="black", fill="orange", bins=20) + 
  theme_bw(base_size=16) + labs(y = expression(italic('Capsella grandiflora'))) +  
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

geneinfofile <- read.table('/Volumes/MyPassport/Capsella_AlexM/JGIfilesCaprub1_0/Josephs_etal._Capsella_GeneData/evx068_Supp/josephsetal_genedata.txt', header=TRUE)
geneinfofile$avgConnec = geneinfofile$sumConnec/nrow(geneinfofile)


genenames <- read.table('/Volumes/MyPassport/Capsella_AlexM/Capsella_GeneName.tsv', header = TRUE)

geneinfofile_names <- merge(genenames, geneinfofile, by.x = 'Pacid', by.y = 'PAC')

Pi_Biology_byGene <- merge(polfile_ordered, geneinfofile_names, by.x = "Gene", by.y = "ProteinName", all.x=T)


###add new merge to include SIFT scores.
SIFTfile <- read.csv('/Volumes/MyPassport/Capsella_AlexM/SIFTScore_PerGene_MeanStdOnly.txt', header = TRUE)
names(SIFTfile)[names(SIFTfile) == 'Mean_sift'] <- 'MeanSIFT'
Pi_Biology_byGene <- merge(Pi_Biology_byGene, SIFTfile, by.x = 'CarubvID', by.y = 'Gene', all.x=T)



reduced_df <- as.data.frame(cbind(log10(Pi_Biology_byGene$pi_S), log10(Pi_Biology_byGene$pi_N/Pi_Biology_byGene$pi_S),
                                    Pi_Biology_byGene$ProteinLength, 
                                    Pi_Biology_byGene$ds, 
                                    Pi_Biology_byGene$norm.exp, 
                                    Pi_Biology_byGene$avgConnec, 
                                    Pi_Biology_byGene$MeanSIFT ))


colnames(reduced_df) <- c('pi4', 'pi0pi4', 'Length', 'dS', 'expression', 'connectivity', 'SIFT')


p.mat <- model.matrix(~0+., data= reduced_df) %>% 
  cor_pmat(use="pairwise.complete.obs")


########################################## Supplementary figure 4 ##########################################

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
Property_pca <- prcomp(na.omit(reduced_df)[,3:7], scale = TRUE)

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

set.seed(1)
model <- pcr(pi4 ~ Length+dS+expression+connectivity+SIFT, data = reduced_df, scale=TRUE, validation="CV")
summary(model)
### suggests that 3 principle components is the correct number to use, with little additional advantage of adding further components.
validationplot(model)
loadings(model)
### performing a principle component regression
set.seed(1)
model <- pcr(pi0pi4 ~ Length+dS+expression+connectivity+SIFT, data = reduced_df, scale=TRUE, validation="CV")
summary(model)
### suggests that 3 principle components is the correct number to use, with little additional advantage of adding further components.
validationplot(model)
loadings(model)



########################################## Load DFE file data - full variables, 20 groups. Generates pop_statistics ##########################################

source('/Applications/polyDFE-master/postprocessing.R')

#### generating a data frame for polyDFE results for every genome feature considered
pop_statistics <- data.frame()

features_list <- c('ProteinLength', 'dnds', 'norm.exp', 'avgConnec', 'MeanSIFT')

for (feature in features_list){
  path = '/Users/jennyjames/Desktop/Capsella_offHarddrive/PolyDFE/'
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
    print(barplot(dfe_for_plot, ylim = c(0, 0.9), col = viridis(7)[4] ))
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
    limit_alpha <- sapply(est, function(e) c("supLimit = 0 " = estimateAlpha(e), "supLimit = 1 " = estimateAlpha(e, supLimit = 1), "supLimit = 5 " = estimateAlpha(e, supLimit = 5), "supLimit = 10 " = estimateAlpha(e, supLimit = 10)))
    alpha_val <- limit_alpha[3]
    
    prime_values <- data.frame(feature, file, pattern_value, dfe, b_val, S_d_val, p_b_val, S_b_val, alpha_val, eps_an_val, criteria_val)
    colnames(prime_values) <- c('GroupingCategory', 'Filename', 'GroupValue', "<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1", 'b', 'S_d', 'p_b', 'S_b', 'alpha', 'eps_an', 'criteria')
    pop_statistics <- rbind(prime_values, pop_statistics)
    
  }
}


#### For figures, labels are taken from these 'pretty names'- however, check that values are labelled as log transformed correctly
names <- pop_statistics$GroupingCategory
pretty_names <- list()

for (name in names){
  if (name == 'norm.exp'){
    pretty_names <- append(pretty_names, 'Expression level')}
  if (name == 'ProteinLength'){ 
    pretty_names <- append(pretty_names, 'Length')}  
  if (name == 'avgConnec'){
    pretty_names <- append(pretty_names, 'Connectivity')} 
  if (name == 'dnds'){
    pretty_names <- append(pretty_names, 'dn/ds')}
  if (name == 'MeanSIFT'){
    pretty_names <- append(pretty_names, 'SIFT score')}    
}


pop_statistics$pretty_names <- pretty_names

pop_statistics$Factor <- sapply(str_split(pop_statistics$Filename, '_'), `[`, 1)
#numbers <- gregexpr("[0-9]+", pop_statistics$Factor)
#result <- regmatches(pop_statistics$Factor, numbers)
#numeric_result <- as.numeric(unlist(result))
#pop_statistics$Factor <- numeric_result



########################################## Perform transform and update values in pop_statistics ##########################################

#### Do not do this- perform transforms later

#### write in log transform function- for features- apart from expression
#pop_statistics[pop_statistics$GroupingCategory == 'ProteinLength', 'GroupValue'] <- log10(pop_statistics[pop_statistics$GroupingCategory == 'ProteinLength', 'GroupValue'])
#pop_statistics[pop_statistics$GroupingCategory == 'avgConnec', 'GroupValue'] <- log10(pop_statistics[pop_statistics$GroupingCategory == 'avgConnec', 'GroupValue'])
#pop_statistics[pop_statistics$GroupingCategory == 'MeanSIFT', 'GroupValue'] <- log10(pop_statistics[pop_statistics$GroupingCategory == 'MeanSIFT', 'GroupValue'])
#pop_statistics[pop_statistics$GroupingCategory == 'dnds', 'GroupValue'] <- log10(pop_statistics[pop_statistics$GroupingCategory == 'dnds', 'GroupValue'])



########################################## Calculate pi over PolyDFE input files, and merge down into pop_statistics ##########################################


#### calculate pi using polyDFE input
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
  path = '/Users/jennyjames/Desktop/Capsella_offHarddrive/PolyDFE/'
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

pop_statistics <- merge(pop_statistics, variable_df)



########################################## Generate mutation_groups- mutation rates per group  ##########################################

mutation_groups <- data.frame()

variable_list <- c('ProteinLength', 'dnds', 'norm.exp', 'avgConnec', 'MeanSIFT')
pop_statistics_list <- c('ProteinLength', 'dnds', 'norm.exp', 'avgConnec', 'MeanSIFT')

for (var_index in 1:length(variable_list)) {
  i <- variable_list[var_index]
  print(i)
  print(summary(Pi_Biology_byGene[[i]]))
  variable <- pop_statistics_list[var_index]
  
  n <- 20
  
  polexp <- Pi_Biology_byGene[which(Pi_Biology_byGene[[i]] != 'NA'),]

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
    len_var <- (length(df$Gene))
    subset_list <- as.list(df$Gene)
    mutation_groups <- rbind(mutation_groups, c(paste(variable, count, sep = ''), SD, mean(df[[i]]), len_var, mean(df$ds), mean(df$dnds)))
    count <- count + 1
    
  }
}

colnames(mutation_groups) <- c('Factor', 'SD', 'FactorMean', 'NumGenes', 'MutationRate', 'dnds')

### Mutation rate here is approximated as ds, but we need to account for the number of generations since the separation from A thaliana, which was how these 
### estimates were calculated, to get the per generation rate- taken from Timetree.org, as 9.9MYA
mutation_groups$NumGenes <- as.numeric(mutation_groups$NumGenes)
mutation_groups$MutationRate <- as.numeric(mutation_groups$MutationRate)
mutation_groups$dnds <- as.numeric(mutation_groups$dnds)
mutation_groups$MutationRate <- mutation_groups$MutationRate/(9.9*10^6)


########################################## Generate df_overall- combine DFE data (pop_statistics) with mutation rate (mutation_groups) estimated per gene group ########################################## 

df_overall <- merge(pop_statistics, mutation_groups, by = "Factor")
df_overall$'Ne' <- df_overall$pi4 / (4 * df_overall$MutationRate)

df_to_write <- apply(df_overall,2,as.character)
write.csv(df_to_write, "/Users/jennyjames/Desktop/Capsella_offHarddrive/PolyDFE/df_overall_DFEGenomeBioMutRates.csv", row.names = FALSE)

########################################## Create overall figures ########################################## 

df_overall <- read.csv("/Users/jennyjames/Desktop/Capsella_offHarddrive/PolyDFE/df_overall_DFEGenomeBioMutRates.csv", header = TRUE,check.names=FALSE)

##### adding a constant value to the negative 'norm.exp' values, so that they can be log transformed
df_overall[df_overall$GroupingCategory == 'norm.exp', 'GroupValue'] <- df_overall[df_overall$GroupingCategory == 'norm.exp', 'GroupValue']+15


### splitting by GroupingCategory
df_overall_split <- split(df_overall, df_overall$GroupingCategory)
                       
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
    theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    coord_cartesian(ylim = c(0,1)) + guides(fill=guide_legend(title=pretty_name))  + theme(legend.position = c(0.85, 0.55))  

  
  p_plot <- ggplot(factor, aes(factor$GroupValue, factor$pi04)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#BFFFAD", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, factor$pi04,log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(pi[0]/pi[4]))+ scale_x_log10() + ylim(-0.1, 0.5)
  
  p4_plot <- ggplot(factor, aes(factor$GroupValue, factor$pi4)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#79E378", size = 2)+
    annotate("text",  x = I(0), y = I(0), label = lm_vals(factor, factor$pi4 ,log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(pi[4])) + scale_x_log10() + ylim(0.005, 0.019)
  
  S_plot <- ggplot(factor, aes(factor$GroupValue, factor$S_d*-1/factor$Ne)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#45425A", size = 2)+
    annotate("text",  x = I(0), y = I(0), label = lm_vals(factor, factor$S_d*-1/factor$Ne,log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(italic(s))) + scale_x_log10() + ylim(-0.003, 0.01)
  
  b_plot <- ggplot(factor, aes(factor$GroupValue, factor$b)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#464A45", size = 2)+
    theme_classic() + labs(x = pretty_name, y = expression(italic(beta))) + 
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor,factor$b,log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    scale_x_log10() + ylim(0, 0.65)
  
  
  se_plot <- ggplot(factor, aes(factor$GroupValue, factor$S_d*-1)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#BFFFAD", size = 2)+
    geom_errorbar(aes(xmin = factor$GroupValue-as.numeric(factor$SD), xmax = factor$GroupValue+as.numeric(factor$SD))) +
    annotate("text",  x = I(0), y = I(0), label = lm_vals(factor, factor$S_d,log10(factor$GroupValue)) , vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = 'standard deviation') + scale_x_log10()
  
  Sd_plot <- ggplot(factor, aes(factor$GroupValue, factor$S_d*-1)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#BFFFAD", size = 2)+
    #    geom_errorbar(aes(xmin = factor$GroupValue-as.numeric(factor$CI95), xmax = factor$GroupValue+as.numeric(factor$CI95))) +
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
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, factor$alpha, log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(italic(alpha)))+ scale_x_log10()
  
  dnds_plot <- ggplot(factor, aes(factor$GroupValue, factor$dnds)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#464A45", size = 2)+
    annotate("text", x = -Inf, y = -Inf, label = lm_vals(factor, factor$dnds, log10(factor$GroupValue)), vjust = -1, hjust = -0.1) + 
    theme_classic() + labs(x = pretty_name, y = expression(italic(d[N])/italic(d[S])) )+ scale_x_log10()
 

  alpha_ds <- ggplot(factor, aes(factor$GroupValue, factor$alpha*factor$dnds)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "#464A45", size = 2)+
    annotate("text", x = -Inf, y = -Inf, label = lm_vals(factor, (factor$alpha*factor$dnds), log10(factor$GroupValue)), vjust = -1, hjust = -0.1) + 
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


























