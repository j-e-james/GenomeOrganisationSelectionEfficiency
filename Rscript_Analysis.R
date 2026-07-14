
####Stats- analysing polyDFE output, generating plots

################# requires postprocessing.R
library(ggplot2)
library(viridis)
library(stringr)
library(gridExtra)
library(dplyr)
library(ape)
library(cowplot)
library(lme4)

source('/Applications/polyDFE-master/postprocessing.R')

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


####basing this on population-specific loops I have generated previously
pop_statistics <- data.frame()
### despite its name, this just needs to be a variable that matches filenames

### A thaliana
species_list <- c('CoeffVarFpkm', 'conneccorrec', 'expresscorrec', 'length', 'NetworkConnectivity', 'fpkm', 'RecombBrazier')
species_list <- c('conneccorrec', 'expresscorrec', 'NetworkConnectivity', 'fpkm')
species_list <- c('GOTermCount')
species_list <- c('NetworkConnectivity', 'fpkm')


### C grandiflora
species_list <- c('conneccorrec', 'expresscorrec', 'length', 'NetworkConnectivity', 'NormExpression')
species_list <- c('GOTermCount')
species_list <- c('NetworkConnectivity', 'NormExpression')
species_list <- c('conneccorrec', 'expresscorrec', 'NetworkConnectivity', 'NormExpression')


for (species in species_list){
  path = '/Users/jennyjames/Desktop/CrossGenomeDFE/GroupingCapsella/PolyDFE10/'
  files <- list.files(path = path, pattern = paste(species, '.+\\.output', sep = ''))
  print(files)
  pattern_df <- as.data.frame(do.call(rbind, str_split(files, pattern = '_')))
  files[with(pattern_df, order(V2, V3))]
  
  for (file in files){
    print(file)
    pattern <- str_split(file, pattern = '_')[[1]][2]
    pattern_value <- as.numeric(pattern)
    print(pattern_value)
    
    est <-parseOutput(paste(path, file, sep = ""))
    dfe_for_plot <- getDiscretizedDFE(est[[1]], c(-100, -10, -1, 0, 1))
    #    print(barplot(dfe_for_plot, ylim = c(0, 0.9), col = viridis(7)[4] ))
    #    text(5, 0.8, labels = pattern_value)
    
    dfe <- data.frame(getDiscretizedDFE(est[[1]], c(-100, -10, -1, 0, 1)))[1,]
    values <- est[[1]]$values[[1]]
    b_val <- values["b"]
    S_d_val <- values["S_d"]
    p_b_val <- values["p_b"]
    S_b_val <- values["S_b"]
    criteria_val <- est[[1]]$criteria
    eps_an_val <- values["eps_an"]
    
    alpha_val<- est[[1]]$alpha[[1]]    
    
    prime_values <- data.frame(species, file, pattern_value, dfe, b_val, S_d_val, p_b_val, S_b_val, alpha_val, eps_an_val, criteria_val)
    colnames(prime_values) <- c('GroupingCategory', 'Filename', 'GroupValue', "<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1", 'b', 'S_d', 'p_b', 'S_b', 'alpha', 'eps_an', 'criteria')
    pop_statistics <- rbind(prime_values, pop_statistics)
    
  }
}

pop_statistics_rename <- pop_statistics
names <- pop_statistics_rename$GroupingCategory
pretty_names <- list()

for (name in names){
  if (name == 'expresscorrec'){
    pretty_names <- append(pretty_names, 'Independent expression level')}
  if (name == 'conneccorrec'){ 
    pretty_names <- append(pretty_names, 'Independent connectivity')}
  if (name == 'length'){ 
    pretty_names <- append(pretty_names, 'length')}
  if (name == 'fpkm'){ 
    pretty_names <- append(pretty_names, 'Expression level (fpkm)')}
  if (name == 'RecombBrazier'){
    pretty_names <- append(pretty_names, 'Recombination rate')}
  if (name == 'NetworkConnectivity'){
    pretty_names <- append(pretty_names, 'Gene network connectivity')}
  if (name == 'GOTermCount'){ 
    pretty_names <- append(pretty_names, 'Number of associated GO terms')}
  if (name == 'CoeffVarFpkm'){ 
    pretty_names <- append(pretty_names, 'CV, expression level (fpkm)')}
  if (name == 'NormExpression'){ 
    pretty_names <- append(pretty_names, 'Expression level, normalised counts')}
  
}


pop_statistics_rename$pretty_names <- pretty_names

pop_statistics_rename$Factor <- sapply(str_split(pop_statistics_rename$Filename, '_'), `[`, 1)




### splitting by GroupingCategory
pop_statistics_split <- split(pop_statistics_rename, pop_statistics_rename$GroupingCategory)

count = 0
for (factor in pop_statistics_split){
  count = count + 1
  print(count)
  factor_name <- names(pop_statistics_split)[count]
  print(factor_name)
  factor_name <-sub("pol.*", "", factor_name)
  pretty_name <- unique(factor$pretty_name)

  
  lambda <- factor$GroupValue[which.max( (1/sqrt(factor$b) )) ]
  print(lambda)
  new_factor <- ( (factor$GroupValue) ^ lambda - 1) / lambda
  
  Sd_plot <- ggplot(factor, aes(factor$GroupValue, log10(factor$S_d*-1))) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "firebrick4", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, log10(factor$S_d*-1), factor$GroupValue), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression('Absolute '*italic(S['d'])*', log transformed')) 
  
  b_plot <- ggplot(factor, aes(factor$GroupValue, factor$b)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "firebrick4", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, factor$b, factor$GroupValue), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(italic(b)))

  coeff_var_plot <- ggplot(factor, aes(factor$GroupValue, 1/sqrt(factor$b))) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "firebrick4", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, 1/sqrt(factor$b), factor$GroupValue), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = 'Coefficient of variation')  
  
  variance_plot <- ggplot(factor, aes(factor$GroupValue, factor$b*factor$S_d^2)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "firebrick4", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, factor$b*factor$S_d^2, factor$GroupValue), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = 'Variance')      
  
  alpha_plot <- ggplot(factor, aes(factor$GroupValue, factor$alpha)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "firebrick4", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, factor$alpha, factor$GroupValue), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = pretty_name, y = expression(alpha))  

  
  log_Sd_plot <- ggplot(factor, aes(log10(factor$GroupValue), log10(factor$S_d*-1))) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "firebrick4", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, log10(factor$S_d*-1), log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = paste(pretty_name, ' log transformed', sep = ','), y = expression('Absolute '*italic(S['d'])*', log transformed')) 
  
  log_coeff_var_plot <- ggplot(factor, aes(log10(factor$GroupValue), 1/sqrt(factor$b))) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "firebrick4", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, 1/sqrt(factor$b), log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = paste(pretty_name, ' log transformed', sep = ','), y = 'Coefficient of variation')    
  
  log_alpha_plot <- ggplot(factor, aes(log10(factor$GroupValue), factor$alpha)) +
    geom_smooth(method = "lm", col = "black") +
    geom_point(colour = "firebrick4", size = 2)+
    annotate("text", x = I(0), y = I(0), label = lm_vals(factor, factor$alpha, log10(factor$GroupValue)), vjust = -1, hjust = -0.1) +
    theme_classic() + labs(x = paste(pretty_name, ' log transformed', sep = ','), y = expression(alpha))   
   
  print(grid.arrange(Sd_plot, b_plot, coeff_var_plot, variance_plot, alpha_plot, ncol = 5, nrow = 1))
  print(grid.arrange(Sd_plot, coeff_var_plot, alpha_plot, Sd_plot, coeff_var_plot, alpha_plot, ncol = 3, nrow = 2))
  
  print(grid.arrange(Sd_plot, coeff_var_plot, alpha_plot, log_Sd_plot, log_coeff_var_plot, log_alpha_plot, ncol = 3, nrow = 2))
  
  
  