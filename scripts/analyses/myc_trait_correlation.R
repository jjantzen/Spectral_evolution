#trait only models - for most general exploratory model predictors
library(ggplot2)
library(tidyr)

library(ape)
library(caper)
library(phytools)
source("./external_code/delta_statistic-master/code.R")

#read data - still missing myc data
complete_data <- readRDS("./data/for_analysis/final_data.rds")
myc_data <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

#get dataframe for 92sp
myc_data_df <- as.data.frame(myc_data[c(2:5)])
rownames(myc_data_df) <- myc_data_df$species

#add woodiness
myc_data_df$woody <- complete_data$Woody[which(complete_data$Species %in% myc_data_df$species)]

#chi squared tests
chisq.test(myc_data_df$woody, myc_data_df$lp) #0.9681
chisq.test(myc_data_df$woody, myc_data_df$gf) #<2.2e-16
chisq.test(myc_data_df$woody, myc_data_df$myc) #7.8e-05

chisq.test(myc_data_df$lp, myc_data_df$gf) #0.5494
chisq.test(myc_data_df$lp, myc_data_df$myc) #0.0158

chisq.test(myc_data_df$gf, myc_data_df$myc) #1.0e-06

#contingency tables

table(myc_data_df$woody, myc_data_df$lp) #4 evergreen nonwoody
table(myc_data_df$woody, myc_data_df$gf) #0 woody herbs, herb shrub/tree
table(myc_data_df$woody, myc_data_df$myc) #0 EM herbs
table(myc_data_df$lp, myc_data_df$gf) #1 evergreen shrub, 4 herb evergreen
table(myc_data_df$lp, myc_data_df$myc) #
table(myc_data_df$gf, myc_data_df$myc) #1 EM shrub, 0 EM herb


#phylosig redo 

myc_trees <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

myc_data_df$woody

#prep output
df_output_woody <- data.frame(matrix(nrow=100, ncol = 4))
colnames(df_output_woody)[1] <- "iteration"
class(df_output_woody[,1]) <- "numeric"
colnames(df_output_woody)[2] <- "Estimated D"
class(df_output_woody[,2]) <- "numeric"
colnames(df_output_woody)[3] <- "Prob of E(D) from no (random) phylo structure"
class(df_output_woody[,3]) <- "numeric"
colnames(df_output_woody)[4] <- "Prob of E(D) from BM phylo structure"
class(df_output_woody[,4]) <- "numeric"

for (i in 1:length(myc_trees)){
  out <- phylo.d(data = myc_data_df, phy = myc_trees[[i]], names.col = species, binvar = woody, permut = 999) #D = -0.5893118, P of E(D) BM = 0.978979, P of E(D) non BM = 0
  df_output_woody$iteration[i] <- i
  df_output_woody[i,2] <- out$DEstimate
  df_output_woody[i,3] <- out$Pval1
  df_output_woody[i,4] <- out$Pval0
}

df_output_woody
saveRDS(df_output_woody, "./analysis/trait_phylosig/woody_phylosig_92sp.rds")

mean(df_output_woody$`Prob of E(D) from BM phylo structure`)

nrow(df_output_woody[which((df_output_woody$`Prob of E(D) from BM phylo structure` > 0.95) == TRUE),])
