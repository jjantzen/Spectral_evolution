#figuring out Raunkier life form

library(ape)
library(caper)
library(phytools)
library(dplyr)
library(plyr)
source("./external_code/delta_statistic-master/code.R")

#read data - still missing myc data
complete_data <- readRDS("./data/for_analysis/final_data.rds")

colnames(complete_data) <- c("Species", "Herb", "Shrub", "Tree", "Woody", "Subclass", "Superorder", "Order", "Family", "Shade", "Drought", "Coarse_soil", "Fine_soil", "Medium_soil", "min_pH", "max_pH", "leaf_persistence", "Raunk_lf", "Habitat", "Raunk_broad")

#get reduced dataset for only relevant traist
complete_data <- complete_data[,c(1:5,8:17,20)]
colnames(complete_data)

# #get only complete cases
# complete_data <- complete_data[complete.cases(complete_data),]

#load myc data

myc_data <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

myc_data_df <- as.data.frame(myc_data$myc)
myc_data_df <- cbind(myc_data$species, myc_data_df, myc_data$gf, myc_data$lp)

colnames(myc_data_df)[1] <- "Species"

#trim raunk to myc 
myc_data_combo <- merge(myc_data_df, complete_data[,c(1,16)], by = "Species", all = TRUE)

myc_data_combo

unique(myc_data_combo$Raunk_broad)

myc_data_combo$Raunk_broad[which(myc_data_combo$Raunk_broad == "h + c")] <- "c + h"
myc_data_combo$Raunk_broad[which(myc_data_combo$Raunk_broad == "h + p")] <- "p + h"
myc_data_combo$Raunk_broad[which(myc_data_combo$Raunk_broad == "ch + p")] <- "p + ch"
myc_data_combo$Raunk_broad[which(myc_data_combo$Raunk_broad == "c + p")] <- "c"
myc_data_combo$Raunk_broad[which(myc_data_combo$Raunk_broad == "p + ch")] <- "p"
myc_data_combo$Raunk_broad[which(myc_data_combo$Raunk_broad == "p + h")] <- "h"
myc_data_combo$Raunk_broad[which(myc_data_combo$Raunk_broad == "p + h")] <- "h"
myc_data_combo$Raunk_broad[which(myc_data_combo$Species == "Apocynum androsaemifolium")] <- "c"
myc_data_combo$Raunk_broad[which(myc_data_combo$Species == "Bromus inermis")] <- "h"
myc_data_combo$Raunk_broad[which(myc_data_combo$Species == "Calamagrostis canadensis")] <- "c"


unique(myc_data_combo$Raunk_broad)

#count number of cases
myc_data_combo %>% dplyr::count(Raunk_broad)

myc_data_combo$Species[which(myc_data_combo$Raunk_broad == "c + h")]
myc_data_combo$Species[which(myc_data_combo$Raunk_broad == "c + p")] #nope, only c
myc_data_combo$Species[which(myc_data_combo$Raunk_broad == "p + ch")] #only p
myc_data_combo$Species[which(myc_data_combo$Raunk_broad == "p + h")] #only h
myc_data_combo$Species[which(myc_data_combo$Raunk_broad == "c + h")] #only c,h,c


#calculate phylogenetic signal for raunk
new_trees <- readRDS("./data/for_analysis/final_trees_matched_spectra.rds")


df_output_rk <- data.frame(matrix(nrow=100, ncol = 3))
colnames(df_output_rk)[1] <- "iteration"
class(df_output_rk[,1]) <- "numeric"
colnames(df_output_rk)[2] <- "deltaA"
class(df_output_rk[,2]) <- "numeric"
colnames(df_output_rk)[3] <- "Pvalue"
class(df_output_rk[,3]) <- "numeric"

#read paper to figure out why those parameters below

for (j in 1:length(new_trees)){
  deltaA_gf <- delta(as.character(myc_data_combo$Raunk_broad),new_trees[[j]],0.1,0.0589,10000,10,100)
  
  random_delta <- rep(NA,100)
  for (i in 1:100){
    rtrait <- sample(myc_data_combo$Raunk_broad)
    random_delta[i] <- delta(rtrait,new_trees[[j]],0.1,0.0589,10000,10,100)
  }
  p_value <- sum(random_delta>deltaA_gf)/length(random_delta)
  
  df_output_rk$iteration[j] <- j
  df_output_rk[j,2] <- deltaA_gf
  df_output_rk[j,3] <- p_value
  print(j)
}



saveRDS(df_output_rk, "./analysis/trait_phylosig/raunk_gf_phylosig_100sp.rds")

sd(df_output_rk$Pvalue, na.rm = TRUE)
mean(df_output_rk$Pvalue, na.rm = TRUE)

nrow(df_output_rk[which((df_output_rk$Pvalue < 0.05) == TRUE),])


#plot raunk on tree
consensus_tree <- readRDS("./data/for_analysis/consensus_tree_full.rds")
colours_touse <- c("#008837", "#a6dba0", "#7b3294", "#F0E442")

Raunk <-setNames(as.factor(myc_data_combo$Raunk_broad),myc_data_combo$Species)
Raunk <- revalue(Raunk, c("p"="Phanerophyte", "h"="Hemicryptophyte", "c" = "Cryptophyte", "t" = "Therophyte"))
#Raunk <- relevel(Raunk, "Shade tolerant")#, "Intermediate shade tolerance", "Shade intolerance"))

jpeg("./output/trait_plots/Raunk_consensus.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse,levels(Raunk)) #palette()[c(4,2,3)]
plotTree(consensus_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(consensus_tree)],
       obj$yy[1:Ntip(consensus_tree)],pch=21,cex=1.2,
       bg=cols[Raunk[consensus_tree$tip.label]])
legend("bottomleft", levels(Raunk), pch=21, pt.bg=colours_touse, pt.cex=2) #palette()[c(4,2,3)]
dev.off()




###################

colnames(myc_data_combo) <- c("Species", "Myc", "GF", "LP", "Raunk")

complete_cases_data <- myc_data_combo[complete.cases(myc_data_combo),]


#contingency tables
table(myc_data_combo$Myc, myc_data_combo$GF) #1 shrub EM
table(myc_data_combo$Myc, myc_data_combo$LP) #
table(myc_data_combo$LP, myc_data_combo$GF) #1 shrub evergreen, 4 herb evergreen

table(myc_data_combo$Myc, myc_data_combo$Raunk) #4 AM t, 0 EM c or t

chisq.test(complete_cases_data$Myc, complete_cases_data$GF) #p= 1e-06
chisq.test(complete_cases_data$Myc, complete_cases_data$LP) #p= 0.0158
chisq.test(complete_cases_data$GF, complete_cases_data$LP) #p= 0.5494

chisq.test(myc_data_combo$GF, myc_data_combo$LP) #p= 0.5494
