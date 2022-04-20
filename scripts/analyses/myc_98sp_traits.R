#4-state myc dataset trait phylosignal and trees
library(ape)
library(caper)
library(phytools)
library(dplyr)
library(plyr)
source("./external_code/delta_statistic-master/code.R")

#load data
myc_data <- readRDS("./data/for_analysis/myc_data_list_98sp_4cat_for_analysis.rds")

#read trees
consensus_tree <- readRDS("./data/for_analysis/consensus_tree_full.rds")

#tips to drop (for myc trees)
exclude <- consensus_tree$tip.label[-which(consensus_tree$tip.label %in% dimnames(myc_data$spectra)[[1]])]

myc_tree <- drop.tip(consensus_tree,tip=exclude)

#plot trees
colours_touse <- c("#008837", "#a6dba0", "#7b3294", "#F0E442")

#plot myc
Myc <-setNames(as.factor(myc_data$myc),myc_data$species)
#Myc <- revalue(Myc, c("AM"="Arbuscular Mycorrhizal", "EM"="Ectomycorrhizal"))

jpeg("./output/trait_plots/Myc_98sp_consensus.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse,levels(Myc)) #palette()[c(4,2,3)]
plotTree(myc_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(myc_tree)],
       obj$yy[1:Ntip(myc_tree)],pch=21,cex=1.2,
       bg=cols[Myc[myc_tree$tip.label]])
legend("bottomleft", levels(Myc), pch=21, pt.bg=colours_touse, pt.cex=2) #palette()[c(4,2,3)]
dev.off()
