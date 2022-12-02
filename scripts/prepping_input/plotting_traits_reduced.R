#plot traits for consensus trees

library(dplyr)
library(phytools)
library(RColorBrewer)

#read data
trait_data <- readRDS("./data/for_analysis/trait_data_100_taxa.rds")

consensus_tree <- read.nexus("./data/for_analysis/consensus_tree.nex")

consensus_tree$tip.label <- gsub("_", " ", consensus_tree$tip.label)

trait_data$Species[-which(trait_data$Species %in% consensus_tree$tip.label)]
consensus_tree$tip.label[-which(consensus_tree$tip.label %in% trait_data$Species)]

consensus_tree$tip.label[which(consensus_tree$tip.label == "'Alnus incana subsp. rugosa'")] <- "Alnus incana"
consensus_tree$tip.label[which(consensus_tree$tip.label == "'Acer saccharum subsp. nigrum'")] <- "Acer nigrum"

#edit tip names in consensus tree to match dataset

plot(consensus_tree)

myc_data <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

#tips to drop 
exclude <- consensus_tree$tip.label[-which(consensus_tree$tip.label %in%myc_data$species)]

myc_tree <- drop.tip(consensus_tree,tip=exclude)

#saveRDS(myc_tree, "./data/for_analysis/myc_consensus_tree.rds")

# #drop tips
# tree <- drop.tip(tree, exclude)
# plot(tree)

#make data into named lists
Shade <-setNames(trait_data$Shade,trait_data$Species)
Drought <-setNames(trait_data$Drought_bin,trait_data$Species)
Woody <-setNames(trait_data$Woody,trait_data$Species)
Coarse_soil <-setNames(trait_data$Coarse_soil,trait_data$Species)
Fine_soil <-setNames(trait_data$Fine_soil,trait_data$Species)
Leaf_persistence <-setNames(trait_data$leaf_persistence,trait_data$Species)


#change level names and orders
library(plyr)
Shade <- revalue(Shade, c("Tolerant"="Shade tolerant", "Intermediate"="Intermediate shade tolerance", "Intolerant" = "Shade intolerance"))
Shade <- relevel(Shade, "Shade tolerant")#, "Intermediate shade tolerance", "Shade intolerance"))

Drought <- revalue(Drought, c("No"="Drought intolerant", "Yes"="Drought tolerant"))
Drought <- relevel(Drought, "Drought tolerant")

Woody <- revalue(Woody, c("0"="Non-woody", "1"="Woody"))
Woody <- relevel(Woody, "Non-woody")

Fine_soil <- revalue(Fine_soil, c("No"="Not adapted to fine soil", "Yes"="Adapted to fine soil"))
Fine_soil <- relevel(Fine_soil, "Adapted to fine soil")

Coarse_soil <- revalue(Coarse_soil, c("No"="Not adapted to coarse soil", "Yes"="Adapted to coarse soil"))
Coarse_soil <- relevel(Coarse_soil, "Adapted to coarse soil")

Leaf_persistence <- revalue(Leaf_persistence, c("d"="Deciduous", "e"="Evergreen"))
Leaf_persistence <- relevel(Leaf_persistence, "Evergreen")

#colours

#three colours
colours_touse <- c("#008837", "#a6dba0", "#7b3294")

#colours_touse <- c("#008837", "#CC79A7", "#009E73")


##         black        orange       skyblue   bluishgreen        yellow 
##     "#000000"     "#E69F00"     "#56B4E9"     "#009E73"     "#F0E442" 
##          blue    vermillion reddishpurple          gray 
##     "#0072B2"     "#D55E00"     "#CC79A7"     "#999999" 

#plot
jpeg("./output/trait_plots/Shade_consensus.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse,levels(Shade)) #palette()[c(4,2,3)]
plotTree(consensus_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(consensus_tree)],
       obj$yy[1:Ntip(consensus_tree)],pch=21,cex=1.2,
       bg=cols[Shade[consensus_tree$tip.label]])
legend("bottomleft", levels(Shade), pch=21, pt.bg=colours_touse, pt.cex=2) #palette()[c(4,2,3)]
dev.off()

jpeg("./output/trait_plots/Drought_consensus.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse[c(1,3)],levels(Drought))#palette()[c(4,2,3)]
plotTree(consensus_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(consensus_tree)],
       obj$yy[1:Ntip(consensus_tree)],pch=21,cex=1.2,
       bg=cols[Drought[consensus_tree$tip.label]])
legend("bottomleft", levels(Drought), pch=21, pt.bg=colours_touse[c(1,3)], pt.cex=2)#palette()[c(4,2,3)]
dev.off()

jpeg("./output/trait_plots/Woody_consensus.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse[c(1,3)],levels(Woody)) #palette()[c(4,2,3)]
plotTree(consensus_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(consensus_tree)],
       obj$yy[1:Ntip(consensus_tree)],pch=21,cex=1.2,
       bg=cols[Woody[consensus_tree$tip.label]])
legend("bottomleft", levels(Woody), pch=21, pt.bg=colours_touse[c(1,3)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Coarse_soil_consensus.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse[c(1,3)],levels(Coarse_soil))
plotTree(consensus_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(consensus_tree)],
       obj$yy[1:Ntip(consensus_tree)],pch=21,cex=1.2,
       bg=cols[Coarse_soil[consensus_tree$tip.label]])
legend("bottomleft", levels(Coarse_soil), pch=21, pt.bg=colours_touse[c(1,3)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Fine_soil_consensus.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse[c(1,3)],levels(Fine_soil))
plotTree(consensus_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(consensus_tree)],
       obj$yy[1:Ntip(consensus_tree)],pch=21,cex=1.2,
       bg=cols[Fine_soil[consensus_tree$tip.label]])
legend("bottomleft", levels(Fine_soil), pch=21, pt.bg=colours_touse[c(1,3)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Leaf_persistence_consensus.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse[c(1,3)],levels(Leaf_persistence))
plotTree(consensus_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(consensus_tree)],
       obj$yy[1:Ntip(consensus_tree)],pch=21,cex=1.2,
       bg=cols[Leaf_persistence[consensus_tree$tip.label]])
legend("bottomleft", levels(Leaf_persistence), pch=21, pt.bg=colours_touse[c(1,3)], pt.cex=2)
dev.off()

Myc <- setNames(as.factor(myc_data$myc), myc_data$species)
Myc
# Drought <- revalue(Drought, c("No"="Drought intolerant", "Yes"="Drought tolerant"))
# Drought <- relevel(Drought, "Drought tolerant")
#make it orange and blue
colours_touse <- c("orange", "blue", "#7b3294")
colours_touse_3 <- c("#008837", "#a6dba0", "#7b3294")

shapes_to_use <- c(21, 22, 23)

jpeg("./output/trait_plots/Myc_consensus_orange.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse[c(1,2)],levels(Myc)) #or 1,3
shaps <- setNames(shapes_to_use[c(1,2,3)], levels(Myc))
plotTree(myc_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(myc_tree)],
       obj$yy[1:Ntip(myc_tree)],pch=shaps[Myc[myc_tree$tip.label]],cex=1.2,
       bg=cols[Myc[myc_tree$tip.label]])
legend("bottomleft", levels(Myc), pch=shaps, pt.bg=colours_touse[c(1,2)], pt.cex=2) #or 1,3
dev.off()

GF <- setNames(as.factor(myc_data$gf), myc_data$species)

jpeg("./output/trait_plots/GF_myc_consensus.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse_3[c(3,2,1)],levels(GF))
shaps <- setNames(shapes_to_use[c(1,2,3)], levels(GF))
plotTree(myc_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(myc_tree)],
       obj$yy[1:Ntip(myc_tree)],pch=shaps[GF[myc_tree$tip.label]],cex=1.2,
       bg=cols[GF[myc_tree$tip.label]])
legend("bottomleft", levels(GF), pch=shaps, pt.bg=colours_touse_3[c(3,2,1)], pt.cex=2)
dev.off()

LP <- setNames(as.factor(myc_data$lp), myc_data$species)
LP <- revalue(LP, c("deciduous"="Deciduous", "evergreen"="Evergreen"))

jpeg("./output/trait_plots/LP_myc_consensus.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse_3[c(3,1)],levels(LP))
shaps <- setNames(shapes_to_use[c(1,2)], levels(LP))
plotTree(myc_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(myc_tree)],
       obj$yy[1:Ntip(myc_tree)],pch=shaps[LP[myc_tree$tip.label]],cex=1.2,
       bg=cols[LP[myc_tree$tip.label]])
legend("bottomleft", levels(LP), pch=shaps, pt.bg=colours_touse_3[c(3,1)], pt.cex=2)
dev.off()


Woody_myc <- Woody[which(names(Woody) %in% names(Myc))]
length(Woody_myc)


jpeg("./output/trait_plots/Woody_myc_consensus.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(colours_touse_3[c(1,3)],levels(Woody_myc))
shaps <- setNames(shapes_to_use[c(1,2)], levels(Woody_myc))
plotTree(myc_tree,ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(myc_tree)],
       obj$yy[1:Ntip(myc_tree)],pch=shaps[Woody_myc[myc_tree$tip.label]],cex=1.2,
       bg=cols[Woody_myc[myc_tree$tip.label]])
legend("bottomleft", levels(Woody_myc), pch=shaps, pt.bg=colours_touse_3[c(1,3)], pt.cex=2)
dev.off()

##########

Myc <- revalue(Myc, c("AM"="Arbuscular mycorrhizal", "EM"="Ectomycorrhizal"))
colours_touse_2 <- c("green4", "salmon4")


#plot into one file
jpeg("./output/trait_plots/4_predictors_distribution_plot.jpg", height = 10, width = 5, units = "in", res = 600)

#make pdf for publication
#pdf("./output/trait_plots/4_predictors_distribution_plot.pdf", height = 10, width = 5)


#par(mfrow=c(1,4))
#layout(matrix(c(1,1,1,1,1,1,1,2), nrow = 1, ncol = 8, byrow = TRUE))

#woodiness
#cols<-setNames(colours_touse_3[c(3,1)],levels(Woody_myc))
#shaps <- setNames(shapes_to_use[c(1,2)], levels(Woody_myc))
plotTree(myc_tree,ftype="i",offset=2,fsize=0.45, xlim = c(0,725))#, mar = c(0.1,0.1,0.1,3.5)
#obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
#points(obj$xx[1:Ntip(myc_tree)],
#       obj$yy[1:Ntip(myc_tree)],pch=shaps[Woody_myc[myc_tree$tip.label]],cex=1,
#       bg=cols[Woody_myc[myc_tree$tip.label]])
#legend(-5,25, levels(Woody_myc), pch=shaps, pt.bg=colours_touse_3[c(3,1)], pt.cex=2, title = "Woodiness\n", bty = "n", cex = 0.5)
#legend(-5,10, levels(Woody_myc), pch=shaps, pt.bg=colours_touse_3[c(3,1)], pt.cex=2, title = "Woodiness\n", bty = "n")

#GF
cols<-setNames(colours_touse_3[c(3,2,1)],levels(GF))
shaps <- setNames(shapes_to_use[c(1,2,3)], levels(GF))
#plotTree(myc_tree,ftype="i",offset=0.6,fsize=0.45)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(myc_tree)],
       obj$yy[1:Ntip(myc_tree)],pch=shaps[GF[myc_tree$tip.label]],cex=1,
       bg=cols[GF[myc_tree$tip.label]])
legend(-5,25, levels(GF), pch=shaps, pt.bg=colours_touse_3[c(3,2,1)], pt.cex=1.5, title = "Growth\nForm", bty = "n", cex = 0.75, y.intersp = 1, title.adj = 0)
#legend(75,10, levels(GF), pch=shaps, pt.bg=colours_touse_3[c(3,2,1)], pt.cex=2, title = "Growth\nForm", bty = "n")

#LP
cols<-setNames(colours_touse_2[c(2,1)],levels(LP))
shaps <- setNames(shapes_to_use[c(1,2)], levels(LP))
#plotTree(myc_tree,ftype="i",offset=0.6,fsize=0.45)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(myc_tree)]+20,
       obj$yy[1:Ntip(myc_tree)],pch=shaps[LP[myc_tree$tip.label]],cex=1,
       bg=cols[LP[myc_tree$tip.label]])
legend(-5,15, levels(LP), pch=shaps, pt.bg=colours_touse_2[c(2,1)], pt.cex=1.5, title = "Leaf\nPersistence", bty = "n", cex = 0.75, y.intersp = 1, title.adj = 0)
#legend(130,10, levels(LP), pch=shaps, pt.bg=colours_touse_3[c(3,1)], pt.cex=2, title = "Leaf\nPersistence", bty = "n")

#myc
cols<-setNames(colours_touse[c(1,2)],levels(Myc)) #or 1,3
shaps <- setNames(shapes_to_use[c(1,2,3)], levels(Myc))
#plotTree(myc_tree,ftype="i",offset=0.6,fsize=0.45)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(myc_tree)]+40,
       obj$yy[1:Ntip(myc_tree)],pch=shaps[Myc[myc_tree$tip.label]],cex=1,
       bg=cols[Myc[myc_tree$tip.label]])
legend(-5,7, levels(Myc), pch=shaps, pt.bg=colours_touse[c(1,2)], pt.cex=1.5, title = "Mycorrhizal\nType", bty = "n", cex = 0.75, y.intersp = 1, title.adj = 0, xjust = 0) #or 1,3
#legend(200,10, levels(Myc), pch=shaps, pt.bg=colours_touse[c(1,2)], pt.cex=2, title = "Mycorrhizal\nType", bty = "n") #or 1,3

legend(357,94, legend = c("GF", "LP", "MT"), bty = "n", horiz = TRUE, text.width=c(0,1.3,0.4), cex = 0.5, x.intersp = 0.7)#c(0,1.5,0.5)

cladelabels(tree = myc_tree, text = "Pinidae", node=findMRCA(myc_tree,c("Thuja occidentalis", "Pinus banksiana")), offset=9.2, wing.length=1, cex=1)
cladelabels(tree = myc_tree, text = "Magnoliidae", node=findMRCA(myc_tree,c("Camassia quamash", "Sorbus americana")), offset=6, wing.length=1, cex=1)

cladelabels(tree = myc_tree, text = "Lilianae", node=findMRCA(myc_tree,c("Camassia quamash", "Dactylis glomerata")), offset=6, wing.length=0.8, cex=0.8)
cladelabels(tree = myc_tree, text = "Rosanae", node=findMRCA(myc_tree,c("Tilia americana", "Sorbus americana")), offset=6.5, wing.length=0.8, cex=0.8)
cladelabels(tree = myc_tree, text = "Asteranae", node=findMRCA(myc_tree,c("Euthamia graminifolia", "Cornus sericea")), offset=5, wing.length=0.8, cex=0.8)
#cladelabels(text = "Asteranae", node=findMRCA(myc_tree,c("Euthamia graminifolia", "Cornus sericea")), offset=NULL, wing.length=0.1, cex=1)


dev.off()

