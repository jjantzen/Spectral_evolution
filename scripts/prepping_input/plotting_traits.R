#plotting trait data to look for distribution and outliers

library(dplyr)
library(phytools)


#read data
combo_data <- readRDS("./data/tidy/new_combo_data_matched_spectra.rds")

new_trees <- readRDS("./data/tidy/new_trees_matched_spectra.rds")


#choose traits to map on tree
#MycType
#NitFix_Capacity
#DroughtTolerance_Qual
#ShadeTolerance_Quant

colnames(combo_data)

#make all traits factors
combo_data$ShadeTolerance_Qual <- as.factor(combo_data$ShadeTolerance_Qual)
combo_data$DroughtTolerance_Qual <- as.factor(combo_data$DroughtTolerance_Qual)
combo_data$NitFix_Capacity <- as.factor(combo_data$NitFix_Capacity)
combo_data$MycType <- as.factor(combo_data$MycType)
combo_data$Myc_AM <- as.factor(combo_data$Myc_AM)
combo_data$Myc_EM <- as.factor(combo_data$Myc_EM)
combo_data$Myc_ErM <- as.factor(combo_data$Myc_ErM)
combo_data$Myc_NM <- as.factor(combo_data$Myc_NM)
combo_data$Woody <- as.factor(combo_data$Woody)
combo_data$Superorder <- as.factor(combo_data$Superorder)
combo_data$Order <- as.factor(combo_data$Order)
combo_data$`Adapted to Coarse Textured Soils` <- as.factor(combo_data$`Adapted to Coarse Textured Soils`)
combo_data$`Adapted to Fine Textured Soils` <- as.factor(combo_data$`Adapted to Fine Textured Soils`)
combo_data$`Adapted to Medium Textured Soils` <- as.factor(combo_data$`Adapted to Medium Textured Soils`)
combo_data$`pH, Minimum` <- as.numeric(combo_data$`pH, Minimum`)
combo_data$`pH, Maximum` <- as.numeric(combo_data$`pH, Maximum`)

#rename some columns
colnames(combo_data)[38] <- "coarse_soil"
colnames(combo_data)[39] <- "fine_soil"
colnames(combo_data)[40] <- "medium_soil"
colnames(combo_data)[41] <- "min_pH"
colnames(combo_data)[42] <- "max_pH"


#make data into named lists
Myc_AM <-setNames(combo_data$Myc_AM,combo_data$Species)
Myc_ErM <-setNames(combo_data$Myc_ErM,combo_data$Species)
Myc_EM <-setNames(combo_data$Myc_EM,combo_data$Species)
Myc_NM <-setNames(combo_data$Myc_NM,combo_data$Species)
MycType <-setNames(combo_data$MycType,combo_data$Species)
NitFix_Capacity <-setNames(combo_data$NitFix_Capacity,combo_data$Species)
ShadeTolerance_Qual <-setNames(combo_data$ShadeTolerance_Qual,combo_data$Species)
DroughtTolerance_Qual <-setNames(combo_data$DroughtTolerance_Qual,combo_data$Species)
Woody <-setNames(combo_data$Woody,combo_data$Species)
Superorder <-setNames(combo_data$Superorder,combo_data$Species)
Order <-setNames(combo_data$Order,combo_data$Species)
coarse_soil <-setNames(combo_data$coarse_soil,combo_data$Species)
fine_soil <-setNames(combo_data$fine_soil,combo_data$Species)
medium_soil <-setNames(combo_data$medium_soil,combo_data$Species)
min_pH <-setNames(combo_data$min_pH,combo_data$Species)
max_pH <-setNames(combo_data$max_pH,combo_data$Species)

#convert to matrices
Myc_AM_matrix <- to.matrix(Myc_AM,levels(Myc_AM))
Myc_ErM_matrix <- to.matrix(Myc_ErM,levels(Myc_ErM))
Myc_EM_matrix <- to.matrix(Myc_EM,levels(Myc_EM))
Myc_NM_matrix <- to.matrix(Myc_NM,levels(Myc_NM))
MycType_matrix <- to.matrix(MycType,levels(MycType))
NitFix_Capacity_matrix <- to.matrix(NitFix_Capacity,levels(NitFix_Capacity))
ShadeTolerance_Qual_matrix <- to.matrix(ShadeTolerance_Qual,levels(ShadeTolerance_Qual))
DroughtTolerance_Qual_matrix <- to.matrix(DroughtTolerance_Qual,levels(DroughtTolerance_Qual))
Woody_matrix <- to.matrix(Woody,levels(Woody))
Superorder_matrix <- to.matrix(Superorder,levels(Superorder))
Order_matrix <- to.matrix(Order,levels(Order))
coarse_soil_matrix <- to.matrix(coarse_soil,levels(Order))
fine_soil_matrix <- to.matrix(fine_soil,levels(Order))
medium_soil_matrix <- to.matrix(medium_soil,levels(Order))

#make continuous data into matrices
row.names(combo_data) <- combo_data$Species
min_pH_matrix <- matrix(combo_data$min_pH, dimnames = list(row.names(combo_data)))
max_pH_matrix <- matrix(combo_data$max_pH, dimnames = list(row.names(combo_data)))

#replace NA with 0 for visualization
min_pH_matrix[which(is.na(min_pH_matrix))] <- 0
max_pH_matrix[which(is.na(max_pH_matrix))] <- 0

#combine Myc into one trait
Myc_data <- combo_data[,c(1,16:19)]
Myc_data$poly <- paste0(Myc_data$Myc_AM, "+", Myc_data$Myc_EM, "+", Myc_data$Myc_NM, "+", Myc_data$Myc_ErM)


str(Myc_data)




#plot
jpeg("./output/trait_plots/Myc_AM.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2)],levels(Myc_AM))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[Myc_AM[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(Myc_AM), pch=21, pt.bg=palette()[c(4,2)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Myc_ErM.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2)],levels(Myc_ErM))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[Myc_ErM[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(Myc_ErM), pch=21, pt.bg=palette()[c(4,2)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Myc_EM.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2)],levels(Myc_EM))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[Myc_EM[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(Myc_EM), pch=21, pt.bg=palette()[c(4,2)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Myc_NM.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2)],levels(Myc_NM))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[Myc_NM[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(Myc_NM), pch=21, pt.bg=palette()[c(4,2)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/NitFix_Capacity.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2)],levels(NitFix_Capacity))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[NitFix_Capacity[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(NitFix_Capacity), pch=21, pt.bg=palette()[c(4,2)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/ShadeTolerance_Qual.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2, 3)],levels(ShadeTolerance_Qual))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[ShadeTolerance_Qual[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(ShadeTolerance_Qual), pch=21, pt.bg=palette()[c(4,2,3)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/DroughtTolerance_Qual.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2,3,1)],levels(DroughtTolerance_Qual))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[DroughtTolerance_Qual[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(DroughtTolerance_Qual), pch=21, pt.bg=palette()[c(4,2,3,1)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Woody.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2,3,1)],levels(Woody))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[Woody[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(Woody), pch=21, pt.bg=palette()[c(4,2,3,1)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Superorder.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2,3,1,6)],levels(Superorder))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[Superorder[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(Superorder), pch=21, pt.bg=palette()[c(4,2,3,1,6)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Order.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(c("red", "blue", "green", "black", "yellow", "orange", "purple","red", "blue", "green", "black", "yellow", "orange", "purple","red", "orange", "green", "black", "yellow", "blue", "purple","black"),levels(Order))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[Order[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(Order), pch=21, pt.bg=c("red", "blue", "green", "black", "yellow", "orange", "purple","red", "blue", "green", "black", "yellow", "orange", "purple","red", "orange", "green", "black", "yellow", "blue", "purple","black"), pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Coarse_soil.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2)],levels(coarse_soil))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[coarse_soil[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(coarse_soil), pch=21, pt.bg=palette()[c(4,2)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Fine_soil.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2)],levels(fine_soil))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[fine_soil[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(fine_soil), pch=21, pt.bg=palette()[c(4,2)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Medium_soil.jpg", height = 10, width = 6, res = 500, units = "in")
cols<-setNames(palette()[c(4,2)],levels(medium_soil))
plotTree(new_trees[[1]],ftype="i",offset=0.6,fsize=0.6)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[1:Ntip(new_trees[[1]])],
       obj$yy[1:Ntip(new_trees[[1]])],pch=21,cex=1.2,
       bg=cols[medium_soil[new_trees[[1]]$tip.label]])
legend("bottomleft", levels(medium_soil), pch=21, pt.bg=palette()[c(4,2)], pt.cex=2)
dev.off()

jpeg("./output/trait_plots/Min_pH.jpg", height = 10, width = 6, res = 500, units = "in")
dotTree(new_trees[[1]],min_pH_matrix,length=10,ftype="i", fsize = 0.5)
dev.off()

jpeg("./output/trait_plots/Max_pH.jpg", height = 10, width = 6, res = 500, units = "in")
dotTree(new_trees[[1]],max_pH_matrix,length=10,ftype="i", fsize = 0.5)
dev.off()


#need finite ylim values error for both
phenogram(new_trees[[1]],min_pH_matrix,spread.labels=TRUE,spread.cost=c(1,0))
plotTree.barplot(new_trees[[1]],min_pH_matrix)



unique(combo_data$min_pH)

#obj<-contMap(new_trees[[1]],min_pH_matrix,plot=FALSE)
#obj<-setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.4,1),outline=FALSE,lwd=c(3,7),leg.txt="log(SVL)")


#polymorphic doesn't work because of NAs
poly_myc <- Myc_data$poly
poly_myc <-setNames(poly_myc,Myc_data$Species)

#get rid of NAs
poly_myc

#for polymorphic traits
plotTree(new_trees[[1]],ftype="off",lwd=1,type="fan")
poly_myc <- strsplit(setNames(as.character(poly_myc),names(poly_myc)),"+",fixed=TRUE)
pies<-matrix(0,Ntip(new_trees[[1]]),4,dimnames=list(new_trees[[1]]$tip.label,c("Myc_AM","Myc_EM","Myc_NM","Myc_ErM")))
for(i in 1:Ntip(new_trees[[1]])){
  pies[new_trees[[1]]$tip.label[i],poly_myc[[new_trees[[1]]$tip.label[i]]]]<-rep(1/length(poly_myc[[new_trees[[1]]$tip.label[i]]]), length(poly_myc[[new_trees[[1]]$tip.label[i]]]))
}

tiplabels(pie=pies,piecol=c("black","yellow","red","blue"),cex=0.35)
legend(x="topleft",legend=c("Myc_AM","Myc_EM","Myc_NM","Myc_ErM"),pt.cex=2,pch=21,
       pt.bg=c("black","yellow","red","blue"))
