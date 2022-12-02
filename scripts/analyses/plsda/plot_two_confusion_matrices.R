#plot two confusion matrices together
library(RColorBrewer)

# Classic palette BuPu, with 4 colors
coul <- brewer.pal(11, "Spectral") 

# Add more colors to this palette :
coul <- colorRampPalette(coul)(2)

coul <- c("#FFFFFF", "#FFFFFF")

pie(rep(1, length(coul)), col = coul , main="") 

#read data
tab_mean_all <- read.csv("./analysis/plsda/PLSDA_all_confumean_19comps.csv", row.names = 1)
tabs_perc_all <- read.csv("./analysis/plsda/PLSDA_all_confuperc_19comps.csv", row.names = 1)

tab_mean_ang <- read.csv("./analysis/plsda/PLSDA_all_confumean_ang_only_20comps.csv", row.names = 1)
tabs_perc_ang <- read.csv("./analysis/plsda/PLSDA_all_confuperc_ang_only_20comps.csv", row.names = 1)

tab_perc_all_matrix <- as.matrix(tabs_perc_all)
tab_perc_ang_matrix <- as.matrix(tabs_perc_ang)

#make plot
jpeg("./output/plsda_confusion_matrix_AM_EM_all_vs_ang.jpg", height = 5, width = 8, res = 600, units = "in")
#pdf("./output/plsda_confusion_matrix_AM_EM_all_vs_ang.pdf", height = 5, width = 8)
par(oma = c(0,1,1,1))
#plot it
layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))#, widths=c(2, 1))
corrplot::corrplot(tab_perc_all_matrix, is.corr = T, col = c(coul, coul), tl.col = 1, cl.pos = "n", #
                   tl.offset =0.7, tl.cex = 0.9, tl.srt = 0, outline = TRUE,
                   addCoef.col ='black', number.font = 1, number.cex = 0.9, addCoefasPercent = T, na.label = " ", mar = c(0.3,2,0.3,1)) 

mtext("Prediction",2, line=3.1, cex=1.2)
mtext("Reference",3, line = 1.2, cex=1.2)
mtext("Classification accuracy (%)", line = -17, cex=0.9)
mtext("a)", 3, line = 3, cex=1.5, adj= 0)
#mtext("A", 3, line = 3, cex=2, adj = 0)

corrplot::corrplot(tab_perc_ang_matrix, is.corr = T,  col = c(coul, coul), tl.col = 1, cl.pos = "n", #
                   tl.offset =0.7, tl.cex = 0.9, tl.srt = 0,outline = TRUE,
                   addCoef.col ='black', number.font = 1, number.cex = 0.9, addCoefasPercent = T, na.label = " ", mar = c(0.3,2,0.3,1)) 

mtext("Prediction",2, line=3.1, cex=1.2)
mtext("Reference",3, line = 1.2, cex=1.2)
mtext("Classification accuracy (%)", line = -17, cex=0.9)
mtext("b)", 3, line = 3, cex=1.5,  adj = 0)
#mtext("B", 3, line = 3, cex=2,)

dev.off()



#################reading data differently
#train_fit <- readRDS("./analysis/plsda/PLSDA_fast_train_T6_ang_only_20comps.rds")
probis <- readRDS("./analysis/plsda/PLSDA_fast_T6_probis_ang_only_20comps.rds")
confus <- readRDS("./analysis/plsda/PLSDA_fast_T6_confus_ang_only_20comps.rds")


tabs <- list()
for(i in 1:length(confus)){
  tabs[[i]] <- confus[[i]]$table
}

tabsi <- Reduce('+', tabs)
tab_mean <- as.data.frame.matrix(tabsi/length(confus))

#write.csv(tab_mean,"./analysis/plsda/PLSDA_all_confumean_ang_only_19comps.csv")

#with percentage values
sums <- colSums(tab_mean)
tabs_perc <- matrix(NA, length(sums),length(sums))
for (i in 1:length(sums)){
  tabs_perc[,i] <- tab_mean[,i]/sums[i]
}

colnames(tabs_perc) <- colnames(confus[[1]]$table)
rownames(tabs_perc) <- rownames(confus[[1]]$table)

#write.csv(tabs_perc,"./analysis/plsda/PLSDA_all_confuperc_ang_only_20comps.csv")

str(tabs_perc_ang)
str(tabs_perc)

