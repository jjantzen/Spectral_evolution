#how to summarize across 100 replicates of PCAs?

#Can plot on consensus tree (with single mapping of states)
#Can I average PC1 scores????? Or compare scores per pc plot to see how much variation
library(phytools)
library(Rmisc)
library(tidyr)
library(ggplot2)
library(mvMORPH)

#saveRDS(model_list, "./analysis/pca_analysis/best_intercept_models_for_pcas_92sp.rds")
model_list <- readRDS("./analysis/pca_analysis/best_intercept_models_for_pcas_92sp.rds")


#get pca scores for all replicates into one dataframe
for (i in 1:100){#length(model_list)){
  #get pca
  print(i)
  pca_best_model <- mvgls.pca(model_list[[i]], plot = FALSE)  
  
  #structure data
  data_from_pcas <- as.data.frame(pca_best_model$scores[,c(1:200)])
  data_from_pcas$species <- rownames(pca_best_model$scores)
  data_from_pcas$iteration <- i
  
  #combine data
  if (i == 1){
    combo <- data_from_pcas
  } else {
    combo <- rbind(combo, data_from_pcas)
  } 
}

nrow(combo)

# #compare scores for pc 1, 2, 3 for trees 1:6
# plot(pca_best_model_1$scores[,1])
# points(pca_best_model_2$scores[,1], col = "blue")
# points(pca_best_model_3$scores[,1], col = "red")
# points(pca_best_model_4$scores[,1], col = "green")
# points(pca_best_model_5$scores[,1], col = "purple")


# #how to structure results for ggplot?
# data_from_pcas <- as.data.frame(pca_best_model_1$scores)
# data_from_pcas$species <- rownames(pca_best_model_1$scores)
# data_from_pcas$iteration <- 1
# 
# data_from_pcas_2 <- as.data.frame(pca_best_model_2$scores)
# data_from_pcas_2$species <- rownames(pca_best_model_2$scores)
# data_from_pcas_2$iteration <- 2
# 
# combo <- rbind(data_from_pcas, data_from_pcas_2)

#prep tree if needed
consensus_tree <- readRDS("./data/for_analysis/myc_consensus_tree.rds")

#read in data
data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")
tree_myc <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

tips_to_drop <- consensus_tree$tip.label[-which(consensus_tree$tip.label %in% data_spectra$species)]

consensus_tree_pruned <- drop.tip(consensus_tree, tip = tips_to_drop)

#reorder species by tree order 
combo$species_factor <- factor(combo$species,
                levels = consensus_tree_pruned$tip.label)



# #save combo output
# saveRDS(combo, "./analysis/pca_analysis/scores_for_all_reps_and_pcs.rds")
combo <- readRDS("./analysis/pca_analysis/scores_for_all_reps_and_pcs.rds")

top_half_sp <- levels(combo$species_factor)[c(1:46)]

combo$plotting_index[which(combo$species_factor %in% top_half_sp)] <- "A"
combo$plotting_index[-which(combo$species_factor %in% top_half_sp)] <- "B"


#plot scores across replicates and species to look at variation
jpeg("./output/PCAs/consensus_plots/scores_by_species_across_tree_reps_jitter.jpg", width = 10, height = 8, units = "in", res = 600)
ggplot(combo)+
  #geom_boxplot(aes(x = species_factor, y = V1))+
  geom_jitter(aes(x = species_factor, y = V1, colour = iteration), size = 0.1)+
  labs(x = "Species sorted by phylogeny", y = "PC scores")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 5))+
  facet_wrap(~ plotting_index, scales = "free_x", nrow = 2)
dev.off()

#get summary stats for scores for plotting on tree (mean and 95CI) as barplots
tail(colnames(combo))
length(colnames(combo))

mean_sd <- list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))

stats_for_all_iterations <- combo %>% 
  dplyr::select(c(1:6, 201:204)) %>% 
  dplyr::group_by(species_factor) %>% 
  dplyr::summarise(across(V1:V6, mean_sd))

for (i in 1:6){
  combo_data_sm <- combo[,c(i,203, 202)]
  data_wide <- as.data.frame(pivot_wider(combo_data_sm, names_from = "species_factor", values_from = colnames(combo_data_sm)[1]))
  
  CI_df <- data.frame(matrix(nrow=92, ncol = 5))
  colnames(CI_df)[1] <- "upper"
  class(CI_df[,1]) <- "numeric"
  colnames(CI_df)[2] <- "mean"
  class(CI_df[,2]) <- "numeric"
  colnames(CI_df)[3] <- "lower"
  class(CI_df[,3]) <- "numeric"
  colnames(CI_df)[4] <- "pc_axis"
  class(CI_df[,4]) <- "numeric"
  colnames(CI_df)[5] <- "species_factor"
  class(CI_df[,5]) <- "character"
  #rownames(CI_df) <- colnames(data_wide)[2:93]
  data_wide <- data_wide[,-1]
  for (j in 1:ncol(data_wide)){
    CI_results <- CI(data_wide[,j], ci = 0.95)
    CI_df$upper[j] <- CI_results[1]
    CI_df$mean[j] <- CI_results[2]
    CI_df$lower[j] <- CI_results[3]
    CI_df$pc_axis[j] <- i
    CI_df$species_factor[j] <- colnames(data_wide)[j]
  }
  
  if (i == 1){
    all_pcs_CIs <- CI_df  
  } else {
    all_pcs_CIs <- rbind(all_pcs_CIs, CI_df)
  }
  
}

all_pcs_CIs

pc1 <- all_pcs_CIs[which(all_pcs_CIs$pc_axis == 1),]

pc2 <- all_pcs_CIs[which(all_pcs_CIs$pc_axis == 2),]

pc3 <- all_pcs_CIs[which(all_pcs_CIs$pc_axis == 3),]


rownames(pc1) <- pc1$species_factor
rownames(pc2) <- pc2$species_factor
rownames(pc3) <- pc3$species_factor

pc1$middle_tick <- pc1$mean - pc1$lower
pc2$middle_tick <- pc2$mean - pc2$lower
pc3$middle_tick <- pc3$mean - pc3$lower

#ended up using this metric (max minus min) to get length of bar to add to min bar for plotting
pc1$last_tick <- pc1$upper - pc1$lower
pc2$last_tick <- pc2$upper - pc2$lower
pc3$last_tick <- pc3$upper - pc3$lower
  
#trying to reformat to get rid of ylim error - yes - need to have underscores not spaces
consensus_tree_pruned$tip.label <- gsub(" ", "_", consensus_tree_pruned$tip.label)
rownames(pc1) <- gsub(" ", "_", rownames(pc1))
#pc1 <- pc1[,c(1:3)]
rownames(pc2) <- gsub(" ", "_", rownames(pc2))
#pc2 <- pc2[,c(1:3)]
rownames(pc3) <- gsub(" ", "_", rownames(pc3))
#pc3 <- pc3[,c(1:3)]

simmap_tree <- readRDS("./analysis/pca_analysis/myc_simmap_trees.rds")

simmap_tree[[1]]$tip.label <- gsub(" ", "_", simmap_tree[[1]]$tip.label)

consensus_tree_pruned$tip.label <- gsub("_", " ", consensus_tree_pruned$tip.label)
myc_named <- setNames(data_spectra$myc, rownames(data_spectra$spectra))

simmap_consensus <- make.simmap(consensus_tree_pruned, myc_named, model="SYM", nsim=100)
summary_simmap <-describe.simmap(simmap_consensus,plot=TRUE,cex=0.7)
plot(summary_simmap)


#plot single map on consensus tree with scores for pc axis 1 and again for pc2
cols_tree<-setNames(c("blue","orange"),unique(myc_named))

colour_bars <- setNames(data_spectra$myc, data_spectra$species)
colour_bars <- gsub("AM", "orange", colour_bars)
colour_bars <- gsub("EM", "blue", colour_bars)
names(colour_bars) <- gsub(" ", "_", names(colour_bars))

jpeg("./output/PCAs/consensus_plots/consensus_myc_pca1_across_reps.jpg", height = 10, width = 6, units = "in", res = 600)
plotTree.wBars(simmap_consensus[[1]],pc1[,c(2)],method="plotSimmap",
               tip.labels=TRUE,fsize=0.7, col = colour_bars, colors = cols_tree)

#plotTree.barplot(simmap_tree[[1]], pc1[,c(3,6,7)], args.barplot=list(beside=FALSE,space=c(0,1.2)))

dev.off()



#get colours for plotting
cols_tree<-setNames(c("blue","orange"),unique(myc_named))

#cols<-brewer.pal(3,"YlGnBu")

ss<-getStates(simmap_consensus[[1]],"tips")
x <- setNames(pc1[,2], rownames(pc1))
colors<-setNames(c("orange","blue"),c("AM","EM"))
barcols<-setNames(sapply(ss,function(x,y) y[which(names(y)==x)],
                         y=colors),names(ss))

#double checking the math
# pc1[,c(1,3,7)]
# 
# pc1[1,3]+pc1[1,2] #1.872009
# 
# (pc1[1,1]-pc1[1,3]) == (pc1[1,7]) #0.1792388
# 
# colnames(pc1)

#barplots for 3 axes - 95% CI
jpeg("./output/PCAs/consensus_plots/consensus_myc_pca123_across_reps_bars_95CI.jpg", height = 10, width = 14, units = "in", res = 600)

par(mfrow=c(1,4)) ## we can also use layout

plotTree.barplot(simmap_consensus[[1]],pc1[,c(3,7)],add=TRUE, args.plotTree = list(fsize= 0.9, xlim = c(0,775)), args.barplot=list(xlab="PC 1", col=c("white", "black"),
                                                   mar=c(5.1,0,2.1,2.1)))

plotTree.barplot(simmap_consensus[[1]],pc2[,c(3,7)],args.barplot=list(xlab="PC 2", col=c("white", "black"),mar=c(5.1,0,2.1,2.1)),args.plotTree=list(plot=FALSE),add=TRUE)
plotTree.barplot(simmap_consensus[[1]],pc3[,c(3,7)],args.barplot=list(xlab="PC 3", col=c("white", "black"),mar=c(5.1,0,2.1,2.1)),args.plotTree=list(plot=FALSE),add=TRUE)

par(mfg=c(1,1))

plot(simmap_consensus[[1]],cols_tree,ftype='off', xlim = c(0,775), mar=c(5.1,1.1,2.1,0)) #method="plotSimmap",
dev.off()
###
#alternative line for setting format of plot
#layout(matrix(c(1,1,1,1,1,2,2,3,3,4,4), nrow = 1, ncol = 11, byrow = TRUE))

#line for adding legend
legend(x="topright",legend=names(colors),pch=22,pt.cex=2,pt.bg=colors,
       box.col="transparent")


#plotting single value coloured by myc type (kind of -  needs fixing)
plotTree.barplot(simmap_tree[[1]],x[simmap_tree[[1]]$tip.label],add=TRUE, args.plotTree = list(fsize= 0.9, xlim = c(0,640)), args.barplot=list(xlab="PC 1",border=barcols, col=c("black"),
                                                                                                                                               mar=c(5.1,0,2.1,2.1)))

#improv boxplot

#get data as named list
pc1_all<-setNames(combo$V1,combo$species_factor)
names(pc1_all) <- gsub(" ", "_", names(pc1_all))

pc2_all<-setNames(combo$V2,combo$species_factor)
names(pc2_all) <- gsub(" ", "_", names(pc2_all))

pc3_all<-setNames(combo$V3,combo$species_factor)
names(pc3_all) <- gsub(" ", "_", names(pc3_all))

#get colour states
ss<-getStates(simmap_consensus[[1]],"tips")
colors<-setNames(c("orange","blue"),c("AM","EM"))

boxcols<-setNames(sapply(ss,function(pc1_all,y) y[which(names(y)==pc1_all)],
                         y=colors),names(ss))

#sort colours to match tree and data
sorted_boxcols <- boxcols[order(match(names(boxcols),simmap_consensus[[1]]$tip.label))]

# #chnage names of tree to match data
# simmap_consensus[[1]]$tip.label <- gsub(" ", "_", simmap_consensus[[1]]$tip.label)


#make plot
jpeg("./output/PCAs/consensus_plots/consensus_myc_pca123_across_reps_boxplot_scores.jpg", height = 10, width = 14, units = "in", res = 600)

#par(mfrow=c(1,4))
layout(matrix(c(1,1,1,2,2,3,3,4,4), nrow = 1, ncol = 9, byrow = TRUE))

plot(simmap_consensus[[1]],cols_tree,ftype='i', fsize = 0.75, xlim = c(0,500), mar=c(5.1,1.1,2.1,0.1))
par(mar=c(5.1,0.1,2.1,1.1))
boxplot(pc1_all~factor(names(pc1_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "")
axis(1)
abline(v = 0)
title(xlab="PC 1 scores")
abline(v = )

boxplot(pc2_all~factor(names(pc2_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "")
axis(1)
abline(v = 0)
title(xlab="PC 2 scores")

boxplot(pc3_all~factor(names(pc3_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "")
axis(1)
abline(v = 0)
title(xlab="PC 3 scores")

dev.off()

#make traitgrams
simmap_consensus[[1]]$tip.label <- gsub(" ", "_", simmap_consensus[[1]]$tip.label)
names(myc_named) <- gsub(" ", "_", names(myc_named))
length(myc_named)

tail(colnames(combo))

#get means for pc1 by name
mean_pc1 <- combo %>% 
  dplyr::select(V1, species_factor) %>% 
  dplyr::group_by(species_factor) %>% 
  dplyr::summarize(mean = mean(V1))

mean_pc1 <- as.data.frame(mean_pc1)
mean_pc1_list <- setNames(mean_pc1$mean, mean_pc1$species_factor)
names(mean_pc1_list) <- gsub(" ", "_", names(mean_pc1_list))

mean_pc2 <- combo %>% 
  dplyr::select(V2, species_factor) %>% 
  dplyr::group_by(species_factor) %>% 
  dplyr::summarize(mean = mean(V2))

mean_pc2 <- as.data.frame(mean_pc2)
mean_pc2_list <- setNames(mean_pc2$mean, mean_pc2$species_factor)
names(mean_pc2_list) <- gsub(" ", "_", names(mean_pc2_list))

mean_pc3 <- combo %>% 
  dplyr::select(V3, species_factor) %>% 
  dplyr::group_by(species_factor) %>% 
  dplyr::summarize(mean = mean(V3))

mean_pc3 <- as.data.frame(mean_pc3)
mean_pc3_list <- setNames(mean_pc3$mean, mean_pc3$species_factor)
names(mean_pc3_list) <- gsub(" ", "_", names(mean_pc3_list))


colors_traitgram<-setNames(c("orange","blue"),c("AM", "EM"))

jpeg("./output/PCAs/consensus_plots/mean_pc1_scores_traitgram_myc_dataset.jpg", height = 8, width = 4, units = "in", res = 600)
par(xaxt="n",yaxt="s")
phenogram(simmap_consensus[[1]], mean_pc1_list, fsize=0.38, spread.labels = TRUE,ftype="i", xlab = "", ylab="Mean PC1 Scores", colors=colors_traitgram)#c("orange", "blue")) spread.cost = c(2,0), 
dev.off()

jpeg("./output/PCAs/consensus_plots/mean_pc2_scores_traitgram_myc_dataset.jpg", height = 8, width = 4, units = "in", res = 600)
par(xaxt="n",yaxt="s")
phenogram(simmap_consensus[[1]], mean_pc2_list, fsize=0.38, spread.labels = TRUE, ftype="i", ylim = c(-2,2), xlab = "", ylab="Mean PC2 Scores", colors=colors_traitgram)#c("orange", "blue"))
dev.off()

jpeg("./output/PCAs/consensus_plots/mean_pc3_scores_traitgram_myc_dataset.jpg", height = 8, width = 4, units = "in", res = 600)
par(xaxt="n",yaxt="s")
phenogram(simmap_consensus[[1]], mean_pc3_list, fsize=0.38, spread.labels = TRUE, ftype="i", ylim=c(-0.75,0.6), ylab="Mean PC3 Scores",xlab = "", colors=colors_traitgram)#c("orange", "blue"))
dev.off()

#plotting in one figure

jpeg("./output/PCAs/consensus_plots/mean_scores_traitgram_myc_dataset_pc123.jpg", height = 8, width = 8, units = "in", res = 800)
pdf("./output/PCAs/consensus_plots/mean_scores_traitgram_myc_dataset_pc123.pdf", height = 8, width = 8)
layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = TRUE))

par(xaxt="n",yaxt="s")
a <- phenogram_edited(simmap_consensus[[1]], mean_pc1_list, fsize=0.45, spread.labels = TRUE, ftype="i", ylab="Mean PC1 Scores", xlab="",#spread.cost = c(0,1),
          colors=colors_traitgram, lwd = 1)#, cex.axis = 0.8)#c("orange", "blue"))
#par(xaxt="n",yaxt="s",font.lab=1)
#axis(1)
#title(xlab="Time since the root",cex.lab=1)
#axis(2)
#title(ylab="Mean PC1 Scores",cex.lab=1.2)

par(xaxt="n",yaxt="s")
b <- phenogram_edited(simmap_consensus[[1]], mean_pc2_list, fsize=0.45, spread.labels = TRUE,  ftype="i", ylab="Mean PC2 Scores", xlab="", ylim=c(-2,2),#spread.cost = c(1,0),
          colors=colors_traitgram, lwd = 1)#c("orange", "blue"))
#par(xaxt="n",yaxt="s",font.lab=1)
#axis(1)
#title(xlab="Time since the root",cex.lab=1)
#axis(2)
#title(ylab="Mean PC2 Scores",cex.lab=1.2)

par(xaxt="n",yaxt="s") #,mar=c(5.1,5.1,2.1,1.1)
c <- phenogram_edited(simmap_consensus[[1]], mean_pc3_list, fsize=0.45, spread.labels = TRUE, ftype="i", ylab="Mean PC3 Scores", xlab="",ylim=c(-0.75,0.6), #spread.cost = c(2,1), 
          colors=colors_traitgram, lwd = 1)#c("orange", "blue"))
#par(xaxt="n",yaxt="s",font.lab=1)
#axis(1)
#title(xlab="Time since the root",cex.lab=1)
#axis(2)
#title(ylab="Mean PC3 Scores",cex.lab=1.2)

dev.off()


#########################

plotTree.wBars(simmap_tree[[1]],pc1[,1], tip.labels=FALSE,fsize=0.7,scale=0.002)

par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
