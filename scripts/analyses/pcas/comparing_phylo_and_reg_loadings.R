#plot loadings of phylo PC vs regular PC
library(phytools)
library(ggplot2)
library(mvMORPH)
library(dplyr)
library(ggpubr)
library(tidyr)
library(gtable)    
library(grid)
library(gridExtra) 

#read in data
data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")
tree_myc <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

consensus_tree <- readRDS("./data/for_analysis/myc_consensus_tree.rds")

tips_to_drop <- consensus_tree$tip.label[-which(consensus_tree$tip.label %in% data_spectra$species)]

consensus_tree_pruned <- drop.tip(consensus_tree, tip = tips_to_drop)

#saveRDS(model_list, "./analysis/pca_analysis/best_intercept_models_for_pcas_92sp.rds")
model_list <- readRDS("./analysis/pca_analysis/best_intercept_models_for_pcas_92sp.rds")


#get pca scores for all replicates into one dataframe
for (i in 11:100){#0){#length(model_list)){
  #get pca
  print(i)
  pca_best_model <- mvgls.pca(model_list[[i]], plot = FALSE)  
  
  vectors_df <- pca_best_model$vectors %>%
    as.data.frame()
  
  vectors_df$wavelength <- c(400:2400)
  vectors_for_plotting <- vectors_df %>% tidyr::gather(., pc_axis, value, -wavelength)
  vectors_for_plotting$pc_axis <- gsub("V", "", vectors_for_plotting$pc_axis)
  class(vectors_for_plotting$pc_axis) <- "numeric"
  
  #structure data
  data_from_pcas <- vectors_for_plotting %>% 
    dplyr::filter(pc_axis %in% c(1:6))
  
  data_from_pcas$iteration <- i
  
  #combine data
  if (i == 1){
    combo <- data_from_pcas
  } else {
    combo <- rbind(combo, data_from_pcas)
  } 
}


#saveRDS(combo, "./analysis/pca_analysis/vectors_all_phylo_pc1to6.rds")

combo <- readRDS("./analysis/pca_analysis/vectors_all_phylo_pc1to6.rds")

spectra.pca <- prcomp(data_spectra$spectra, center = TRUE, scale. = TRUE)
prcomp_keeprs <- as.data.frame(spectra.pca$rotation[,c(1:6)], row.names = dimnames(spectra.pca$rotation)[1])

str(prcomp_keeprs)

prcomp_keeprs$wavelength <- rownames(prcomp_keeprs)
str(prcomp_keeprs)

prcomp_keeprs$iteration <- "0"

prcomp_long <- pivot_longer(prcomp_keeprs, cols = c(1:6))
colnames(prcomp_long)[3] <- "pc_axis"

prcomp_long$pc_axis <- gsub("PC", "", prcomp_long$pc_axis)
class(prcomp_long$pc_axis) <- "numeric"
class(prcomp_long$wavelength) <- "integer"
class(prcomp_long$iteration) <- "integer"

#combine non phylo with phylo 
combo_2 <- rbind(combo, prcomp_long)

#saveRDS(combo_2, "./analysis/pca_analysis/vectors_all_phylo_pc1to6_inc_nonphylo.rds")




#make labels for plotting
pc_labs <- paste("PC", c(1:6))
names(pc_labs) <- c(1:6)

jpeg(paste0("./output/PCAs/consensus_plots/myc_loadings_top10_all_iteration.jpg"), res = 400, width = 15, height = 15, units = "in")
v <- ggplot(combo_2[which(combo_2$iteration == 0),], aes(x=wavelength,y=value)) +
  geom_bar(stat = "identity")+
  facet_wrap(vars(pc_axis), ncol = 3, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
  xlab("Wavelength (nm)")+
  ylab("")+ #Vector unit
  theme_bw()+
  theme(strip.text.x = element_text(size = 15), axis.title=element_text(size=15))

v

jpeg(paste0("./output/PCAs/consensus_plots/myc_loadings_top10_all_iteration.jpg"), res = 400, width = 15, height = 15, units = "in")
#make boxplot for variation in loadings across iterations
ggplot(combo_2[which(combo_2$iteration != 0),]) +
  geom_boxplot(aes(x = as.factor(wavelength), y = value))+
  facet_wrap(vars(pc_axis), ncol = 3, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
  xlab("Wavelength (nm)")+
  ylab("")+ #Vector unit
  geom_hline(yintercept = 0)+
  theme(strip.text.x = element_text(size = 15), axis.title=element_text(size=15), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  theme_bw()

ggplot(combo_2[which(combo_2$iteration == 0),], aes(x = as.factor(wavelength), y = value)) +
  geom_bar(stat = "identity")+
  facet_wrap(vars(pc_axis), ncol = 3, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
  xlab("Wavelength (nm)")+
  ylab("")+ #Vector unit
  geom_hline(yintercept = 0)+
  theme(strip.text.x = element_text(size = 15), axis.title=element_text(size=15), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  theme_bw()
dev.off()


#plot both together

combo_subset_0 <- subset(combo_2, iteration == "0")

combo_subset_1 <- subset(combo_2, iteration != "0")

combo_subset_pc1 <- subset(combo_subset_1, pc_axis == "1")
keep_pc1_iterations <- combo_subset_pc1$iteration[which(combo_subset_pc1$wavelength == 1000 & combo_subset_pc1$value <= 0)]
keep_pc1_combo <- combo_subset_pc1[which(combo_subset_pc1$iteration %in% keep_pc1_iterations),]

combo_subset_pc2 <- subset(combo_subset_1, pc_axis == "2")
keep_pc2_iterations <- combo_subset_pc2$iteration[which(combo_subset_pc2$wavelength == 500 & combo_subset_pc2$value >= 0)]
keep_pc2_combo <- combo_subset_pc2[which(combo_subset_pc2$iteration %in% keep_pc2_iterations),]

#combo_subset_pc2 <- subset(combo_subset_pc2, value >= "0")

combo_subset_pc3 <- subset(combo_subset_1, pc_axis == "3")
keep_pc3_iterations <- combo_subset_pc3$iteration[which(combo_subset_pc3$wavelength == 1000 & combo_subset_pc3$value >= 0)]
keep_pc3_combo <- combo_subset_pc3[which(combo_subset_pc3$iteration %in% keep_pc3_iterations),]

#combo_subset_pc3 <- subset(combo_subset_pc3, value <= "0")

combo_subset_pc4 <- subset(combo_subset_1, pc_axis == "4")
keep_pc4_iterations <- combo_subset_pc4$iteration[which(combo_subset_pc4$wavelength == 1000 & combo_subset_pc4$value >= 0)]
keep_pc4_combo <- combo_subset_pc4[which(combo_subset_pc4$iteration %in% keep_pc4_iterations),]

#combo_subset_pc4 <- subset(combo_subset_pc4, value >= "0")

combo_subset_pc5 <- subset(combo_subset_1, pc_axis == "5")
keep_pc5_iterations <- combo_subset_pc5$iteration[which(combo_subset_pc5$wavelength == 1000 & combo_subset_pc5$value <= 0)]
keep_pc5_combo <- combo_subset_pc5[which(combo_subset_pc5$iteration %in% keep_pc5_iterations),]

#combo_subset_pc5 <- subset(combo_subset_pc5, value >= "0")

combo_subset_pc6 <- subset(combo_subset_1, pc_axis == "6")
keep_pc6_iterations <- combo_subset_pc6$iteration[which(combo_subset_pc6$wavelength == 1000 & combo_subset_pc6$value >= 0)]
keep_pc6_combo <- combo_subset_pc6[which(combo_subset_pc6$iteration %in% keep_pc6_iterations),]

#combo_subset_pc6 <- subset(combo_subset_pc6, value <= "0")

#combine into one dataframe again for plotting
all_kept_pc1_to_6 <- rbind(keep_pc1_combo, keep_pc2_combo, keep_pc3_combo, keep_pc4_combo, keep_pc5_combo, keep_pc6_combo)

#edit aesthetics - axis ticks/labels, size of facet labels, colour of background
#jpeg(paste0("./output/PCAs/consensus_plots/myc_loadings_top6_all_iteration_with_nonphylo.jpg"), res = 400, width = 20, height = 15, units = "in")
jpeg(paste0("./output/PCAs/consensus_plots/myc_loadings_top6_single_direction_with_nonphylo.jpg"), res = 400, width = 20, height = 15, units = "in")

p <- ggplot(combo_subset_0, aes(x = as.factor(wavelength), y = value)) +
  geom_bar(stat = "identity", colour = "blue")+
  facet_wrap(vars(pc_axis), ncol = 6, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
  xlab("Wavelength (nm)")+
  ylab("")+ #Vector unit
  geom_hline(yintercept = 0)+
  theme(strip.text.x = element_text(size = 15), axis.title=element_text(size=15), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  theme_bw()


q <- ggplot(all_kept_pc1_to_6) +
  geom_point(aes(x = as.factor(wavelength), y = value), alpha = 0.2)+
  facet_wrap(vars(pc_axis), ncol = 6, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
  xlab("Wavelength (nm)")+
  ylab("")+ #Vector unit
  geom_hline(yintercept = 0)+
  theme_bw()+
  theme(strip.text.x = element_text(size = 15), axis.title=element_text(size=15), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

ggarrange(p,q, nrow = 2)
  

dev.off()



ggplot() +
  geom_boxplot(combo_2, aes(x = as.factor(wavelength), y = value))+
  geom_bar(stat = "identity")+
  facet_wrap(vars(pc_axis), ncol = 3, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
  xlab("Wavelength (nm)")+
  ylab("")+ #Vector unit
  geom_hline(yintercept = 0)+
  theme(strip.text.x = element_text(size = 15), axis.title=element_text(size=15), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  theme_bw()

#prep tree if needed
consensus_tree <- readRDS("./data/for_analysis/myc_consensus_tree.rds")

#reorder species by tree order 
combo$species_factor <- factor(combo$species,
                               levels = consensus_tree_pruned$tip.label)

########
#plotting loadings on pca ouwie plots

aics_ouwie <- readRDS("./analysis/pca_analysis/pc_univariate_model_aics_92sp.rds")
aovs_univariate <- readRDS("./analysis/pca_analysis/pc_univariate_model_aovs_92sp.rds")

#plot the aic weights for each model for each pc axis (repeated iterations)
long_aics <- pivot_longer(aics_ouwie, c(3:9))
long_aics

#write function for labeling facets
pc_names <- as_labeller(c(`1` = "PC 1", `2` = "PC 2", `3` = "PC 3", `4` = "PC 4", `5` = "PC 5",
                          `6` = "PC 6", `7` = "PC 7", `8` = "PC 8", `9` = "PC 9", `10` = "PC 10"))

p <- ggplot(long_aics[which(long_aics$PC_axis %in% c(1:3)),], aes(x = name, y = round(value, 2)))+
  #geom_point()+
  geom_boxplot(aes(fill = as.factor(name)))+
  #geom_jitter(color = "black", size = 0.5)+
  labs(y = "AICw", x = "Model")+
  facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  scale_fill_discrete(name = "Model")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 12), axis.title=element_text(size=12), axis.text = element_text(size = 6),
        strip.background = element_rect(fill = "white", colour = "grey50"), legend.title = element_text(size=12), legend.text = element_text(size=10))
  #theme_bw()

pc_labs <- paste("PC", c(1:3))
names(pc_labs) <- c(1:3)

q <- ggplot(combo[which(combo$pc_axis %in% c(1:3)),], aes(x = wavelength, y = value)) + #as.factor()
  geom_bar(stat = "identity", aes(fill = as.factor(pc_axis)))+#, show.legend = FALSE
  scale_x_continuous(breaks = c(400, 900, 1400, 1900, 2400), labels = c("400", "900", "1400", "1900", "2400"))+
  facet_wrap(vars(pc_axis), ncol = 3, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
  xlab("Wavelength (nm)")+
  ylab("PC Loadings")+ #Vector unit
  scale_fill_manual(name = "Axis", values = c("#619CFF", "#619CFF", "#619CFF"))+
  geom_hline(yintercept = 0)+
  theme(strip.text.x = element_text(size = 12), strip.background = element_rect(fill = "white", colour = "grey50"), axis.title=element_text(size=12), #axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey50"), legend.key = element_rect(colour = "grey80"), axis.text = element_text(size = 6),
        legend.title = element_text(size=12), legend.text = element_text(size=10), panel.background = element_rect(fill = "white", colour = NA))
q
# # Get the gtables
# gA <- ggplotGrob(p)
# gB <- ggplotGrob(q)
# 
# # Set the widths
# gA$widths <- gB$widths
# 
# # Arrange the two charts.
# # The legend boxes are centered
# grid.newpage()
# grid.arrange(gA, gB, nrow = 2)

jpeg("./output/PCAs/ouwie/boxplot_facet_by_pc_model_with_loadings.jpg", width = 8, height = 6, units = "in", res = 600)
#pdf("./output/PCAs/ouwie/boxplot_facet_by_pc_model_with_loadings.pdf", width = 8, height = 6)
ggarrange(p,q, ncol = 1, common.legend = TRUE, legend = "right", labels = c("A", "B"))# = 2, widths = c(0.5, 1))
dev.off()

