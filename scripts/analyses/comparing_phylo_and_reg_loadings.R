#plot loadings of phylo PC vs regular PC
library(phytools)
library(ggplot2)
library(mvMORPH)
library(dplyr)

#read in data
data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")
tree_myc <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

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

spectra.pca <- prcomp(data_spectra$trait, center = TRUE, scale. = TRUE)





pc_labs <- paste("PC", c(1:9))
names(pc_labs) <- c(1:9)

jpeg(paste0("./output/PCAs/iterations/myc_loadings_top10_iteration", i, ".jpg"), res = 400, width = 15, height = 15, units = "in")
v <- ggplot(vectors_for_plotting[which(vectors_for_plotting$pc_axis %in% c(1:9)),], aes(x=wavelength,y=value)) +
  geom_bar(stat = "identity")+
  facet_wrap(vars(pc_axis), ncol = 3, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
  xlab("Wavelength (nm)")+
  ylab("")+ #Vector unit
  theme_bw()+
  theme(strip.text.x = element_text(size = 15), axis.title=element_text(size=15))
print(v)
dev.off()


#prep tree if needed
consensus_tree <- readRDS("./data/for_analysis/myc_consensus_tree.rds")

#reorder species by tree order 
combo$species_factor <- factor(combo$species,
                               levels = consensus_tree_pruned$tip.label)


