#prepping exploratory data

data <- readRDS("./data/for_analysis/final_data.rds")
data

trees <- readRDS("./data/for_analysis/final_trees_matched_spectra.rds")
trees

reduced_data <- readRDS("./data/for_analysis/reduced_trait_data.rds")
reduced_data

tree_drops <- trees[[1]]$tip.label[-which(trees[[1]]$tip.label %in% reduced_data$Species)]

reduced_trees <- lapply(trees,drop.tip,tip=tree_drops)
reduced_trees

#saveRDS(reduced_trees, "./data/for_analysis/reduced_trees_exploratory.rds")

data_drops <- reduced_data$Species[-which(reduced_data$Species %in% trees[[1]]$tip.label)]

trees[[1]]
nrow(reduced_data)

spectra <- readRDS("./data/for_analysis/spectra_not_reordered_to_tree.rds")

spectra_df <- as.data.frame(spectra)
spectra_df_big <- as.matrix(spectra_df[which(rownames(spectra_df) %in% reduced_data$Species),])

#combine trait and spectra
data_spectra <- list(spectra=spectra_df_big, woody = reduced_data$Woody, shade = reduced_data$Shade, coarse = reduced_data$Coarse_soil, fine = reduced_data$Fine_soil, lp = reduced_data$leaf_persistence, drought = reduced_data$Drought_bin, species = reduced_data$Species)



