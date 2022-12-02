#calculating number of samples per species

new_spectra <- readRDS("./data/tidy/new_spectra_matched_trees.rds")

spectra_df <- as.data.frame(new_spectra, row.names = new_spectra$meta$species_names, metadata = TRUE)#, colnames = bands)

colnames(spectra_df)[1:200]
ncol(spectra_df)

#read trees
new_trees_sm <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")
species_names <- new_trees_sm[[1]]$tip.label

#get mean of spectra
counts <- spectra_df[,c(29,166:2166)] %>% 
  dplyr::filter(species_names %in% species_names) %>% 
  group_by(species_names) %>% 
  summarize(n = n())

mean(counts$n)
sd(counts$n)

counts[which(counts$n <= 5),]

spectra_means <- spectra_df[,c(29,166:2166)] %>% 
  #dplyr::select(-c(sample_name)) %>% 
  group_by(species_names) %>% 
  summarize_all(mean, na.rm = TRUE)

#3 format data
spectra_means_df <- data.frame(spectra_means, stringsAsFactors = FALSE)

row.names(spectra_means_df) <- as.character(spectra_means_df$species_names)
spectra_means_df <- spectra_means_df[,-1]
colnames(spectra_means_df) <- gsub("X", "", colnames(spectra_means_df))

spec_list <- as.matrix(spectra_means_df)



#read spectra
spectra <- readRDS("./data/for_analysis/spectra_not_reordered_to_tree.rds")
