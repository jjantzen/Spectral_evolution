#brightness normalization code
library(spectrolab)
library(dplyr)
library(tidyr)

#load spectral data
pre_spectra <- readRDS("./data/for_analysis/final_spectra_matched_trees.rds")

#convert to df
spectra_df <- as.data.frame(pre_spectra)
colnames(spectra_df)

spectra_df <- spectra_df[,c(1,29,166:2166)]

str(spectra_df)

spectra_long <- pivot_longer(spectra_df, cols = -c(1:2))

# normalise all spectra
### Brightness normalization ----
bright_norm <- function(x) {
  x$value_norm <- x$value / sqrt(sum(x$value^2))
  return(x)  
}


# apply normalization function
spectra_avg_plant_norm <- spectra_long %>% 
  group_by(sample_name) %>% 
  do(bright_norm(.))

saveRDS(spectra_avg_plant_norm, "./data/for_analysis/brightness_normalized_sample_all_taxa.rds")

str(spectra_avg_plant_norm)

str(spectra_avg_plant_norm$value_norm)

spectra_avg_plant_norm_df <- as.data.frame(spectra_avg_plant_norm)

#spread normalized data
wide_normalized <- spectra_avg_plant_norm_df %>% 
  dplyr::select(-c(value, value_norm)) %>% 
  pivot_wider(id_cols = c(sample_name, species_names), names_from = name, values_from = value_norm_num)

wide_normalized

#get mean of spectra
spectra_means <- wide_normalized%>% 
  dplyr::select(-c(sample_name)) %>% 
  dplyr::group_by(species_names) %>% 
  dplyr::summarize_all(mean, na.rm = TRUE)

saveRDS(spectra_means, "./data/for_analysis/brightness_normalized_species_means_all_taxa.rds")

#compare with other formats
spectra_means

spectra_means_df <- data.frame(spectra_means, stringsAsFactors = FALSE)

row.names(spectra_means_df) <- as.character(spectra_means_df$species_names)
spectra_means_df <- spectra_means_df[,-1]
colnames(spectra_means_df) <- gsub("X", "", colnames(spectra_means_df))

spec_list <- as.matrix(spectra_means_df)
str(spec_list)

saveRDS(spec_list, "./data/for_analysis/brightness_normalized_spectra_not_reordered.rds")


#make new data objects for myc

#read in final objects

myc_data_list_92sp_binary_for_analysis <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

str(myc_data_list_92sp_binary_for_analysis)

spec_df <- as.data.frame(spec_list)
spectra_df_big <- as.matrix(spec_df[which(rownames(spec_df) %in% myc_bright_norm$species),])

myc_bright_norm <- list(spectra = spectra_df_big,
                        myc = myc_data_list_92sp_binary_for_analysis$myc, gf = myc_data_list_92sp_binary_for_analysis$gf, lp = myc_data_list_92sp_binary_for_analysis$lp, species = myc_data_list_92sp_binary_for_analysis$species)

saveRDS(myc_bright_norm, "./data/for_analysis/brightness_normalized_myc_data_list_92sp.rds")

#other data versions
final_spectra_matched_trees

myc_data_list_for_analysis

str(spectra_not_reordered_to_tree)
