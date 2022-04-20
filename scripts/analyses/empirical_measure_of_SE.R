#calculate SE and variance of spectra for means

library(mvMORPH)
library(spectrolab)
library(dplyr)

#read spectra
new_spectra <- readRDS("./data/tidy/new_spectra_matched_trees.rds")

#make into dataframe
spectra_df <- as.data.frame(new_spectra, row.names = new_spectra$meta$species_names, metadata = TRUE)#, colnames = bands)

#get mean of spectra
spectra_se <- spectra_df[,c(29,166:2166)] %>% 
  #dplyr::select(-c(sample_name)) %>% 
  group_by(species_names) %>% 
  summarise_all(list(sd = ~sd(.), se = ~sd(./sqrt(.)))) %>% 
  dplyr::select(-contains("sd"))#mean = ~mean(.), 
  #summarize_all(mean, na.rm = TRUE) 

colnames(spectra_means)[2003:4003]
spectra_means$species_names

#calculate mean se across species
spectra_mean_se <- spectra_se %>% 
  dplyr::select(-species_names) %>% 
  summarise_all(list(mean), na.rm = TRUE)
 #summarise(across(everything(), list(mean)))

spectra_max_se <- spectra_se %>% 
  dplyr::select(-species_names) %>% 
  summarise_all(list(max), na.rm = TRUE)

max(spectra_max_se) #0.0682

#summarise(across(everything(), list(mean)))

#summarise(across(everything(), list(mean)))


saveRDS(spectra_mean_se, "./data/empirical_spectral_stats/spectra_standard_error.rds")

rowMeans(spectra_mean_se) #0.01520941
