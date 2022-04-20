#Figure out which taxa missing from shade, drought and soil trait data

library(phytools)

#big trees
new_trees <- readRDS("./data/tidy/new_trees_matched_spectra.rds")

#data
trait_data <- readRDS("./data/for_analysis/reduced_trait_data.rds")

trait_data

#species missing

missing_species <- new_trees[[1]]$tip.label[-which(new_trees[[1]]$tip.label %in% trait_data$Species)]

missing_species

#original metadata to figure out where species are from
metadata <- readRDS("./data/tidy/new_spectra_matched_trees.rds")

#matching dataset
metadata_df <- as.data.frame(metadata)
colnames(metadata_df_matched)

metadata_df_matched <- metadata_df[which(metadata_df$species_names %in% missing_species),]
metadata_matched <- as_spectra(metadata_df_matched, name_idx = 1, meta_idxs = c(2:165))

colnames(metadata_df_matched)

#keep just species and projects?
project_info <- metadata_df_matched[,c(1,29,89,90)]

#group by species 
species_by_project <- project_info %>% 
  dplyr::select(-sample_name) %>% 
  dplyr::distinct(species_names, project) %>% 
  dplyr::arrange(species_names)

species_by_site <- project_info %>% 
  dplyr::select(-sample_name) %>% 
  dplyr::distinct(species_names, site) %>% 
  dplyr::arrange(species_names)

species_by_project
species_by_site

#getting lists of species by project
CGOP_species <- species_by_project$species_names[which(species_by_project$project == "2018-Hacker-PhD-UBC")]

nonCGOP_species <- species_by_project$species_names[-which(species_by_project$project == "2018-Hacker-PhD-UBC")]

CGOP_species
nonCGOP_species

nonCGOP_species %in% CGOP_species

#getting lists of species by sites
CGOP_site_species <- species_by_site$species_names[which(species_by_site$site == "CGOP_1")]

nonCGOP_site_species <- species_by_site$species_names[-which(species_by_site$site == "CGOP_1")]

nonCGOP_site_species %in% CGOP_site_species

####

