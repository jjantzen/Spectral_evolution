#prepping data for myc models
library(phytools)
library(dplyr)
library(tidyr)
library(spectrolab)

#read data and trees for myc data

#data_spectra <- readRDS("./data/for_analysis/myc_data_list_98sp_4cat_for_analysis.rds")
data_spectra_mini <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

#new_trees <- readRDS("./data/for_analysis/myc_tree_98sp_for_analysis.rds")
new_trees_sm <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

#get list of taxa to exclude: conifers and fern
new_trees_sm$STATE_241090000$tip.label

seed_plants <- c("Polystichum munitum")

ang_only <- c("Abies balsamea", "Larix laricina", "Picea abies", "Picea glauca", "Picea mariana", "Pinus banksiana", "Pinus resinosa", "Pinus rigida", "Pinus strobus", 
              "Polystichum munitum", "Thuja occidentalis", "Tsuga canadensis")

#drop tips and taxa from data
str(data_spectra_mini)

#excluding fern
seed_data <- data.frame(data_spectra_mini, stringsAsFactors = FALSE) %>% 
  dplyr::filter(species != seed_plants)
                
tail(colnames(seed_data))
nrow(seed_data)

seed_data_spectra <- seed_data[,c(1:2001)]
colnames(seed_data_spectra) <- gsub("spectra.", "", colnames(seed_data_spectra))
rownames(seed_data_spectra) <- seed_data_species

str(seed_data_spectra)

ncol(seed_data_spectra)

seed_data_myc <- seed_data[,2002]
seed_data_gf <- seed_data[,2003]
seed_data_lp <- seed_data[,2004]
seed_data_species <- seed_data[,2005]

seed_data_spectra_list <- as.matrix(seed_data_spectra)

str(seed_data_spectra_list)

data_spectra_seeds <- list(spectra=seed_data_spectra_list, myc = seed_data_myc, gf = seed_data_gf, lp = seed_data_lp, species = seed_data_species)

str(data_spectra_seeds)

seed_trees <- new_trees <- lapply(new_trees_sm,drop.tip,tip=seed_plants)

#only ang
ang_data <- as.data.frame(data_spectra_mini, stringsAsFactors = FALSE) %>% 
  dplyr::filter(!species %in% ang_only)

ang_data_spectra <- ang_data[,c(1:2001)]
colnames(ang_data_spectra) <- gsub("spectra.", "", colnames(ang_data_spectra))
rownames(ang_data_spectra) <- ang_data_species


ang_data_myc <- ang_data[,2002]
ang_data_gf <- ang_data[,2003]
ang_data_lp <- ang_data[,2004]
ang_data_species <- ang_data[,2005]

ang_data_spectra_list <- as.matrix(ang_data_spectra)
str(ang_data_spectra_list)

data_spectra_ang <- list(spectra=ang_data_spectra_list, myc = ang_data_myc, gf = ang_data_gf, lp = ang_data_lp, species = ang_data_species)

str(data_spectra_ang)

ang_trees <- new_trees <- lapply(new_trees_sm,drop.tip,tip=ang_only)

#save output for running on cluster
saveRDS(ang_trees, "./data/for_analysis/ang_only_trees_for_myc.rds")
saveRDS(seed_trees, "./data/for_analysis/seed_plant_trees_for_myc.rds")

saveRDS(data_spectra_seeds, "./data/for_analysis/seed_plant_data_for_myc.rds")
saveRDS(data_spectra_ang, "./data/for_analysis/ang_only_data_for_myc.rds")

str(data_spectra_seeds)
