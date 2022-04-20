#prepping data for myc models
library(phytools)
library(dplyr)
library(tidyr)

#read trees
trees <- readRDS("./data/for_analysis/final_trees_matched_spectra.rds")

#read spectra
spectra <- readRDS("./data/for_analysis/spectra_not_reordered_to_tree.rds")

#read myc data
data <- read.csv("./data/predictors/myc_data_species_30_03.csv", stringsAsFactors = FALSE)
data <- data[c(1:100),c(1,3)]

colnames(data) <- c("Species", "Myc")

#data$Myc[which(data$Myc == "AM or NM")] <- "AM"
data$Myc[which(data$Myc == "NM")] <- NA

#another version without the multistate species
data_reduced <- data
data_reduced$Myc[which(data_reduced$Myc == "AM-EM")] <- NA
data_reduced$Myc[which(data_reduced$Myc == "ErM")] <- NA


data <- na.omit(data)
data_reduced <- na.omit(data_reduced)
nrow(data)
nrow(data_reduced)

#read gf data
data_2 <- readRDS("./data/tidy/new_growth_forms_matched_spectra.rds")
data_2 <- data_2[,c(2,3,4,5)]

#combine into one growth form column
data_3 <- data_2 %>% 
  pivot_longer(cols = c(2:4)) %>% 
  dplyr::filter(value == 1)

#read lp data
lp <- read.csv("./data/predictors/lp_and_rk.csv", stringsAsFactors = FALSE)
lp <- lp[,c(1,29)]
colnames(lp)[1] <- "species"

lp$Leaf_persistence[which(lp$Leaf_persistence == "deciduous?")] <- "deciduous"

#match data sets
gf_big <- data_3[which(data_3$species_names %in% data$Species),]
lp_big <- lp[which(lp$species %in% data$Species),]
nrow(lp_big)
nrow(gf_big)

unique(data$Species) %in% data_3$species_names

unique(data_3$species_names)

gf_sm <- data_3[which(data_3$species_names %in% data_reduced$Species),]
lp_sm <- lp[which(lp$species %in% data_reduced$Species),]
nrow(gf_sm)
nrow(lp_sm)


#get unique entry per species
data_3 #start here combine shrub and tree etc for those with dups

dups <- gf_big %>% group_by(species_names) %>% dplyr::summarize(n=n()) %>% dplyr::filter(n == 2)

dups_sm <- gf_sm %>% group_by(species_names) %>% dplyr::summarize(n=n()) %>% dplyr::filter(n == 2)

#manually edit these?
gf_big$value[which(gf_big$species_names == "Betula populifolia" & gf_big$name == "Shrub")] <- "0"
gf_big$value[which(gf_big$species_names == "Frangula alnus" & gf_big$name == "Tree")] <- "0"
gf_big$value[which(gf_big$species_names == "Sorbus decora" & gf_big$name == "Shrub")] <- "0"
gf_big$value[which(gf_big$species_names == "Sorbus americana" & gf_big$name == "Shrub")] <- "0"
gf_big$value[which(gf_big$species_names == "Rhamnus cathartica" & gf_big$name == "Shrub")] <- "0"
gf_big$value[which(gf_big$species_names == "Prunus pensylvanica" & gf_big$name == "Shrub")] <- "0"
gf_big$value[which(gf_big$species_names == "Prunus nigra" & gf_big$name == "Shrub")] <- "0"
gf_big$value[which(gf_big$species_names == "Crataegus monogyna" & gf_big$name == "Shrub")] <- "0"
gf_big$value[which(gf_big$species_names == "Acer pensylvanicum" & gf_big$name == "Shrub")] <- "0"
gf_big$value[which(gf_big$species_names == "Oemleria cerasiformis" & gf_big$name == "Tree")] <- "0"

gf_sm$value[which(gf_sm$species_names == "Betula populifolia" & gf_sm$name == "Shrub")] <- "0"
gf_sm$value[which(gf_sm$species_names == "Frangula alnus" & gf_sm$name == "Tree")] <- "0"
gf_sm$value[which(gf_sm$species_names == "Sorbus decora" & gf_sm$name == "Shrub")] <- "0"
gf_sm$value[which(gf_sm$species_names == "Sorbus americana" & gf_sm$name == "Shrub")] <- "0"
gf_sm$value[which(gf_sm$species_names == "Rhamnus cathartica" & gf_sm$name == "Shrub")] <- "0"
gf_sm$value[which(gf_sm$species_names == "Prunus pensylvanica" & gf_sm$name == "Shrub")] <- "0"
gf_sm$value[which(gf_sm$species_names == "Prunus nigra" & gf_sm$name == "Shrub")] <- "0"
gf_sm$value[which(gf_sm$species_names == "Crataegus monogyna" & gf_sm$name == "Shrub")] <- "0"
gf_sm$value[which(gf_sm$species_names == "Acer pensylvanicum" & gf_sm$name == "Shrub")] <- "0"
gf_sm$value[which(gf_sm$species_names == "Oemleria cerasiformis" & gf_sm$name == "Tree")] <- "0"

data_reduced$Species[-which(data_reduced$Species %in% gf_big$species_names)]

#the species missing from the datasets for some reason "Vitis riparia"
#because it's neither herb, shrub or tree

nrow(data_3)

tail(gf_sm)
#combine into one growth form column
data_4 <- gf_big %>% 
  #pivot_longer(cols = c(2:4)) %>% 
  dplyr::filter(value == 1) %>% 
  dplyr::select(species_names, name)

data_4_sm <- gf_sm %>% 
  #pivot_longer(cols = c(2:4)) %>% 
  dplyr::filter(value == 1) %>% 
  dplyr::select(species_names, name)

nrow(data_4)

colnames(data_4) <- c("species_names", "GF")

colnames(data_4_sm) <- c("species_names", "GF")

#prune to match (missing data in myc)
tree_keep <- trees[[1]]
tree_drops <- tree_keep$tip.label[-which(tree_keep$tip.label %in% data$Species)]

tree_drops_sm <- tree_keep$tip.label[-which(tree_keep$tip.label %in% data_reduced$Species)]

new_trees <- lapply(trees,drop.tip,tip=tree_drops)

new_trees_sm <- lapply(trees,drop.tip,tip=tree_drops_sm)
new_trees

tree_sm <- drop.tip(tree_keep, tree_drops)
tree_sm

#prune spectra
spectra_df <- as.data.frame(spectra)
spectra_df_big <- as.matrix(spectra_df[which(rownames(spectra_df) %in% data$Species),])
spectra_df_big <- spectra_df_big[-which(rownames(spectra_df_big) == "Vitis riparia"),]
spectra_df_sm <- as.matrix(spectra_df[which(rownames(spectra_df) %in% data_reduced$Species),])
spectra_df_sm <- spectra_df_sm[-which(rownames(spectra_df_sm) == "Vitis riparia"),]
#spectra_sm <- as_spectra(spectra_df_sm, name_idx = 1, meta_idxs = c(2:8))

#objects to keep
data_4 #gf
tree_sm
spectra_df_sm
lp
data #myc

#drop vitis riparia from myc, lp and species
data[which(data$Species == "Vitis riparia"),]

rownames(spectra_df_big)[-which(rownames(spectra_df_big) %in% data$Species)]

length(rownames(spectra_df_big))
length(data$Species)
length(data_4$species_names)

#I can't figure out which species is missing (98 vs 99 species)

length(rownames(spectra_df_big)) #98 right number
names_final <- rownames(spectra_df_big)

long_names <- data$Species

rownames(spectra_df_big)[-which(rownames(spectra_df_big) %in% data$Species)]

data_4$species_names[which(data_4$species_names %in% data$Species == FALSE)]

data$Species[which(data$Species %in% data_4$species_names == FALSE)]

#trim myc data
data <- data[-which(data$Species %in% "Vitis riparia" == TRUE),]
nrow(data)

data_reduced2 <- data_reduced[-which(data_reduced$Species %in% "Vitis riparia" == TRUE),]
nrow(data_reduced2)

#trim lp data
lp_big$species

lp_big <- lp_big[-which(lp_big$species %in% "Vitis riparia" == TRUE),]
nrow(lp_big)

lp_sm <- lp_sm[-which(lp_sm$species %in% "Vitis riparia" == TRUE),]
nrow(lp_sm)

#trim trees
new_trees_sm <- lapply(new_trees_sm,drop.tip,tip="Vitis riparia")
new_trees <- lapply(new_trees,drop.tip,tip="Vitis riparia")

new_trees

#make list of dataframes and model
#make data list for model
data_spectra <- list(spectra=spectra_df_big, myc = data$Myc, gf = data_4$GF, lp = lp_big$Leaf_persistence, species = data$Species)

data_spectra_mini <- list(spectra=spectra_df_sm, myc = data_reduced2$Myc, gf = data_4_sm$GF, lp = lp_sm$Leaf_persistence, species = data_reduced2$Species)

str(data_spectra_mini)

#save data and trees for myc data

saveRDS(data_spectra, "./data/for_analysis/myc_data_list_98sp_4cat_for_analysis.rds")
saveRDS(data_spectra_mini, "./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

saveRDS(new_trees, "./data/for_analysis/myc_tree_98sp_for_analysis.rds")
saveRDS(new_trees_sm, "./data/for_analysis/myc_tree_92sp_for_analysis.rds")

length(tree_sm)
#98 species - includes Ericoid and AM-EM species - only excludes vitis (non h/t/s) and one NM species
#92 species - only include EM or AM species, excludes AM-EM, ErM, NM and Vitis
