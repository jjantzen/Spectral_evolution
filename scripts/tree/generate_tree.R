#Generate species tree from tol for spectral data

#load libraries
library(dplyr)
library(ape)
library(taxize)
library(phytools)
library(rotl)

#read in species list
species_list <- read.csv("./data/tidy/species_names_for_phylo.csv", stringsAsFactors = FALSE)

#make list into vector of characters
species_binomials <- as.vector(species_list$species_name)

species_binomials_clean <- gsub("_", " ", species_binomials)

#resolve taxonomy - have to remove Kalmia var and Phragmites subsp 
species_names_resolved <- tnrs_match_names(species_binomials_clean[-c(55,63)])

length(species_binomials)
length(species_binomials_clean)
nrow(species_names_resolved)

species_names_resolved$unique_name[which(species_names_resolved$is_synonym == TRUE)] %in% species_binomials_clean

#check which ones missing
length(species_binomials_clean)
nrow(species_names_resolved)
colnames(species_names_resolved)

#get OTT IDs for those missing
missing_taxa <- species_binomials_clean[c(55,63)]

missing_taxa <- gsub("var. angustifolia", "", missing_taxa)
missing_taxa <- gsub("subsp. australis", "", missing_taxa)
missing_taxa

#try again with missing ones
missing_ids <- tnrs_match_names(missing_taxa)
species_names_resolved

#join dataframes
joint_resolved_names <- rbind(species_names_resolved, missing_ids)
nrow(joint_resolved_names)

write.csv(joint_resolved_names, "./data/tidy/resolved_taxonomy.csv", row.names = FALSE)

#some flags: incertae_sedis 
#what are sibling_higher?

#check if ott ids will work
notin_tree <- joint_resolved_names$search_string[which(is_in_tree(joint_resolved_names$ott_id) == FALSE)]
#joint_resolved_names[160:173,]

#get tree from OTOL
spectra_tree <- tol_induced_subtree(ott_ids = joint_resolved_names$ott_id)
spectra_tree

#plot phylogeny
plot(spectra_tree, cex=0.5)

write.tree(spectra_tree, "./data/tol_tree.tre")

#add branch lengths
spectra_tree_br <- compute.brlen(spectra_tree, method = "Grafen", power = 1)

#make ultrametric
spectra_tree_ultra <- chronos(spectra_tree_br, lambda=0)  

plot(spectra_tree_ultra, cex = 0.5)

write.tree(spectra_tree_ultra, "./data/tol_brlen_grafen1_lam0.tre")
