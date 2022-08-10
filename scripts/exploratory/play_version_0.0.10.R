##play with data in version 0.0.10 of spectrolab
library(spectrolab)
library(rotl)
library(dplyr)
library(ape)
library(taxize)
library(phytools)
library(mvMORPH)

##check the version
packageVersion("spectrolab")

##because saved from dataframe, rownames in column 1
spec_csv <- read.csv("./data/raw/Shan_spectral_data_saved.csv", check.names = FALSE, row.names = 1)

#can read in rds now 
CABO_spectra <- readRDS("./data/all_spectra.rds")

##checking the right columns as metadata
colnames(spec_csv)
#spec_csv[,1]

##make into class spectra with names and metadata columns specified
data_spectral <- as_spectra(spec_csv, name_idx = 1, meta_idxs = c(2:4))

str(data_spectral)
dim(data_spectral)
str(CABO_spectra)
dim(CABO_spectra)

##get metadata separately
metadata <- meta(CABO_spectra)

n <- names(CABO_spectra)

w <- bands(CABO_spectra)

r <- value(CABO_spectra)

s <- meta(CABO_spectra, "species", simplify = TRUE)

##subset spectra examples

#just bands labeled with integers
CABO_spectra[,800:900]

#all bands within the range even non-integers
CABO_spectra[,bands(CABO_spectra, 800, 900)]

#get just a couple individuals for that range
eg1 <- CABO_spectra[1:3,400:2400]

#plot spectra with lines coloured by species
plot(eg1)
plot(CABO_spectra)
plot_quantile(CABO_spectra, total_prob = 0.95, col = rgb(1,0,0,0.5), lwd = 0.5, border = TRUE, add = TRUE)
plot_regions(CABO_spectra, regions = default_spec_regions(), add = TRUE)
title("90% spectral quantile")

plot_interactive(CABO_spectra)

#subset by metadata using split function

unique(metadata$project)
head(metadata)

16340856 %in% metadata$sample_id

metadata[which(metadata$sample_id == 16340856),]

#assign SWA-Warren species IDs to Agonis flexuosa (Willd.) Sweet
metadata$species[which(metadata$project == "SWA-Warren")]

#first add new name as factor level
metadata$species <- factor(metadata$species, levels=c(levels(metadata$species), "Agonis flexuosa (Willd.) Sweet"))
#then change SWA-Warren values (NAs) to right name - just in metadata object not in spectral object
metadata$species[which(metadata$project == "SWA-Warren")] <- "Agonis flexuosa (Willd.) Sweet"

#change in spectral object too
meta(CABO_spectra)$species <- factor(meta(CABO_spectra)$species, levels=c(levels(meta(CABO_spectra)$species), "Agonis flexuosa (Willd.) Sweet"))
meta(CABO_spectra)$species[which(meta(CABO_spectra)$project == "SWA-Warren")] <- "Agonis flexuosa (Willd.) Sweet"
tail(meta(CABO_spectra))


metadata[which(metadata$project == "2017-Dessain-MSc"),]

#still need to add species IDs for these two projects
metadata[which(metadata$project == "2019-Pardo-MSc-UdeM"),]
metadata[which(metadata$project == "2019-Phragmites-temporal"),]


#split spectra by species
spec_list <- split(CABO_spectra, as.character(meta(CABO_spectra)$species)) #why does this need to be as character?
spec_list <- split(data_spectral, meta(data_spectral)$species)

unique(meta(data_spectral)$species)
unique(meta(CABO_spectra)$species)


length(CABO_spectra)
plot(spec_list$`Abies balsamea (Linnaeus) Miller`) #, col = c("red", "green", "blue", "yellow", "purple", "black")

length(spec_list)

saveRDS(spec_list, "./data/spectra_by_species.rds")

#alternative way to subset - looks like 4 samples have duplicates?
spec_subset_dups <- subset_by(CABO_spectra, by = names(CABO_spectra), n_min = 2, n_max = Inf)


jpeg("./figures/website_plot.jpg", width = 600, height = 400, units = "px")
plot(mean(CABO_spectra[which(meta(CABO_spectra, "species", simplify = TRUE) == "Populus balsamifera Linnaeus"),]), col = "purple", xlab = "Wavelength (nm)", ylab = "Reflectance")
plot(mean(CABO_spectra[which(meta(CABO_spectra, "species", simplify = TRUE) == "Prunus serotina Ehrhart"),]), col = "red", add = TRUE)
plot(CABO_spectra[which(meta(CABO_spectra, "species", simplify = TRUE) == "Claytonia perfoliata Donn ex Willdenow"),], col = "gold", add = TRUE)
plot(mean(CABO_spectra[which(meta(CABO_spectra, "species", simplify = TRUE) == "Alnus incana subsp. rugosa (Du Roi) R.T. Clausen"),]), col = "blue", add = TRUE)
plot(mean(CABO_spectra[which(meta(CABO_spectra, "species", simplify = TRUE) == "Picea rubens Sargent"),]), col = "green", add = TRUE)
plot(mean(CABO_spectra[which(meta(CABO_spectra, "species", simplify = TRUE) == "Cornus sericea Linnaeus"),]), col = "brown", add = TRUE)
plot(mean(CABO_spectra[which(meta(CABO_spectra, "species", simplify = TRUE) == "Bromus sterilis Linnaeus"),]), col = "dark green", add = TRUE)
plot(mean(CABO_spectra[which(meta(CABO_spectra, "species", simplify = TRUE) == "Camassia quamash (Pursh) Greene"),]), col = "orange", add = TRUE)
dev.off()


unique(meta(CABO_spectra, "species", simplify = TRUE))
spec_subset_singles <- subset_by(CABO_spectra, by = names(CABO_spectra), n_min = 1, n_max = 1)

plot(spec_subset_singles)
#spectra 449-451 are odd also 1392

odd_vis <- CABO_spectra[449:451,]

meta(CABO_spectra)[c(449:451,1392),]

unique(meta(CABO_spectra)$species)

no_sp_ids_spectra <- CABO_spectra[which(is.na(meta(CABO_spectra)$species))]

all_sp_ids_spectra <- CABO_spectra[which(!is.na(meta(CABO_spectra)$species))]

#why do these ones not have species IDs?
meta(no_sp_ids_spectra)

species_names <- as.character(unique(meta(all_sp_ids_spectra)$species))

split_names <- strsplit(species_names, " ")

species_binomials <- list()

#trying to include subspecific taxa
for (i in 1:length(split_names)){
  if (i == 28){
    species_binomials[i] <- paste0(split_names[[i]][[1]])  
  } else if (i == 42){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[3]], " ", split_names[[i]][[4]])  
  } else if (i == 57){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[7]], " ", split_names[[i]][[8]])  
  } else if (i == 68){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[3]], " ", split_names[[i]][[4]])
  } else if (i == 87){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[3]], " ", split_names[[i]][[4]])
  } else if (i == 95){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[3]], " ", split_names[[i]][[4]])
  } else if (i == 96){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[3]], " ", split_names[[i]][[4]], " ", split_names[[i]][[5]], " ", split_names[[i]][[6]])
  } else if (i == 97){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[4]], " ", split_names[[i]][[5]])
  } else if (i == 98){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[4]], " ", split_names[[i]][[5]])
  } else if (i == 99){
    species_binomials[i] <- paste0(split_names[[i]][[1]])
  } else if (i == 130){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[4]], " ", split_names[[i]][[5]])
  } else if (i == 131){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[3]], " ", split_names[[i]][[4]])
  } else if (i == 132){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[4]], " ", split_names[[i]][[5]])
  } else if (i == 141){
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]], " ", split_names[[i]][[3]], " ", split_names[[i]][[4]])
  } else {
    species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]])
  }
}

#no subspecific taxa included
for (i in 1:length(split_names)){
  species_binomials[i] <- paste0(split_names[[i]][[1]], " ", split_names[[i]][[2]])
}

#make list into vector of characters
species_binomials <- as.vector(as.character(species_binomials))

#resolve taxonomy (not at subspecific level) - only with subspecies not included
species_names_resolved <- tnrs_match_names(species_binomials)

#synonyms 34, 56; 94 incertae_sedis 
#what are sibling_higher?

#get tree from OTOL
spectra_tree <- tol_induced_subtree(ott_ids = species_names_resolved$ott_id[which(!is.na(species_names_resolved$ott_id))])

#plot phylogeny
plot(spectra_tree, cex=0.5)

write.tree(spectra_tree, "./data/tol_tree.tre")

