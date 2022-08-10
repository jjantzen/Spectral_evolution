#### Setting up session ####
#from online session with Florence

# install.packages(c("devtools", "tidyverse", "ape", "vegan", "ade4", "phytools", "mvMORPH", "geiger", "geomorph", 
#                   "hsdar", "corrplot", "reshape", "agricolae"))

library(devtools)
install_version("caret", version = "6.0-84", repos = "http://cran.us.r-project.org")
install_version("pls", version = "2.7-2", repos = "http://cran.us.r-project.org")

sourceurl<-"https://cran.r-project.org/src/contrib/Archive/prospectr/prospectr_0.2.0.tar.gz"
install.packages(sourceurl, repos = NULL, type = "source")

install_version("spectrolab", version = "0.0.13", repos = "http://cran.us.r-project.org") #skip prospectr update
#install_github("meireles/spectrolab")
rm(sourceurl)

#### Load libraries ####
library(mvMORPH)
library(spectrolab)
library(tidyverse)

#read data and tree
spectra <- readRDS("./data/tidy/spectra_spnames.rds")
spectra_FL <- readRDS("./data/Florence/allclean_all.rds")
tree <- read.tree("./data/tol_brlen_grafen1_lam0.tre")
taxonomy <- read.csv("./data/tidy/resolved_taxonomy.csv")

#resolve names for spectra and tree by making table for conversions

#make spectra dataframe
spec_rename <- as.data.frame(spectra)

#create species codes
names_conversion <- data.frame(tree_names = tree$tip.label)

names_conversion <- names_conversion %>% 
  mutate(genus_code = toupper(substr(tree_names, start = 1, stop = 2)), species_code = toupper(substr((sapply(strsplit(as.character(tree_names), "_", fixed = TRUE), "[", 2)), start = 1, stop = 2))) %>% 
  mutate(code = paste0(genus_code, species_code))

length(unique(names_conversion$tree_names))

#rename tip labels
tree$tip.label <- names_conversion$code

plot(tree, cex = 0.4)

#add code column to spectral

splitting_species_names <- function(species_names_dataframe){
  sample_meta_names <- species_names_dataframe %>% 
    separate(species, c("genus", "specific_epithet", "extra0", "extra1", "extra2", "extra3", "extra4", "extra5", "extra6", "extra7", "extra8", "extra9", "extra10"), sep = " ", fill = "right", remove = FALSE)
  
  #deal with subspecific names
  sample_meta_names$subspecific <- NA
  sample_meta_names$subspecific[which(sample_meta_names$extra5 == "australis")] <- paste0(sample_meta_names$extra4[which(sample_meta_names$extra5 == "australis")], "_", sample_meta_names$extra5[which(sample_meta_names$extra5 == "australis")])
  
  sample_meta_names$subspecific[which(sample_meta_names$extra2 == "angustifolia")] <- paste0(sample_meta_names$extra1[which(sample_meta_names$extra2 == "angustifolia")], "_", sample_meta_names$extra2[which(sample_meta_names$extra2 == "angustifolia")])
  sample_meta_names$subspecific[which(sample_meta_names$extra2 == "vaginatum")] <- paste0(sample_meta_names$extra1[which(sample_meta_names$extra2 == "vaginatum")], "_", sample_meta_names$extra2[which(sample_meta_names$extra2 == "vaginatum")]) 
  
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "spissum")] <- paste0(sample_meta_names$extra0[which(sample_meta_names$extra1 == "spissum")], "_", sample_meta_names$extra1[which(sample_meta_names$extra1 == "spissum")])
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "strigosus")] <- paste0(sample_meta_names$extra0[which(sample_meta_names$extra1 == "strigosus")], "_", sample_meta_names$extra1[which(sample_meta_names$extra1 == "strigosus")]) 
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "rugosa")] <- paste0(sample_meta_names$extra0[which(sample_meta_names$extra1 == "rugosa")], "_", sample_meta_names$extra1[which(sample_meta_names$extra1 == "rugosa")])
  return(sample_meta_names)
}

spec_rename_split <- splitting_species_names(spec_rename)

spec_rename_split <- spec_rename_split[,-c(6:16)]

species_names_spectral <- spec_rename_split %>% 
  mutate(genus_code = toupper(substr(genus, start = 1, stop = 2)), species_code = toupper(substr(specific_epithet, start = 1, stop = 2))) %>% 
  mutate(code = paste0(genus_code, species_code)) %>% 
  unite(species_name, genus, specific_epithet, subspecific, sep = " ", na.rm = TRUE)

species_names_spectral$code

spec_rename_split[,1:5]

#deal with subspecies
acers <- species_names_spectral[which(species_names_spectral$code == "ACSA")] 

species_names_spectral$code[which(species_names_spectral$species_name == "Acer saccharum")] <- "ACSA1"

#trim spectral data to match tree
spectral_trimmed <- species_names_spectral[which(species_names_spectral$code %in% tree$tip.label),]
length(unique(spectral_trimmed$code))
length(unique(tree$tip.label))
length(tree$tip.label)

tree$tip.label[duplicated(tree$tip.label)]

tree$tip.label[which(tree$tip.label == "ACSA")] <- c("ACSA1", "ACSA2")
tree$tip.label[which(tree$tip.label == "VILA")] <- c("VILA1", "VILA2")
tree$tip.label[which(tree$tip.label == "CACA")] <- c("CACA1", "CACA2")
tree$tip.label[which(tree$tip.label == "ERVA")] <- c("ERVA1", "ERVA2")



missing <- tree$tip.label[-which(unique(tree$tip.label) %in% unique(spectral_trimmed$code))]

names_conversion$tree_names[which(names_conversion$code %in% missing)]
unique(spectral_trimmed$code[which(spectral_trimmed$code %in% missing)])

sort(unique(spectral_trimmed$code))

spectral_trimmed[which(spectral_trimmed$code == "POMU"),1:5]
spectral_trimmed[which(spectral_trimmed$code == "OSRE"),1:5]
spectral_trimmed[which(spectral_trimmed$code == "ONSE"),1:5]
spectral_trimmed[which(spectral_trimmed$code == "TUFA"),1:5]
spectral_trimmed[which(spectral_trimmed$code == "GIBI"),1:5]
spectral_trimmed[which(spectral_trimmed$code == "CLCL"),1:5]




colnames(spec_rename_split)
spec_rename_split$
str(spec_rename)
unique(spec_rename$species)

#remove 6:16
#fix two with subsp names

#rename tips of tree and dataset to match


taxonomy$unique_name %in% tolower(spec_rename$species)

#compare to florence data
str(spectra)
str(spectra_FL)



tree$tip.label
spec_rename$species


#trim dataset by taxa in tree
spec_t <- specs_s[which(specs_s$code %in% tree$tip.label),]

#first filter by property to get just reflectance not transmittance
spec_keep <- spec_t[which(spec_t$property == "reflectance"),]
spec_keep
nrow(spec_keep)
str(spec_keep)

str(spec_keep$vascan_taxon)

spec_keep_df <- as.data.frame(spec_keep)

max(table(spec_keep$vascan_taxon))

#get spectral data
spec_data <- as.matrix(spec_keep[,16:2016])

nrow(spec_data)

#get trait data
genus <- spec_keep$genus
length(genus)

nrow(spec_data)
str(genus)

#combine into list of dataframes
data_list <- list(spectra=spec_data, genus = as.factor(genus))

str(spec_data)
str(trait)
str(data_list$spectra)
str(data)
str(data_list)

# Fit the multivariate linear model
site_spectra_model <- mvgls(spectra ~ genus, data=data_list, tree=tree, model="lambda")


#have an issue with sampling - individual vs species level (tips (sp) don't match the data (ind))


# Print the model fit
print(site_spectra_model)
summary(site_spectra_model)
# The overall MANOVA test:
aov <- manova.gls(site_spectra_model,
                  nperm=999, test="Wilks", verbose=TRUE)

# display the test statistic
plot(aov)

# Fit with parallel computing on 4 cores for 9999 permutations
aov2 <- manova.gls(fit,
                   nperm=9999, nbcores=4L, test="Wilks", verbose=TRUE)

# display the test statistic
plot(aov2)

#2 variables
# Fit the MANCOVA model
fit2 <- mvgls(trait ~ size + habitat, data=data, tree=tree, model="lambda")

#The overall MANOVA test:
aov3 <- manova.gls(fit2, nperm=999, test="Wilks", verbose=TRUE)

#The overall (type II) MANOVA test:
aov4 <- manova.gls(fit2, nperm=999, test="Wilks", type="II", verbose=TRUE)

#Fit the MANCOVA model with interaction term
fit3 <- mvgls(trait ~ size + habitat + size*habitat, data=data, tree=tree, model="lambda")

# The overall (type II) MANOVA test:
aov5 <- manova.gls(fit3, nperm=999, test="Wilks", type="II", verbose=TRUE)

#add a random effect for species (intraspecific) 
#do mean for species
#do polytomy of length 0 for individuals per species
#make species identity a random effect with phylogeny while looking at relationship between traits of interest
#mcmcglmm - bayesian package



