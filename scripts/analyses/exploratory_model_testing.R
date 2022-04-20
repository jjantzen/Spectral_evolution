#exploratory model - most complex version?

library(mvMORPH)
library(phytools)
library(dplyr)

#read data
combo_data <- readRDS("./data/for_analysis/reduced_trait_data.rds")

#read tree
new_trees <- readRDS("./data/tidy/new_trees_matched_spectra.rds")

#prune trees
tree <- new_trees[[1]]

#tips to drop 
exclude <- new_trees[[1]]$tip.label[-which(new_trees[[1]]$tip.label %in%combo_data$Species)]

#drop tips
tree <- drop.tip(tree, exclude)

#read spectra
spec_data <- readRDS("./data/tidy/spectra_not_reordered_to_tree.rds")

#prune spec to data
spectra <- spec_data[which(rownames(spec_data) %in% combo_data$Species),]
str(spectra)

#make data list for model
data_spectra <- list(spectra=spectra, woody = combo_data$Woody, shade = combo_data$Shade, drought = combo_data$Drought_bin, leaf_pers = combo_data$leaf_persistence,
                     fine = combo_data$Fine_soil, coarse = combo_data$Coarse_soil)

colnames(combo_data)

str(data_spectra)

#prepare set of models to run
models <- c("BM", "OU", "EB", "lambda")

#additive model
model_fit_add_ou <- mvgls(spectra ~ woody + leaf_pers + shade + drought + coarse + fine, data=data_spectra, tree=tree, model=models[2], error = TRUE) 

#manova
ou_add_manova <- manova.gls(model_fit_add_ou, test = "Wilks", type = "II", nperm = 100, nbcores = 4L) #need to request 4 cores from cluster

#read in from cluster
ou_add_manova <- readRDS("./analysis/testing_models/testing_exploratory/ou_add_manova.rds")

model_fit_add_ou <- readRDS("./analysis/testing_models/testing_exploratory/model_fit_add_ou.rds")

ou_add_manova #nothing significant?

model_fit_add_ou$param

ou_add_manova_2 <- readRDS("./analysis/testing_models/testing_exploratory/ou_add_manova_2.rds")
ou_add_manova_2$stat
