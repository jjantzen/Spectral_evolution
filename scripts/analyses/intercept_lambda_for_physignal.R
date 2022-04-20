#Phylosignal of spectra

library(geomorph)
library(mvMORPH)
library(phytools)

#read objects
spectra_matrix <- readRDS("./spectra_not_reordered_to_tree.rds")
data_spectra <- list(trait=spectra_matrix)

new_trees <- readRDS("./new_trees_matched_spectra.rds")
tree1 <- new_trees[[1]]

#get lambda for entire model
fit_lambda_intercept <- mvgls(trait ~ 1, data=data_spectra, tree=tree1, model="lambda") 
saveRDS(fit_lambda_intercept, "./fit_lambda_intercept.rds")

