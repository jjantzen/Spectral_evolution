#Rescaling phylogeny to height of 1 for interpretation of OU model

library(phytools)

#read trees
new_trees <- readRDS("./data/tidy/new_trees_matched_spectra.rds")

#rescale trees

#total height of 1
scale <- 1

#for each of the trees
for (i in 1:length(new_trees)){
  new_trees[[i]]$edge.length <- new_trees[[i]]$edge.length/max(nodeHeights(new_trees[[i]])[,2])*scale
}

#checking results
tree1 <- new_trees[[10]]
max(nodeHeights(tree1))
plot(tree1)

#save output
saveRDS(new_trees, "./data/tidy/rescaled_trees_matched_spectra.rds")
