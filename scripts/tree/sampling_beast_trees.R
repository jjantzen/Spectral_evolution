#randomly sample 100 trees from BEAST posterior

#load libraries
library(ape)
library(phytools)

#function
sample.trees<-function(trees, burnin, final.number, format){           
  
  #NUMBER OF TREES IN ORIGINAL FILE
  original.number<-as.numeric(length(trees))                           
  
  #THIS CREATES THE POST BURNIN PORTION
  post.burnin.trees<-trees[(burnin*original.number):original.number]   
  
  #THIS DOWNSAMPLES THE COLLECTION OF TREES
  final.trees<-sample(post.burnin.trees, final.number)                 
  
  #THIS SAVES THEM AS NEWICK FORMAT
  if(format=="new"){write.tree(final.trees, file="trees.nwk")}         
  
  #THIS SAVES THEM AS NEXUS FORMAT
  if(format=="nex"){write.nexus(final.trees, file="trees.nex")}
  
  return(final.trees)
}

#read trees
trees <- read.nexus("../Analysis non R/Beast/Redo_beast_partitioned_models/supermatrix_trimal80_3_cluster1546-out.trees")
str(trees)
#sample from trees

#also saves as trees.nwk in home folder
sampled_trees <- sample.trees(trees, 0.3, 100, "new")

#save as rds
saveRDS(sampled_trees, "./data/for_analysis/beast_sampled_100_trees_final.rds")
