#Phylosignal of spectra

library(geomorph)
library(mvMORPH)
library(phytools)

#calculate Kmulti using physignal function - for procrustes shape variables

#read spectra
spectra_matrix <- readRDS("./data/for_analysis/spectra_not_reordered_to_tree.rds")

str(spectra_matrix)

#data_spectra <- list(trait=spectra_matrix)

#read trees
new_trees <- readRDS("./data/for_analysis/final_trees_matched_spectra.rds")

#rescaled_trees <- readRDS("./data/tidy/rescaled_trees_matched_spectra.rds")

#A is matrix or 3D array
#phy is tree of class phylo
#physignal(A, phy, iter = 9, seed = NULL, print.progress = TRUE) #999

#calculates p value with significance

#full spectra


df_output <- data.frame(matrix(nrow=100, ncol = 3))
colnames(df_output)[1] <- "iteration"
class(df_output[,1]) <- "numeric"
colnames(df_output)[2] <- "K"
class(df_output[,2]) <- "numeric"
colnames(df_output)[3] <- "pvalue"
class(df_output[,3]) <- "numeric"


for (i in 1:length(new_trees)){
  sig <- physignal(spectra_matrix, new_trees[[i]], iter = 999, seed = NULL, print.progress = TRUE) 
  df_output$iteration[i] <- i
  df_output$K[i] <- sig$phy.signal
  df_output$pvalue[i] <- sig$pvalue
}

#extract K and pvalue from output lists
df_output2 <- df_output %>% add_row(iteration = NA, K = mean(df_output$K[c(1:100)]), pvalue =  mean(df_output$pvalue[1:100]))

df_output2 <- df_output2 %>% add_row(iteration = NA, K = max(df_output2$K[c(1:100)]), pvalue =  max(df_output2$pvalue[1:100]))

# df_output$iteration[101] <- "mean"
# df_output$K[101] <-  mean(df_output$K[c(1:100)]) #0.1130031
# df_output$pvalue[101] <-  mean(df_output$pvalue[1:100]) #0.00144
# df_output$iteration[102] <- "max"
# df_output$K[102] <-  max(df_output$K[c(1:100)]) #0.1754643
# df_output$pvalue[102] <-  max(df_output$pvalue[1:100]) #0.03
# 
# 
# #calculate mean
# #old values
# mean(df_output$K) #0.1754853 
# mean(df_output$pvalue) #0.00158
# max(df_output$pvalue) #max 0.029 so all significant

#save output
write.csv(df_output2, "./analysis/phylosig/Kmulti_100trees_final.csv")


#prune phylogeny to clade and family and calculate phylosignal for each across all 100 trees
tree_names <- new_trees[[1]]$tip.label

#import taxonomic data
data_set <- readRDS("./data/for_analysis/final_data.rds")

colnames(data_set)

taxonomy <- data_set[,c(1,6,7,8,9)] #fill in right columns - species names, order and subclass and superorder and family
taxonomy

example_tree <- new_trees[[99]]

findMRCA(example_tree, tips=taxonomy$Species[which(taxonomy$Order == "Vitales")], type="node")

#~subclass
conifer_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Subclass == "Pinidae")], type="node") #specify nodes c(165)
conifer_clade <- c()
for (i in 1:length(new_trees)){
  conifer_clade[[i]] <- extract.clade(new_trees[[i]], node=conifer_nodes[[i]]) #specify nodes c(165)
}
#conifer_clade <- lapply(new_trees,extract.clade,node=c(103)) #specify nodes
ang_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Subclass == "Magnoliidae")], type="node") #specify nodes c(165)
ang_clade <- c()
for (i in 1:length(new_trees)){
  ang_clade[[i]] <- extract.clade(new_trees[[i]], node=ang_nodes[[i]]) #specify nodes c(165)
}
#ang_clade <- lapply(new_trees,extract.clade,node=c(113)) #specify nodes

#superorder (ang only)
rosanae_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Superorder == "Rosanae")], type="node") #specify nodes c(165)
rosanae_clade <- c()
for (i in 1:length(new_trees)){
  rosanae_clade[[i]] <- extract.clade(new_trees[[i]], node=rosanae_nodes[[i]]) #specify nodes c(165)
}
#rosanae_clade <- lapply(new_trees,extract.clade,node=c(116)) #specify nodes
asteranae_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Superorder == "Asteranae")], type="node") #specify nodes c(165)
asteranae_clade <- c()
for (i in 1:length(new_trees)){
  asteranae_clade[[i]] <- extract.clade(new_trees[[i]], node=asteranae_nodes[[i]]) #specify nodes c(165)
}
#asteranae_clade <- lapply(new_trees,extract.clade,node=c(171)) #specify nodes
#ranunculanae_clade <- lapply(new_trees,extract.clade,node=c()) #specify nodes - no mrca
#caryophyllanae_clade <- lapply(new_trees,extract.clade,node=c()) #specify nodes - no mrca
lilianae_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Superorder == "Lilianae")], type="node") #specify nodes c(165)
lilianae_clade <- c()
for (i in 1:length(new_trees)){
  lilianae_clade[[i]] <- extract.clade(new_trees[[i]], node=lilianae_nodes[[i]]) #specify nodes c(165)
}
#lilianae_clade <- lapply(new_trees,extract.clade,node=c(187)) #specify nodes

#order
pinales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Pinales")], type="node") #specify nodes c(165)
pinales_clade <- c()
for (i in 1:length(new_trees)){
  pinales_clade[[i]] <- extract.clade(new_trees[[i]], node=pinales_nodes[[i]]) #specify nodes c(165)
}
#pinales_clade <- lapply(new_trees,extract.clade,node=c(104)) #specify nodes
sapindales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Sapindales")], type="node") #specify nodes c(165)
sapindales_clade <- c()
for (i in 1:length(new_trees)){
  sapindales_clade[[i]] <- extract.clade(new_trees[[i]], node=sapindales_nodes[[i]]) #specify nodes c(165)
}
#sapindales_clade <- lapply(new_trees,extract.clade,node=c(119)) #specify nodes
fagales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Fagales")], type="node") #specify nodes c(165)
fagales_clade <- c()
for (i in 1:length(new_trees)){
  fagales_clade[[i]] <- extract.clade(new_trees[[i]], node=fagales_nodes[[i]]) #specify nodes c(165)
}
#fagales_clade <- lapply(new_trees,extract.clade,node=c(131)) #specify nodes
gentianales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Gentianales")], type="node") #specify nodes c(165)
gentianales_clade <- c()
for (i in 1:length(new_trees)){
  gentianales_clade[[i]] <- extract.clade(new_trees[[i]], node=gentianales_nodes[[i]]) #specify nodes c(165)
}
#gentianales_clade <- lapply(new_trees,extract.clade,node=c(175)) #specify nodes
#ranunculales_clade <- lapply(new_trees,extract.clade,node=c()) #specify nodes - no mrca
poales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Poales")], type="node") #specify nodes c(165)
poales_clade <- c()
for (i in 1:length(new_trees)){
  poales_clade[[i]] <- extract.clade(new_trees[[i]], node=poales_nodes[[i]]) #specify nodes c(165)
}
#poales_clade <- lapply(new_trees,extract.clade,node=c(188)) #specify nodes
asparagales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Asparagales")], type="node") #specify nodes c(165)
asparagales_clade <- c()
for (i in 1:length(new_trees)){
  asparagales_clade[[i]] <- extract.clade(new_trees[[i]], node=asparagales_nodes[[i]]) #specify nodes c(165)
}
#asparagales_clade <- lapply(new_trees,extract.clade,node=c(199)) #specify nodes
rosales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Rosales")], type="node") #specify nodes c(165)
rosales_clade <- c()
for (i in 1:length(new_trees)){
  rosales_clade[[i]] <- extract.clade(new_trees[[i]], node=rosales_nodes[[i]]) #specify nodes c(165)
}
#rosales_clade <- lapply(new_trees,extract.clade,node=c(147)) #specify nodes
ericales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Ericales")], type="node") #specify nodes c(165)
ericales_clade <- c()
for (i in 1:length(new_trees)){
  ericales_clade[[i]] <- extract.clade(new_trees[[i]], node=ericales_nodes[[i]]) #specify nodes c(165)
}
#ericales_clade <- lapply(new_trees,extract.clade,node=c(185)) #specify nodes
asterales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Asterales")], type="node") #specify nodes c(165)
asterales_clade <- c()
for (i in 1:length(new_trees)){
  asterales_clade[[i]] <- extract.clade(new_trees[[i]], node=asterales_nodes[[i]]) #specify nodes c(165)
}
#asterales_clade <- lapply(new_trees,extract.clade,node=c(180)) #specify nodes
#caryophyllales_clade <- lapply(new_trees,extract.clade,node=c()) #specify nodes - no mrca
fabales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Fabales")], type="node") #specify nodes c(165)
fabales_clade <- c()
for (i in 1:length(new_trees)){
  fabales_clade[[i]] <- extract.clade(new_trees[[i]], node=fabales_nodes[[i]]) #specify nodes c(165)
}
#fabales_clade <- lapply(new_trees,extract.clade,node=c(162)) #specify nodes
#cornales_clade <- lapply(new_trees,extract.clade,node=c()) #specify nodes - no mrca
lamiales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Lamiales")], type="node") #specify nodes c(165)
lamiales_clade <- c()
for (i in 1:length(new_trees)){
  lamiales_clade[[i]] <- extract.clade(new_trees[[i]], node=lamiales_nodes[[i]]) #specify nodes c(165)
}
#lamiales_clade <- lapply(new_trees,extract.clade,node=c(176)) #specify nodes
apiales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Apiales")], type="node") #specify nodes c(165)
apiales_clade <- c()
for (i in 1:length(new_trees)){
  apiales_clade[[i]] <- extract.clade(new_trees[[i]], node=apiales_nodes[[i]]) #specify nodes c(165)
}
#apiales_clade <- lapply(new_trees,extract.clade,node=c(183)) #specify nodes
#myrtales_clade <- lapply(new_trees,extract.clade,node=c()) #specify nodes - no mrca
dipsacales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Dipsacales")], type="node") #specify nodes c(165)
dipsacales_clade <- c()
for (i in 1:length(new_trees)){
  dipsacales_clade[[i]] <- extract.clade(new_trees[[i]], node=dipsacales_nodes[[i]]) #specify nodes c(165)
}
#dipsacales_clade <- lapply(new_trees,extract.clade,node=c(184)) #specify nodes
#polypodiales_clade <- lapply(new_trees,extract.clade,node=c()) #specify nodes - no mrca
malpighiales_nodes <- lapply(new_trees,findMRCA,tips=taxonomy$Species[which(taxonomy$Order == "Malpighiales")], type="node") #specify nodes c(165)
malpighiales_clade <- c()
for (i in 1:length(new_trees)){
  malpighiales_clade[[i]] <- extract.clade(new_trees[[i]], node=malpighiales_nodes[[i]]) #specify nodes c(165)
}

#cupressales_clade <- lapply(new_trees,extract.clade,node=c()) #specify nodes - no mrca
#malvales_clade <- lapply(new_trees,extract.clade,node=c()) #specify nodes - no mrca
#vitales_clade <- lapply(new_trees,extract.clade,node=c()) #specify nodes - no mrca

#not going to calculate by family because too many (unless seems interesting for bigger orders)

#extract.clade(phy, node, root.edge = 0, collapse.singles = TRUE, interactive = FALSE)

#run function for each set of trees

#need to fix being able to call the specific name to match
#error when trying to run function but works on its own

physig_iterations <- function(spectra_matrix, new_trees, name){#, taxonomy){ #exclude taxonomy for spectral regions
  #exclude trimming for spectral regions
  #spectra_matrix_sm <- spectra_matrix[which(dimnames(spectra_matrix)[[1]] %in% taxonomy$Species[which(taxonomy$Subclass == name)]),] #modify this according to clade level
  #str(spectra_matrix_sm)
  str(new_trees[[2]])
  df_output <- data.frame(matrix(nrow=102, ncol = 3))
  colnames(df_output)[1] <- "iteration"
  class(df_output[,1]) <- "numeric"
  colnames(df_output)[2] <- "K"
  class(df_output[,2]) <- "numeric"
  colnames(df_output)[3] <- "pvalue"
  class(df_output[,3]) <- "numeric"
  
  
  #do physignal calc
  for (i in 1:length(new_trees)){
    sig <- physignal(spectra_matrix, new_trees[[i]], iter = 999, seed = NULL, print.progress = TRUE) 
    df_output$iteration[i] <- i
    df_output$K[i] <- sig$phy.signal
    df_output$pvalue[i] <- sig$pvalue
  }
  
  df_output$iteration[101] <- NA
  df_output$K[101] <-  mean(df_output$K[c(1:100)])
  df_output$pvalue[101] <-  mean(df_output$pvalue[1:100])
  df_output$iteration[102] <- NA
  df_output$K[102] <-  max(df_output$K[c(1:100)])
  df_output$pvalue[102] <-  max(df_output$pvalue[1:100])
  
  saveRDS(df_output, paste0("./analysis/phylosig/Kmulti_100trees_" , name, "_clade.rds"))
  write.csv(df_output, paste0("./analysis/phylosig/Kmulti_100trees_" , name, "_clade.csv"))
}


#run for different clades
#(spectra_matrix, new_trees, name, level, taxonomy)

conifer_K <- physig_iterations(spectra_matrix, conifer_clade, "Pinidae", taxonomy)
ang_K <- physig_iterations(spectra_matrix, ang_clade, "Magnoliidae", taxonomy)

rosanae_K <- physig_iterations(spectra_matrix, rosanae_clade, "Rosanae", taxonomy)
asteranae_K <- physig_iterations(spectra_matrix, asteranae_clade, "Asteranae", taxonomy)
#ranunculanae_K <- physig_iterations(spectra_matrix, ranunculanae_clade, "ranunculanae")
#caryophyllanae_K <- physig_iterations(spectra_matrix, caryophyllanae_clade, "caryophyllanae")
lilianae_K <- physig_iterations(spectra_matrix, lilianae_clade, "Lilianae", taxonomy)

pinales_K <- physig_iterations(spectra_matrix, pinales_clade, "Pinales", taxonomy)
sapindales_K <- physig_iterations(spectra_matrix, sapindales_clade, "Sapindales", taxonomy)
fagales_K <- physig_iterations(spectra_matrix, fagales_clade, "Fagales", taxonomy)
#gentianales_K <- physig_iterations(spectra_matrix, gentianales_clade, "Gentianales", taxonomy) #too few taxa
#ranunculales_K <- physig_iterations(spectra_matrix, ranunculales_clade, "ranunculales")
poales_K <- physig_iterations(spectra_matrix, poales_clade, "Poales", taxonomy)
#asparagales_K <- physig_iterations(spectra_matrix, asparagales_clade, "Asparagales", taxonomy) #too few taxa
rosales_K <- physig_iterations(spectra_matrix, rosales_clade, "Rosales", taxonomy)
ericales_K <- physig_iterations(spectra_matrix, ericales_clade, "Ericales", taxonomy)
asterales_K <- physig_iterations(spectra_matrix, asterales_clade, "Asterales", taxonomy)
#caryophyllales_K <- physig_iterations(spectra_matrix, caryophyllales_clade, "caryophyllales")
fabales_K <- physig_iterations(spectra_matrix, fabales_clade, "Fabales", taxonomy)
#cornales_K <- physig_iterations(spectra_matrix, cornales_clade, "cornales")
lamiales_K <- physig_iterations(spectra_matrix, lamiales_clade, "Lamiales", taxonomy)
#apiales_K <- physig_iterations(spectra_matrix, apiales_clade, "Apiales", taxonomy) #too few taxa
#myrtales_K <- physig_iterations(spectra_matrix, myrtales_clade, "myrtales")
#dipsacales_K <- physig_iterations(spectra_matrix, dipsacales_clade, "Dipsacales", taxonomy) #too few taxa
#polypodiales_K <- physig_iterations(spectra_matrix, polypodiales_clade, "polypodiales")
malpighiales_K <- physig_iterations(spectra_matrix, malpighiales_clade, "Malpighiales", taxonomy)
#cupressales_K <- physig_iterations(spectra_matrix, cupressales_clade, "cupressales")
#malvales_K <- physig_iterations(spectra_matrix, malvales_clade, "malvales")
#vitales_K <- physig_iterations(spectra_matrix, vitales_clade, "vitales")



#do the same but for different wavelengths

#do visible, NIR and SWIR
str(spectra_matrix)

#visible 400-699 nm
#NIR 700-1399 nm
#SWIR 1400-2400 nm

#doublecheck wavelengths correspond to right columns
visible_K <- physig_iterations(spectra_matrix[,1:300], new_trees, "visible")
NIR_K <- physig_iterations(spectra_matrix[,301:1000], new_trees, "NIR")
SWIR_K <- physig_iterations(spectra_matrix[,1001:2001], new_trees, "SWIR")

#getting stats

mean(Kmulti_100trees_SWIR_clade$K[c(1:100)])
sd(Kmulti_100trees_SWIR_clade$K[c(1:100)])
mean(Kmulti_100trees_SWIR_clade$pvalue[c(1:100)])
sd(Kmulti_100trees_SWIR_clade$pvalue[c(1:100)])

nrow(Kmulti_100trees_SWIR_clade[c(1:100),][which((Kmulti_100trees_SWIR_clade$pvalue[c(1:100)] < 0.05) == TRUE),])

Kmulti_100trees_final <- read.csv("./analysis/phylosig/Kmulti_100trees_final.csv")

Kmulti_100trees_final_2 <- Kmulti_100trees_final[c(1:100),]
