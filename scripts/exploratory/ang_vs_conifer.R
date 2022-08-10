#Analysis for Angiosperms vs conifers across all wavelengths
library(ape)
library(phytools)
library(spectrolab)
library(stringr)
library(dplyr)
library(mvMORPH)

#1 get tree (single tree for now)
tree <- read.tree("./data/trees/trimal_80_concat.treefile") #read.beast for future
tree$tip.label
#root tree - first drop fern
tree <- drop.tip(tree, 117)
rooted_tree <- midpoint.root(tree) #, 99, resolve.root = TRUE) #, position=0.5*tree$edge.length[which(tree$edge[,2]==117)])
plot.phylo(rooted_tree2, show.node.label= FALSE)
rooted_tree2 <- reroot(tree, 95, resolve.root = TRUE, position=0.5*tree$edge.length[which(tree$edge[,2]==95)]) 

#2 get binary lineage variable (A vs C) and spectral data
ang_conifer <- read.csv("./data/tidy/ang_vs_conifer.csv", stringsAsFactors = FALSE, header = TRUE)

spectra <- readRDS("./data/tidy/trimmed_spectra_by_meta.rds")

#reformat names
spectra_df <- as.data.frame(spectra, row.names = spectra$names, metadata = TRUE)#, colnames = bands)

CABO_names <- word(spectra_df$species, 1,2, sep=" ")

spectra_df <- cbind(species_names = CABO_names, spectra_df)

#get mean of spectra
spectra_means <- spectra_df %>% 
  dplyr::select(-c(sample_name, sample_id, species, project)) %>% 
  group_by(species_names) %>% 
  summarize_all(mean, na.rm = TRUE)

#check for spectral taxa in tree and vice versa
rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

tree_taxa <- rooted_tree$tip.label

spectra_taxa <- as.character(unique(spectra_means$species_names))

# tree_taxa[-which(tree_taxa %in% spectra_taxa)]
# spectra_taxa[-which(spectra_taxa %in% tree$tip.label)]

rooted_tree$tip.label[which(rooted_tree$tip.label == "Eriophorum vaginatum subsp. spissum")] <- "Eriophorum vaginatum"
rooted_tree$tip.label[which(rooted_tree$tip.label == "Acer saccharum subsp. nigrum")] <- "Acer nigrum"
rooted_tree$tip.label[which(rooted_tree$tip.label == "Alnus incana subsp. rugosa")] <- "Alnus incana"

#drop tips and data
pruned_tree <- drop.tip(rooted_tree, rooted_tree$tip.label[-which(rooted_tree$tip.label %in% spectra_taxa)])
# pruned_tree$tip.label
#pruned_tree <- reroot(pruned_tree, 98, position=0.5*tree$edge.length[which(tree$edge[,2]==117)])
plot.phylo(pruned_tree, show.node.label= TRUE)

#make names character not factor
spectra_means$species_names <- as.character(spectra_means$species_names)

spectra_means <- spectra_means[which(spectra_means$species_names %in% pruned_tree$tip.label),]

#spectra_taxa2 <- as.character(unique(spectra_means$species_names))

#check for names
unique(spectra_means$species_names)
pruned_tree$tip.label

#3 format data
spectra_means_df <- data.frame(spectra_means, stringsAsFactors = FALSE)

row.names(spectra_means_df) <- as.character(spectra_means_df$species_names)
spectra_means_df <- spectra_means_df[,-1]
colnames(spectra_means_df) <- gsub("X", "", colnames(spectra_means_df))

colnames(spectra_means_df)

str(spectra_means_df)

# row.names(spectra_means_df)[-which(row.names(spectra_means_df) %in% pruned_tree$tip.label)]
# pruned_tree$tip.label[-which(pruned_tree$tip.label %in% row.names(spectra_means_df))]
# pruned_tree$tip.label[which(duplicated(pruned_tree$tip.label) == TRUE)]
ang_conifer$species
ang_conifer$species <- word(ang_conifer$species, 1,2, sep=" ")
#nrow(ang_conifer)
ang_conifer <-  ang_conifer[which(ang_conifer$species %in% pruned_tree$tip.label),]
#nrow(ang_conifer)
ang_conifer <- ang_conifer[-which(duplicated(ang_conifer) == TRUE),]

ang_conifer$angiosperm <- gsub("0", "Con", ang_conifer$angiosperm)
ang_conifer$angiosperm <- gsub("1", "Ang", ang_conifer$angiosperm)

# pruned_tree$tip.label[-which(pruned_tree$tip.label %in% ang_conifer$ï..species)]
lineage <- as.matrix(ang_conifer[,3])
rownames(lineage) <- ang_conifer[,1]

spec_list <- as.matrix(spectra_means_df)

dimnames(spec_list)[[1]]

#str(spec_list)
#str(lineage)
#str(pruned_tree)

#reorder lineage to match spectra
reorder_lin <- match(dimnames(spec_list)[[1]], dimnames(lineage)[[1]])
lineage_reordered <- as.matrix(lineage[reorder_lin])
rownames(lineage_reordered) <- dimnames(spec_list)[[1]]
str(lineage_reordered)

#combine data into one dataframe
data_spectra <- list(trait=spec_list, lineage=as.factor(lineage_reordered)) #matrix(spectra_means_df)

str(lineage_reordered)

write.csv(data_spectra, "./data/intermediate_objects/spectral_dataframe_ang_conifer.csv")
saveRDS(pruned_tree, "./data/intermediate_objects/pruned_ang_con_tree.rds")
saveRDS(data_spectra, "./data/intermediate_objects/data_spectra_ang_conifer.rds")

# data_spectra[which(is.na(data_spectra$trait))]
# missing <- data_spectra$trait[!complete.cases(data_spectra$trait),]
# data_spectra$trait[!is.finite(data_spectra$trait)]
# str(missing)
# missing <- spectra_means_df[!complete.cases(spectra_means_df),]
# 
# nans <- spectra_means_df[!is.nan(spectra_means_df),]
# 
# list_nans <- c()
# for (i in 1:ncol(spectra_means_df)){
#   nan <- is.nan(spectra_means_df[,i])
#  # nan <- paste0(nan, i)
#   list_nans <- c(list_nans, nan)
# }
# 
# TRUE %in% list_nans
# 
# 
# list_inf <- c()
# for (i in 1:ncol(spectra_means_df)){
#   infinite <- spectra_means_df[,i][!is.finite(spectra_means_df[,i])]
#   infinite <- paste0(infinite, i)
#   list_inf <- c(list_inf, infinite)
# }
# 
# list_zero <- c()
# for (i in 1:ncol(spectra_means_df)){
#   zero <- all(spectra_means_df[i,] == 0)
#   #zero <- paste0(zero, i)
#   list_zero <- c(list_zero, zero)
# }
# i <- 1
# TRUE %in% list_zero
# 
# length(list_zero)
# 
# spec_list[!is.finite(spec_list)]
# 
# spectra_means_df[,1]
# 
# class(spectra_means_df$`400`)
# 
# unique(spectra_means_df$`400`)
# 
# missing
# str(data$mandible)
# str(data_spectra$trait)
# str(data_spectra)
# str(pruned_tree)
# 


#4 model evolution (4 models)
fit_1_ac <- mvgls(trait ~ lineage, data=data_spectra, tree=pruned_tree, model="BM") 
saveRDS(fit_1_ac, "./data/intermediate_objects/fit_1_ac.rds")

#fit_1b_ac <- mvgls(trait ~ lineage, data=data_spectra, tree=pruned_tree, model="BM") #was going to try approx but only works for intercept
fit_2_ac <- mvgls(trait ~ lineage, data=data_spectra, tree=pruned_tree, model="OU") 
saveRDS(fit_2_ac, "./data/intermediate_objects/fit_2_ac.rds")

fit_3_ac <- mvgls(trait ~ lineage, data=data_spectra, tree=pruned_tree, model="EB") 
saveRDS(fit_3_ac, "./data/intermediate_objects/fit_3_ac.rds")

fit_4_ac <- mvgls(trait ~ lineage, data=data_spectra, tree=pruned_tree, model="lambda") 
saveRDS(fit_4_ac, "./data/intermediate_objects/fit_4_ac.rds")

#fit_1_int <- mvgls(trait ~ 1, data=data_spectra, tree=pruned_tree, model="BM") 

head(sort(unique(data_spectra$lineage)))

#5 model test
GIC(fit_1_ac) #
GIC(fit_2_ac) #
GIC(fit_3_ac) #
GIC(fit_4_ac) #

gics <- list(GIC(fit_1_ac)$GIC,GIC(fit_2_ac)$GIC,GIC(fit_3_ac)$GIC, GIC(fit_4_ac)$GIC)

capture.output(gics, file = "./gic_ang_conifer.txt")

# #do manovas
# aov_1_ac <- manova.gls(fit_1_ac, nperm=999, test="Wilks", verbose=TRUE, nbcores = 4)
# aov_2_ac <- manova.gls(fit_2_ac, nperm=999, test="Wilks", verbose=TRUE, nbcores = 4)
# aov_3_ac <- manova.gls(fit_3_ac, nperm=999, test="Wilks", verbose=TRUE, nbcores = 4)
# aov_4_ac <- manova.gls(fit_4_ac, nperm=999, test="Wilks", verbose=TRUE, nbcores = 4)

#6 test impact of lineage on spectra

#7 PCA

#8 univariate models
