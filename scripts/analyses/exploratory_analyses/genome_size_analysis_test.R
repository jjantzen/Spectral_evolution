#genome size model trial
library(phytools)
library(spectrolab)

#read trees
trees <- readRDS("./data/for_analysis/final_trees_matched_spectra.rds")
trees
consensus_tree <- read.nexus("./data/for_analysis/consensus_tree.nex")
# plot(consensus_tree)
consensus_tree$tip.label <- gsub("_", " ", consensus_tree$tip.label)
consensus_tree$tip.label[which(consensus_tree$tip.label == "'Alnus incana subsp. rugosa'")] <- "Alnus incana"
consensus_tree$tip.label[which(consensus_tree$tip.label == "'Acer saccharum subsp. nigrum'")] <- "Acer nigrum"

#saveRDS(consensus_tree, "./data/for_analysis/consensus_tree_full.rds")

#read spectral data
spectra <- readRDS("./data/for_analysis/spectra_not_reordered_to_tree.rds")
spectra

#read genome size
ploidy <- read.csv("./data/predictors/c-value.csv", stringsAsFactors = FALSE)
ploidy

#get small dataset
ploidy_sm <- ploidy[which(!is.na(ploidy$cvalue)),c(1:2)]

#choose one value for multi-state species - chose the smallest value
ploidy_sm$cvalue[which(ploidy_sm$species == "Poa pratensis")] <- 4.24
ploidy_sm$cvalue[which(ploidy_sm$species == "Rubus idaeus")] <- 0.29
ploidy_sm$cvalue[which(ploidy_sm$species == "Phragmites australis")] <- 2.0
ploidy_sm$cvalue[which(ploidy_sm$species == "Claytonia perfoliata")] <- 1.47

ploidy_sm$cvalue <- as.numeric(ploidy_sm$cvalue)

#prune tree to list of species in ploidy dataset
tree <- trees[[1]]

exclude <- tree$tip.label[-which(tree$tip.label %in% ploidy_sm$species)]


genome_trees <- lapply(trees,drop.tip,tip=exclude)

genome_consensus_tree <- drop.tip(consensus_tree, tip = exclude)
genome_consensus_tree$tip.label

#tree_sm <-  drop.tip(tree, exclude)

#trim spectra to match species in ploidy dataset
spectra_sm <- spectra[which(dimnames(spectra)[[1]] %in% ploidy_sm$species),]

dimnames(spectra_sm)[1]

#save genome size dataframes
#saveRDS(spectra_sm, "./data/for_analysis/genome_size_spectra.rds")
saveRDS(genome_trees, "./data/for_analysis/genome_size_trees.rds")
#saveRDS(ploidy_sm, "./data/for_analysis/genome_size_data.rds")
saveRDS(genome_consensus_tree, "./data/for_analysis/genome_consensus_tree.rds")


spectra_sm <- readRDS("./data/for_analysis/genome_size_spectra.rds")
ploidy_sm <- readRDS("./data/for_analysis/genome_size_data.rds")



#calculate phylogenetic signal in genome size
library(phytools)

#make named list
genome_size <- setNames(as.numeric(ploidy_sm$cvalue), ploidy_sm$species)

class(genome_size)

genome_size <- genome_size[!is.na(genome_size)]

df_output_genome <- data.frame(matrix(nrow=100, ncol = 3))
colnames(df_output_genome)[1] <- "iteration"
class(df_output_genome[,1]) <- "numeric"
colnames(df_output_genome)[2] <- "K"
class(df_output_genome[,2]) <- "numeric"
colnames(df_output_genome)[3] <- "Pvalue"
class(df_output_genome[,3]) <- "numeric"


for (i in 1:length(genome_trees)){
  out <- phylosig(genome_trees[[i]], genome_size, method="K", test=TRUE, nsim=999)
  df_output_genome$iteration[i] <- i
  df_output_genome[i,2] <- out$K
  df_output_genome[i,3] <- out$P
}

# gs_K <- phylosig(genome_trees[[1]], genome_size, method="K", test=TRUE, nsim=999) #0.411066 pvalue = 0.001001
# gs_lambda <- phylosig(tree_sm, genome_size, method="lambda", test=TRUE, nsim=999) #0.938108 pvalue = ~0

saveRDS(df_output_genome, "./analysis/trait_phylosig/genome_size_phylosig_K.rds")
saveRDS(df_output_genome_l, "./analysis/trait_phylosig/genome_size_phylosig_lambda.rds")


sd(df_output_genome$K)

nrow(df_output_genome[which((df_output_genome$Pvalue < 0.05) == TRUE),])


#plot gs on tree
dotTree(tree_sm,genome_size,length=10,ftype="i")

min(genome_size)

#doesn't work ylim error
# jpeg("./output/trait_plots/Genome_size_barplot.jpg", width = 10, height = 10, res = 400, units = "in")
# plotTree.wbars(tree_sm,ploidy_sm$cvalue) #args.barplot=list(beside=TRUE,xlim=c(0,30), height = )
# dev.off()

jpeg("./output/trait_plots/Genome_size_phenogram.jpg", width = 10, height = 10, res = 400, units = "in")
phenogram(genome_consensus_tree,genome_size,spread.labels=TRUE,spread.cost=c(10,0), fsize = 0.7)
dev.off()

#get list of conifers/fern
conifers <- c("Abies balsamea", "Larix laricina", "Picea abies", "Picea glauca", "Picea mariana", "Pinus banksiana", 
              "Pinus resinosa", "Pinus rigida", "Pinus strobus", "Polystichum munitum", "Thuja occidentalis", "Tsuga canadensis")

ang_tree <- drop.tip(genome_consensus_tree, conifers)
plot(ang_tree)

jpeg("./output/trait_plots/Genome_size_phenogram_ang_only.jpg", width = 10, height = 10, res = 400, units = "in")
phenogram(ang_tree,genome_size[-which(names(genome_size) %in% conifers)],spread.labels=TRUE,spread.cost=c(10,0), fsize = 0.7)
dev.off()

#phylosignal for angiosperms only 
gs_K_ang <- phylosig(ang_tree, genome_size[-which(names(genome_size) %in% conifers)], method="K", test=TRUE, nsim=999) #0.286828 pvalue = 0.00500501
gs_lambda_ang <- phylosig(ang_tree, genome_size[-which(names(genome_size) %in% conifers)], method="lambda", test=TRUE, nsim=999) #0.667931 pvalue = ~0

#conifer only 
non_conifers <- tree_sm$tip.label[-which(tree_sm$tip.label %in% conifers)]
con_tree <- drop.tip(tree_sm, non_conifers)
plot(con_tree)

jpeg("./output/trait_plots/Genome_size_phenogram_conifer_only.jpg", width = 10, height = 10, res = 400, units = "in")
phenogram(con_tree,genome_size[-which(names(genome_size) %in% non_conifers)],spread.labels=TRUE,spread.cost=c(10,0), fsize = 0.7)
dev.off()

#more plots
jpeg("./output/trait_plots/Genome_size_heatmap.jpg", width = 10, height = 12, res = 400, units = "in")
obj<-contMap(genome_consensus_tree,genome_size,plot=FALSE)
plot(obj,lwd=7)#,xlim=c(-0.2,3.6)
#errorbar.contMap(obj)
dev.off()


#models
library(mvMORPH)

#make data list
gs_spectra <- list(spectra=spectra_sm, gs = ploidy_sm$cvalue)
class(gs_spectra$gs) <- "numeric"

#run models
# gs_bm_model <- mvgls(spectra ~ gs, data = gs_spectra, tree = tree_sm, model = "BM", error = TRUE)
# gs_ou_model <- mvgls(spectra ~ gs, data = gs_spectra, tree = tree_sm, model = "OU", error = TRUE)
# gs_eb_model <- mvgls(spectra ~ gs, data = gs_spectra, tree = tree_sm, model = "EB", error = TRUE)
# gs_lambda_model <- mvgls(spectra ~ gs, data = gs_spectra, tree = tree_sm, model = "lambda", error = TRUE)

#create output dataframe
df_output <- data.frame(matrix(nrow=4, ncol = 3))
colnames(df_output) <- c("model", "GIC", "parameter")

class(df_output[,c(2:3)]) <- "numeric"
class(df_output[,1]) <- "character"

#run for each of 4 models of evolution
models <- c("BM", "OU", "EB", "lambda")

for (j in 1:length(models)){
  model <- mvgls(spectra ~ gs, data=gs_spectra, tree=tree_sm, model=models[j], error = TRUE)
  #assign output
  df_output$model[j] <- models[j]
  df_output$parameter[j] <- model$param
  df_output$GIC[j] <- GIC(model)$GIC
  #save model itself
  saveRDS(model, paste0("./gs_models/gs_model_", models[j], ".rds"))
}

saveRDS(df_output, "./gs_models/model_parameters_gs.rds")

#look at output
gics <- readRDS("./analysis/testing_models/testing_gs/model_parameters_gs.rds")
gics

gs_model_BM <- readRDS("./analysis/testing_models/testing_gs/gs_model_BM.rds")

gs_model_BM


