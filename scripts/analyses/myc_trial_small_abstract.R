#Do trial of myc data on small dataset

library(phytools)
library(mvMORPH)
library(spectrolab)
library(dplyr)
library(reshape)
library(OUwie)
library(ggplot2)

#read trees
trees <- readRDS("./data/tidy/new_trees_matched_spectra.rds")

#read spectra
spectra <- readRDS("./data/tidy/spectra_not_reordered_to_tree.rds")

#read myc data
data <- read.csv("./data/predictors/myc_prelim_version.csv", stringsAsFactors = FALSE)
data <- data[,c(1,2)]

colnames(data)

data <- na.omit(data)
nrow(data)

#read gf data
data_2 <- readRDS("./data/tidy/new_growth_forms_matched_spectra.rds")
data_2 <- data_2[,c(2,3,4,5)]

#combine into one growth form column
data_3 <- data_2 %>% 
  pivot_longer(cols = c(2:4)) %>% 
  dplyr::filter(value == 1)

#read lp data
lp <- read.csv("./data/predictors/lp_and_rk.csv", stringsAsFactors = FALSE)
lp <- lp[,c(1,29)]
colnames(lp)[1] <- "species"

#match data sets
data_3 <- data_3[which(data_3$species_names %in% data$Species),]
lp <- lp[which(lp$species %in% data$Species),]
nrow(lp)
nrow(data_3)

#get unique entry per species
data_3 #start here combine shrub and tree etc for those with dups

dups <- data_3 %>% group_by(species_names) %>% dplyr::summarize(n=n()) %>% dplyr::filter(n == 2)

#manually edit these?
data_3$value[which(data_3$species_names == "Betula populifolia" & data_3$name == "Shrub")] <- "0"
data_3$value[which(data_3$species_names == "Frangula alnus" & data_3$name == "Tree")] <- "0"
data_3$value[which(data_3$species_names == "Sorbus decora" & data_3$name == "Shrub")] <- "0"
data_3$value[which(data_3$species_names == "Sorbus americana" & data_3$name == "Shrub")] <- "0"
data_3$value[which(data_3$species_names == "Rhamnus cathartica" & data_3$name == "Shrub")] <- "0"
data_3$value[which(data_3$species_names == "Prunus pensylvanica" & data_3$name == "Shrub")] <- "0"
data_3$value[which(data_3$species_names == "Prunus nigra" & data_3$name == "Shrub")] <- "0"
data_3$value[which(data_3$species_names == "Crataegus monogyna" & data_3$name == "Shrub")] <- "0"
data_3$value[which(data_3$species_names == "Acer pensylvanicum" & data_3$name == "Shrub")] <- "0"

#combine into one growth form column
data_4 <- data_3 %>% 
  #pivot_longer(cols = c(2:4)) %>% 
  dplyr::filter(value == 1) %>% 
  dplyr::select(species_names, name)

nrow(data_4)

#prune to match (missing data in myc)
tree_keep <- trees[[1]]
tree_drops <- tree_keep$tip.label[-which(tree_keep$tip.label %in% data$Species)]

tree_sm <- drop.tip(tree_keep, tree_drops)
tree_sm

#prune spectra
spectra_df <- as.data.frame(spectra)
spectra_df_sm <- as.matrix(spectra_df[which(rownames(spectra_df) %in% data$Species),])
#spectra_sm <- as_spectra(spectra_df_sm, name_idx = 1, meta_idxs = c(2:8))

#objects to keep
data_4
tree_sm
spectra_df_sm
lp
data

#make list of dataframes and model
#make data list for model
data_spectra <- list(spectra=spectra_df_sm, myc = data$Myc_binary, gf = data_4$name, lp = lp$Leaf_persistence)

str(data_spectra)

#save object
saveRDS(data_spectra, "./data/tidy/myc_trial_list_for_analysis.rds")
saveRDS(tree_sm, "./data/tidy/myc_trial_tree_for_analysis.rds")

# #make model
# myc_bm <- mvgls(spectra ~ myc*gf*lp, data = data_spectra, tree = tree_sm, model = "BM", error = TRUE)
# saveRDS(myc_bm, "./myc_model_trial_bm.rds")
# 
# myc_ou <- mvgls(spectra ~ myc*gf*lp, data = data_spectra, tree = tree_sm, model = "OU", error = TRUE)
# saveRDS(myc_ou, "./myc_model_trial_ou.rds")
# 
# myc_eb <- mvgls(spectra ~ myc*gf*lp, data = data_spectra, tree = tree_sm, model = "EB", error = TRUE)
# saveRDS(myc_eb, "./myc_model_trial_eb.rds")
# 
# myc_lam <- mvgls(spectra ~ myc*gf*lp, data = data_spectra, tree = tree_sm, model = "lambda", error = TRUE)
# saveRDS(myc_lam, "./myc_model_trial_lambda.rds")


####Do on cluster - see script there
# #create output dataframe
# df_output <- data.frame(matrix(nrow=4, ncol = 3))
# colnames(df_output) <- c("model", "GIC", "parameter")
# 
# #run for each of 4 models of evolution
# models <- c("BM", "OU", "EB", "lambda")
# 
# for (j in 1:length(models)){
#   myc_bm <- mvgls(spectra ~ myc*gf*lp, data = data_spectra, tree = tree_sm, model=models[j], error = TRUE)
#   #assign output
#   df_output$model[j] <- models[j]
#   df_output$parameter[j] <- model$param
#   df_output$GIC[j] <- GIC(model)$GIC
#   #save model itself
#   saveRDS(model, paste0("./myc_models/myc_model_trial_", models[j], ".rds"))
# }
# 
# saveRDS(df_output, "./myc_models/model_parameters_myc_trial.rds")
# 

data_spectra <- readRDS("./data/tidy/myc_trial_list_for_analysis.rds")
tree_sm <- readRDS("./data/tidy/myc_trial_tree_for_analysis.rds")

#Read in model results

myc_model_trial_BM <- readRDS("./analysis/testing_models/testing_myc/myc_model_trial_BM.rds")
myc_model_trial_OU <- readRDS("./analysis/testing_models/testing_myc/myc_model_trial_OU.rds")
myc_model_trial_EB <- readRDS("./analysis/testing_models/testing_myc/myc_model_trial_EB.rds")
myc_model_trial_lam <- readRDS("./analysis/testing_models/testing_myc/myc_model_trial_lambda.rds")
fit_ou_intercept <- readRDS("./data/physignal/fit_ou_intercept.rds")

myc_gics <- readRDS("./analysis/testing_models/testing_myc/model_parameters_myc_trial.rds")
myc_gics  #best is lambda then BM

#Map trait to tree - make simmap tree for predictor
myc_named <- setNames(data_spectra$myc, rownames(data_spectra$spectra))

simmap_tree <- make.simmap(tree_sm, myc_named, model="SYM", nsim=1)
plot(simmap_tree)

#assign states to nodes
simmap_tree$node.label<-getStates(simmap_tree,"nodes")

#Do PCA for model
pca_best_model <- mvgls.pca(myc_model_trial_BM, plot = TRUE) #check which model to use - is BM/lambda

pca_lambda_model <- mvgls.pca(myc_model_trial_lam, plot = TRUE) #check which model to use - is BM/lambda

pca_intercept_model <- mvgls.pca(fit_ou_intercept, plot = TRUE) #check which model to use - is BM/lambda

#Extract PCs
data_pc1 <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named,
                       X=as.numeric(pca_best_model$scores[,1]))

data_pc2 <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named,
                       X=as.numeric(pca_best_model$scores[,2]))

data_pc3 <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named,
                       X=as.numeric(pca_best_model$scores[,3]))

data_pc4 <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named,
                       X=as.numeric(pca_best_model$scores[,4]))

data_pc20 <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named,
                        X=as.numeric(pca_best_model$scores[,20]))

#from lambda model (actual best)
data_pc1_lam <- data.frame(species=rownames(pca_lambda_model$scores),myc=myc_named,
                       X=as.numeric(pca_lambda_model$scores[,1]))

data_pc2_lam <- data.frame(species=rownames(pca_lambda_model$scores),myc=myc_named,
                       X=as.numeric(pca_lambda_model$scores[,2]))

data_pc3_lam <- data.frame(species=rownames(pca_lambda_model$scores),myc=myc_named,
                       X=as.numeric(pca_lambda_model$scores[,3]))

data_pc4_lam <- data.frame(species=rownames(pca_lambda_model$scores),myc=myc_named,
                       X=as.numeric(pca_lambda_model$scores[,4]))

#from lambda model (actual best)
data_pc1_OU <- data.frame(species=rownames(pca_OU_model$scores),myc=myc_named,
                           X=as.numeric(pca_OU_model$scores[,1]))

data_pc2_OU <- data.frame(species=rownames(pca_OU_model$scores),myc=myc_named,
                           X=as.numeric(pca_OU_model$scores[,2]))

data_pc3_OU <- data.frame(species=rownames(pca_OU_model$scores),myc=myc_named,
                           X=as.numeric(pca_OU_model$scores[,3]))

data_pc4_OU <- data.frame(species=rownames(pca_OU_model$scores),myc=myc_named,
                           X=as.numeric(pca_OU_model$scores[,4]))

# #from intercept model (actual best) - not the same number of taxa 
# #would need intercept model for 87 taxa
# data_pc1_int <- data.frame(species=rownames(pca_intercept_model$scores),myc=myc_named,
#                           X=as.numeric(pca_intercept_model$scores[,1]))
# 
# data_pc2_int <- data.frame(species=rownames(pca_intercept_model$scores),myc=myc_named,
#                           X=as.numeric(pca_intercept_model$scores[,2]))
# 
# data_pc3_int <- data.frame(species=rownames(pca_intercept_model$scores),myc=myc_named,
#                           X=as.numeric(pca_intercept_model$scores[,3]))
# 
# data_pc4_int <- data.frame(species=rownames(pca_intercept_model$scores),myc=myc_named,
#                           X=as.numeric(pca_intercept_model$scores[,4]))
#Model PCs (univariate)

#model options: model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA","TrendyM","TrendyMS")
pc1_BM1_test <- OUwie(simmap_tree, data_pc1, model = "BM1", simmap.tree = TRUE) #single rate
pc1_BMS_test <- OUwie(simmap_tree, data_pc1, model = "BMS", simmap.tree = TRUE) #multiple rates
pc1_OUM_test <- OUwie(simmap_tree, data_pc1, model = "OUM", simmap.tree = TRUE) #multiple optima
pc1_OU1_test <- OUwie(simmap_tree, data_pc1, model = "OU1", simmap.tree = TRUE) #single optimum

aic <- setNames(c(pc1_BM1_test$AIC, pc1_BMS_test$AIC, pc1_OU1_test$AIC, pc1_OUM_test$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic 
aic.w(aic) #no, single optimum, weight = 0.729

#second PC axis 
pc2_BM1_test <- OUwie(simmap_tree, data_pc2, model = "BM1", simmap.tree = TRUE) #single rate
pc2_BMS_test <- OUwie(simmap_tree, data_pc2, model = "BMS", simmap.tree = TRUE) #multiple rates
pc2_OUM_test <- OUwie(simmap_tree, data_pc2, model = "OUM", simmap.tree = TRUE) #multiple optima
pc2_OU1_test <- OUwie(simmap_tree, data_pc2, model = "OU1", simmap.tree = TRUE) #single optimum

aic_2 <- setNames(c(pc2_BM1_test$AIC, pc2_BMS_test$AIC,pc2_OU1_test$AIC, pc2_OUM_test$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_2
aic.w(aic_2) #yes, corresponds to myc oum weight = 0.817

#third PC axis 
pc3_BM1_test <- OUwie(simmap_tree, data_pc3, model = "BM1", simmap.tree = TRUE) #single rate
pc3_BMS_test <- OUwie(simmap_tree, data_pc3, model = "BMS", simmap.tree = TRUE) #multiple rates
pc3_OUM_test <- OUwie(simmap_tree, data_pc3, model = "OUM", simmap.tree = TRUE) #multiple optima
pc3_OU1_test <- OUwie(simmap_tree, data_pc3, model = "OU1", simmap.tree = TRUE) #single optimum

aic_3 <- setNames(c(pc3_BM1_test$AIC, pc3_BMS_test$AIC, pc3_OU1_test$AIC, pc3_OUM_test$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_3
aic.w(aic_3) #also oum - weight = 0.821

#fourth PC axis 
pc4_BM1_test <- OUwie(simmap_tree, data_pc4, model = "BM1", simmap.tree = TRUE) #single rate
pc4_BMS_test <- OUwie(simmap_tree, data_pc4, model = "BMS", simmap.tree = TRUE) #multiple rates
pc4_OUM_test <- OUwie(simmap_tree, data_pc4, model = "OUM", simmap.tree = TRUE) #multiple optima
pc4_OU1_test <- OUwie(simmap_tree, data_pc4, model = "OU1", simmap.tree = TRUE) #single optimum

aic_4 <- setNames(c(pc4_BM1_test$AIC, pc4_BMS_test$AIC, pc4_OU1_test$AIC, pc4_OUM_test$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_4
aic.w(aic_4) #back to ou1 - weight = 0.711


#and for lambda
#model options: model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA","TrendyM","TrendyMS")
pc1_BM1_test_lam <- OUwie(simmap_tree, data_pc1_lam, model = "BM1", simmap.tree = TRUE) #single rate
pc1_BMS_test_lam <- OUwie(simmap_tree, data_pc1_lam, model = "BMS", simmap.tree = TRUE) #multiple rates
pc1_OUM_test_lam <- OUwie(simmap_tree, data_pc1_lam, model = "OUM", simmap.tree = TRUE) #multiple optima
pc1_OU1_test_lam <- OUwie(simmap_tree, data_pc1_lam, model = "OU1", simmap.tree = TRUE) #single optimum

aic <- setNames(c(pc1_BM1_test_lam$AIC, pc1_BMS_test_lam$AIC, pc1_OU1_test_lam$AIC, pc1_OUM_test_lam$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic 
aic.w(aic) #no, single optimum, weight = 0.724

#second PC axis 
pc2_BM1_test_lam <- OUwie(simmap_tree, data_pc2_lam, model = "BM1", simmap.tree = TRUE) #single rate
pc2_BMS_test_lam <- OUwie(simmap_tree, data_pc2_lam, model = "BMS", simmap.tree = TRUE) #multiple rates
pc2_OUM_test_lam <- OUwie(simmap_tree, data_pc2_lam, model = "OUM", simmap.tree = TRUE) #multiple optima
pc2_OU1_test_lam <- OUwie(simmap_tree, data_pc2_lam, model = "OU1", simmap.tree = TRUE) #single optimum

aic_2 <- setNames(c(pc2_BM1_test_lam$AIC, pc2_BMS_test_lam$AIC,pc2_OU1_test_lam$AIC, pc2_OUM_test_lam$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_2
aic.w(aic_2) #no, corresponds to myc ou1 weight = 0.589 (oum weight 0.411)

#third PC axis 
pc3_BM1_test_lam <- OUwie(simmap_tree, data_pc3_lam, model = "BM1", simmap.tree = TRUE) #single rate
pc3_BMS_test_lam <- OUwie(simmap_tree, data_pc3_lam, model = "BMS", simmap.tree = TRUE) #multiple rates
pc3_OUM_test_lam <- OUwie(simmap_tree, data_pc3_lam, model = "OUM", simmap.tree = TRUE) #multiple optima
pc3_OU1_test_lam <- OUwie(simmap_tree, data_pc3_lam, model = "OU1", simmap.tree = TRUE) #single optimum

aic_3 <- setNames(c(pc3_BM1_test_lam$AIC, pc3_BMS_test_lam$AIC, pc3_OU1_test_lam$AIC, pc3_OUM_test_lam$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_3
aic.w(aic_3) #also ou1 - weight = 0.639

#fourth PC axis 
pc4_BM1_test_lam <- OUwie(simmap_tree, data_pc4_lam, model = "BM1", simmap.tree = TRUE) #single rate
pc4_BMS_test_lam <- OUwie(simmap_tree, data_pc4_lam, model = "BMS", simmap.tree = TRUE) #multiple rates
pc4_OUM_test_lam <- OUwie(simmap_tree, data_pc4_lam, model = "OUM", simmap.tree = TRUE) #multiple optima
pc4_OU1_test_lam <- OUwie(simmap_tree, data_pc4_lam, model = "OU1", simmap.tree = TRUE) #single optimum

aic_4 <- setNames(c(pc4_BM1_test_lam$AIC, pc4_BMS_test_lam$AIC, pc4_OU1_test_lam$AIC, pc4_OUM_test_lam$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_4
aic.w(aic_4) #back to ou1 - weight = 0.701

#and for OU
#model options: model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA","TrendyM","TrendyMS")
pc1_BM1_test_OU <- OUwie(simmap_tree, data_pc1_OU, model = "BM1", simmap.tree = TRUE) #single rate
pc1_BMS_test_OU <- OUwie(simmap_tree, data_pc1_OU, model = "BMS", simmap.tree = TRUE) #multiple rates
pc1_OUM_test_OU <- OUwie(simmap_tree, data_pc1_OU, model = "OUM", simmap.tree = TRUE) #multiple optima
pc1_OU1_test_OU <- OUwie(simmap_tree, data_pc1_OU, model = "OU1", simmap.tree = TRUE) #single optimum

aic <- setNames(c(pc1_BM1_test_OU$AIC, pc1_BMS_test_OU$AIC, pc1_OU1_test_OU$AIC, pc1_OUM_test_OU$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic 
aic.w(aic) #no, single optimum, weight = 0.705

#second PC axis 
pc2_BM1_test_OU <- OUwie(simmap_tree, data_pc2_OU, model = "BM1", simmap.tree = TRUE) #single rate
pc2_BMS_test_OU <- OUwie(simmap_tree, data_pc2_OU, model = "BMS", simmap.tree = TRUE) #multiple rates
pc2_OUM_test_OU <- OUwie(simmap_tree, data_pc2_OU, model = "OUM", simmap.tree = TRUE) #multiple optima
pc2_OU1_test_OU <- OUwie(simmap_tree, data_pc2_OU, model = "OU1", simmap.tree = TRUE) #single optimum

aic_2 <- setNames(c(pc2_BM1_test_OU$AIC, pc2_BMS_test_OU$AIC,pc2_OU1_test_OU$AIC, pc2_OUM_test_OU$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_2
aic.w(aic_2) #yes, corresponds to myc oum weight = 0.693

#third PC axis 
pc3_BM1_test_OU <- OUwie(simmap_tree, data_pc3_OU, model = "BM1", simmap.tree = TRUE) #single rate
pc3_BMS_test_OU <- OUwie(simmap_tree, data_pc3_OU, model = "BMS", simmap.tree = TRUE) #multiple rates
pc3_OUM_test_OU <- OUwie(simmap_tree, data_pc3_OU, model = "OUM", simmap.tree = TRUE) #multiple optima
pc3_OU1_test_OU <- OUwie(simmap_tree, data_pc3_OU, model = "OU1", simmap.tree = TRUE) #single optimum

aic_3 <- setNames(c(pc3_BM1_test_OU$AIC, pc3_BMS_test_OU$AIC, pc3_OU1_test_OU$AIC, pc3_OUM_test_OU$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_3
aic.w(aic_3) #also oum - weight = 0.543

#fourth PC axis 
pc4_BM1_test_OU <- OUwie(simmap_tree, data_pc4_OU, model = "BM1", simmap.tree = TRUE) #single rate
pc4_BMS_test_OU <- OUwie(simmap_tree, data_pc4_OU, model = "BMS", simmap.tree = TRUE) #multiple rates
pc4_OUM_test_OU <- OUwie(simmap_tree, data_pc4_OU, model = "OUM", simmap.tree = TRUE) #multiple optima
pc4_OU1_test_OU <- OUwie(simmap_tree, data_pc4_OU, model = "OU1", simmap.tree = TRUE) #single optimum

aic_4 <- setNames(c(pc4_BM1_test_OU$AIC, pc4_BMS_test_OU$AIC, pc4_OU1_test_OU$AIC, pc4_OUM_test_OU$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_4
aic.w(aic_4) #back to ou1 - weight = 0.650


#plot pc2 vs pc3
col.group <- data_spectra$myc
col.group <- gsub("AM", "orange", col.group)
col.group <- gsub("EM", "blue", col.group)

jpeg("./output/PCAs/myc_phylo_pc2_3.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(pca_best_model, myc_model_trial_BM, axes = c(2,3), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
dev.off()

#plot pc1 vs pc4
jpeg("./output/PCAs/myc_phylo_pc1_4.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(pca_best_model, myc_model_trial_BM, axes = c(1,4), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
dev.off()

#get wavelengths for pc2 and pc3
pca_best_model

####conduct anova on pc 2 and 3 for myc
str(data_pc2)
pc2_aov <- aov(X ~ myc, data = data_pc2)

pc3_aov <- aov(X ~ myc, data = data_pc3)

s_aov <- summary(pc2_aov)
s_aov[[1]]$`Pr(>F)`[1]
s_aov[[1]]$`F value`[1]
s_aov[[1]]$Df
s_aov[[1]]$`Sum Sq`
s_aov[[1]]$`Mean Sq`

pc1_aov <- aov(X ~ myc, data = data_pc1)

pc4_aov <- aov(X ~ myc, data = data_pc4)

#do multiple anovas (need multiple test correction) for many pcs (beyond top 4)
data_spectra$lp[which(data_spectra$lp == "deciduous?")] <- "deciduous"
model_lm_v1 <- lm(data_pc1$X ~ data_spectra$myc*data_spectra$gf*data_spectra$lp)

model_lm_v2 <- lm(data_pc2$X ~ data_spectra$myc*data_spectra$gf*data_spectra$lp)

model_lm_v2b <- lm(data_pc2$X ~ data_spectra$myc*data_spectra$gf)

model_lm_v2c <- lm(data_pc2$X ~ data_spectra$myc*data_spectra$lp)

model_lm_v3 <- lm(data_pc3$X ~ data_spectra$myc*data_spectra$gf*data_spectra$lp)

model_lm_v4 <- lm(data_pc4$X ~ data_spectra$myc*data_spectra$gf*data_spectra$lp)

summary(model_lm_v1)

summary(lm(data_pc1$X ~ data_spectra$myc))

#choose how many pcs to keep - how?
var_explained <- 100 * pca_best_model$values / sum(pca_best_model$values)

# #one example axis for verifying
# tot <- sum(pca_best_model$values)
# valX <- round(pca_best_model$values[1] * 100/tot, digits = 2)
# valX

#create scree plot
jpeg("./output/PCAs/myc_trial_scree_plot.jpg", res = 300, width = 10, height = 10, units = "in")
qplot(c(1:20), var_explained[c(1:20)]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 100)
dev.off()

var_explained[c(1:20)]

#top 6 axes explain 99% of variance
var_explained[1]+var_explained[2]+var_explained[3]+var_explained[4]+var_explained[5]+var_explained[6]

#####calculate loadings
#loading can be defined as Eigenvector multiplied by the square root of the Eigen value - for standardized, or just eigenvector
pca_best_model_loadings <- pca_best_model
pca_best_model_loadings$loadings_calc <- pca_best_model_loadings$vectors * sqrt(pca_best_model_loadings$values)
pca_best_model_loadings$loadings_calc[1:10]

###plot loadings

#add names to matrices - need wavelength not species
dimnames(pca_best_model_loadings$loadings_calc) <- list(c(400:2400),c(400:2400))
dimnames(pca_best_model_loadings$vectors) <- list(c(400:2400),c(1:2001))

n.pc1 <- ifelse(pca_best_model_loadings$vectors[,2] > 0, yes=-0.01, no=pca_best_model_loadings$vectors[,2]-0.01)

#just plot eigenvectors
jpeg("./output/PCAs/myc_trial_pc2_vectors.jpg", res = 300, width = 10, height = 10, units = "in")
b2 <- barplot(pca_best_model_loadings$vectors[,2], main="PC 2 Eigenvectors Plot", las=2)#
abline(h=0)
dev.off()

jpeg("./output/PCAs/myc_trial_pc3_vectors.jpg", res = 300, width = 10, height = 10, units = "in")
b2 <- barplot(pca_best_model_loadings$vectors[,3], main="PC 3 Eigenvectors Plot", las=2)#
abline(h=0)
dev.off()

jpeg("./output/PCAs/myc_trial_pc6_vectors.jpg", res = 300, width = 10, height = 10, units = "in")
b2 <- barplot(pca_best_model_loadings$vectors[,6], main="PC 6 Eigenvectors Plot", las=2)#
abline(h=0)
dev.off()

#make one plot with series of eigenvector plots

#create dataframe from vectors with columns as pc axis and rows as wavelength

vectors_df <- as.data.frame(pca_best_model_loadings$vectors)
str(vectors_df)
colnames(vectors_df) <- c(1:2001)
rownames(vectors_df)

# ggplot(vectors_df[,c(1,2)], aes(x = rownames(vectors_df), y = vectors_df[,1]))+
#   geom_bar(stat = "identity")# +
#   facet_wrap(colnames(vectors_df), ncol = 2, scales = "free")
# 
# str(pca_best_model_loadings$vectors)

  
vectors_df <- pca_best_model_loadings$vectors %>% 
  as.data.frame() 

vectors_df$wavelength <- c(400:2400)
colnames(vectors_df)

vectors_for_plotting <- vectors_df %>% tidyr::gather(., pc_axis, value, -wavelength)

class(vectors_for_plotting$pc_axis) <- "numeric"

pc_labs <- paste("PC", c(1:12))
names(pc_labs) <- c(1:12)

jpeg("./output/PCAs/myc_trial_top12_vectors.jpg", res = 400, width = 15, height = 15, units = "in")
ggplot(vectors_for_plotting[which(vectors_for_plotting$pc_axis %in% c(1:12)),], aes(x=wavelength,y=value)) + 
  geom_bar(stat = "identity")+
  facet_wrap(vars(pc_axis), ncol = 3, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
  theme(strip.text.x = element_text(size = 15))+
  xlab("Wavelength (nm)")+
  ylab("Vector unit")
dev.off()

pc_labs_30 <- paste("PC", c(1:30))
names(pc_labs_30) <- c(1:30)

jpeg("./output/PCAs/myc_trial_top30_vectors.jpg", res = 400, width = 20, height = 20, units = "in")
ggplot(vectors_for_plotting[which(vectors_for_plotting$pc_axis %in% c(1:30)),], aes(x=wavelength,y=value)) + 
  geom_bar(stat = "identity")+
  facet_wrap(vars(pc_axis), ncol = 5, labeller = labeller(pc_axis = pc_labs_30))+# scales = "free"
  theme(strip.text.x = element_text(size = 15))+
  xlab("Wavelength (nm)")+
  ylab("Vector unit")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
dev.off()



#text(x=b2, y=n.pc1, labels=names(pca_best_model_loadings$vectors[,2]), adj=1, srt=90, xpd=FALSE)

# 
# b2 <- barplot(pca_best_model_loadings$loadings_calc[,2], main="PC 2 Loadings Plot", las=2)
# abline(h=0.0001)


# #extract loadings with highest values
# str(pca_best_model_loadings$vectors)
# 
# #get highest value eigenvectors and names
# dimnames(pca_best_model_loadings$vectors)[[1]][which(pca_best_model_loadings$vectors > 0.01)] #???
# 

#eigenvalues - what do they mean?
b3 <- barplot(pca_best_model_loadings$values, main="PC 2 Eigenvalues Plot", las=2)

#plot actual spectra by myc association
spec_for_plot <- as_spectra(data_spectra$spectra)

meta(spec_for_plot)$myc <- data_spectra$myc

plot(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], lwd = 0.75, lty = 1, col = "blue", main = "Spectra by mycorrhizal association")
plot(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], lwd = 0.75, lty = 1, col = "red", add = TRUE)

jpeg("./output/spectra_by_myc.jpg", res = 400, width = 10, height = 10, units = "in")
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], total_prob = 0.8, col = rgb(1,0,0,0.25), border = TRUE, main = "Spectra by mycorrhizal association (80% quantile)")
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], total_prob = 0.8, col = rgb(0,0,1,0.25), border = TRUE, add = TRUE)
legend(2000, 0.5, legend=c("AM", "EM"), fill=c(rgb(1,0,0,0.25), rgb(0,0,1,0.25)), cex=2)#lty=1,col
dev.off()


##### look at manovas for high dimensional stuff
model_parameters_myc_simple <- readRDS("./analysis/testing_models/testing_myc/model_parameters_myc_simple.rds")

OU_myc_only_manova <- readRDS("./analysis/testing_models/testing_myc/OU_myc_only_manova.rds")

OU_lp_only_manova <- readRDS("./analysis/testing_models/testing_myc/OU_lp_only_manova.rds")

OU_gf_only_manova <- readRDS("./analysis/testing_models/testing_myc/OU_gf_only_manova.rds")

myc_only_model <- readRDS("./analysis/testing_models/testing_myc/myc_model_simple.rds")
myc_only_model
lp_only_model <- readRDS("./analysis/testing_models/testing_myc/lp_model_simple.rds")
gf_only_model <- readRDS("./analysis/testing_models/testing_myc/gf_model_simple.rds")
lp_only_model
gf_only_model
OU_myc_only_manova

#summary(myc_model_trial_BM)


#compare two different models (one vs three predictors) - doesn't work with these models
anova(gf_only_model, myc_only_model)

stationary(myc_only_model)

#none are significant - why?
#does coding matter (text vs binary 0/1)
#which var-covar matrix is best for pgls
#what is pgls actually doing - look at textbook again
#look into how to compare hierarchical models?
#signal in trait vs residuals

#calculate phylosignal in myc, gf and lp traits
model_parameters_myc_trial

summary(myc_model_trial_lambda)
str(myc_model_trial_lambda)
myc_model_trial_lambda$terms

library(mvMORPH)
bm_myc_manova_2
OU_myc_only_manova
OU_gf_only_manova
OU_lp_only_manova
model_parameters_myc_trial



colnames(combo_data)
combo_data$Species
library(mvMORPH)
summary(myc_model_lambdaiteration1)

myc_model_BMiteration1
myc_model_lambdaiteration1$param

