#testing OUWie

library(OUwie)
library(mvMORPH)
library(factoextra)
library(phytools)
library(l1ou) #which one do I want to use? - probably OUwie? - or mvMORPH mvOU function?

#load model for pcas
ou_fit <- readRDS("./data/physignal/fit_ou_intercept.rds")

#run pca
pca_ou <- mvgls.pca(ou_fit, plot = TRUE)

#read any predictors
combo_data <- readRDS("./data/tidy/new_combo_data_matched_spectra.rds")

#get predictors as named list - do this for another predictor as well (leaf persistence and drought_bin)
wood <- combo_data$Woody
names(wood) <- combo_data$Species

#read trees
new_trees <- readRDS("./data/tidy/new_trees_matched_spectra.rds")

tree <- new_trees[[1]]

#make simmap tree for predictor
simmap_tree <- make.simmap(tree, wood, model="SYM", nsim=1)
plot(simmap_tree)

## make analysis input data.frame
wood_2 <- setNames(combo_data$Woody, combo_data$Species)

#need to make sure order lines up by species
#pc1 <- as.data.frame(pca_ou$scores[,1], row.names = rownames(pca_ou$scores)) #I think - use scores not values or vectors??

data_pc1 <- data.frame(species=rownames(pca_ou$scores),woody=wood_2,
                 X=as.numeric(pca_ou$scores[,1]))

data_pc2 <- data.frame(species=rownames(pca_ou$scores),woody=wood_2,
                       X=as.numeric(pca_ou$scores[,2]))

data_pc3 <- data.frame(species=rownames(pca_ou$scores),woody=wood_2,
                       X=as.numeric(pca_ou$scores[,3]))

data_pc4 <- data.frame(species=rownames(pca_ou$scores),woody=wood_2,
                       X=as.numeric(pca_ou$scores[,4]))

data_pc20 <- data.frame(species=rownames(pca_ou$scores),woody=wood_2,
                       X=as.numeric(pca_ou$scores[,20]))

#assign states to nodes
simmap_tree$node.label<-getStates(simmap_tree,"nodes")

#fit OUwie models
#do I need to put trees into simmap.tree format first? (using predictors?) - I think so 
#data is pc axis

#model options: model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA","TrendyM","TrendyMS")
pc1_BM1_test <- OUwie(simmap_tree, data_pc1, model = "BM1", simmap.tree = TRUE) #single rate
pc1_BMS_test <- OUwie(simmap_tree, data_pc1, model = "BMS", simmap.tree = TRUE) #multiple rates
pc1_OUM_test <- OUwie(simmap_tree, data_pc1, model = "OUM", simmap.tree = TRUE) #multiple optima
pc1_OU1_test <- OUwie(simmap_tree, data_pc1, model = "OU1", simmap.tree = TRUE) #single optimum

pc1_BM1_test #AICc = 444.5138, lnL = -220.1951, BIC = 449.6004
pc1_BMS_test #AICc = 420.3369, lnL = -207.0435, BIC = 427.9024
pc1_OUM_test #AICc = 363.2179, lnL = -177.3984, BIC = 373.2176
pc1_OU1_test #AICc = 363.9947, lnL = -178.8723, BIC = 371.5602

aic <- setNames(c(pc1_BM1_test$AIC, pc1_BMS_test$AIC, pc1_OU1_test$AIC, pc1_OUM_test$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic 
aic.w(aic) # almost all weight favours OUM model

#second PC axis 
pc2_BM1_test <- OUwie(simmap_tree, data_pc2, model = "BM1", simmap.tree = TRUE) #single rate
pc2_BMS_test <- OUwie(simmap_tree, data_pc2, model = "BMS", simmap.tree = TRUE) #multiple rates
pc2_OUM_test <- OUwie(simmap_tree, data_pc2, model = "OUM", simmap.tree = TRUE) #multiple optima
pc2_OU1_test <- OUwie(simmap_tree, data_pc2, model = "OU1", simmap.tree = TRUE) #single optimum

pc2_BM1_test #AICc = 357.4211, lnL = -176.6487, BIC = 362.5077
pc2_BMS_test #AICc = 358.3393, lnL = -176.0447, BIC = 365.9048
pc2_OUM_test #AICc = 249.1638, lnL = -120.5819, BIC = 259.5845
pc2_OU1_test #AICc = 249.8748, lnL = -121.9374, BIC = 257.6903

aic_2 <- setNames(c(pc2_BM1_test$AIC, pc2_BMS_test$AIC,pc2_OU1_test$AIC, pc2_OUM_test$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_2
aic.w(aic_2) # almost all weight favours OUM model

#third PC axis 
pc3_BM1_test <- OUwie(simmap_tree, data_pc3, model = "BM1", simmap.tree = TRUE) #single rate
pc3_BMS_test <- OUwie(simmap_tree, data_pc3, model = "BMS", simmap.tree = TRUE) #multiple rates
pc3_OUM_test <- OUwie(simmap_tree, data_pc3, model = "OUM", simmap.tree = TRUE) #multiple optima
pc3_OU1_test <- OUwie(simmap_tree, data_pc3, model = "OU1", simmap.tree = TRUE) #single optimum

pc3_BM1_test #AICc = 197.1416, lnL = -96.57081, BIC = 202.352
pc3_BMS_test #AICc = 192.2644, lnL = -93.13219, BIC = 200.0799
pc3_OUM_test #AICc = 23.50955, lnL = -7.544248, BIC = 33.50918
pc3_OU1_test #AICc = 66.71343, lnL = -30.23171, BIC = 74.27894

aic_3 <- setNames(c(pc3_BM1_test$AIC, pc3_BMS_test$AIC, pc3_OU1_test$AIC, pc3_OUM_test$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_3
aic.w(aic_3) # almost all weight favours OUM model

#fourth PC axis 
pc4_BM1_test <- OUwie(simmap_tree, data_pc4, model = "BM1", simmap.tree = TRUE) #single rate
pc4_BMS_test <- OUwie(simmap_tree, data_pc4, model = "BMS", simmap.tree = TRUE) #multiple rates
pc4_OUM_test <- OUwie(simmap_tree, data_pc4, model = "OUM", simmap.tree = TRUE) #multiple optima
pc4_OU1_test <- OUwie(simmap_tree, data_pc4, model = "OU1", simmap.tree = TRUE) #single optimum

pc4_BM1_test #AICc = 82.34715, lnL = -39.11172, BIC = 87.43378
pc4_BMS_test #AICc = 75.83595, lnL = -34.79297, BIC = 83.40146
pc4_OUM_test #AICc = -116.9693, lnL = 62.6952, BIC = -106.9697
pc4_OU1_test #AICc = -112.6766, lnL = 59.46331, BIC = -105.1111

aic_4 <- setNames(c(pc4_BM1_test$AIC, pc4_BMS_test$AIC, pc4_OU1_test$AIC, pc4_OUM_test$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_4
aic.w(aic_4) # almost all weight favours OUM model

#20th PC axis 
pc20_BM1_test <- OUwie(simmap_tree, data_pc20, model = "BM1", simmap.tree = TRUE) #single rate
pc20_BMS_test <- OUwie(simmap_tree, data_pc20, model = "BMS", simmap.tree = TRUE) #multiple rates
pc20_OUM_test <- OUwie(simmap_tree, data_pc20, model = "OUM", simmap.tree = TRUE) #multiple optima
pc20_OU1_test <- OUwie(simmap_tree, data_pc20, model = "OU1", simmap.tree = TRUE) #single optimum

pc20_BM1_test #AICc = -546.0862, lnL = 275.105, BIC = -540.9996
pc20_BMS_test #AICc = -558.0392, lnL = 282.1446, BIC = -550.4737
pc20_OUM_test #AICc = -731.1345, lnL = 369.7778, BIC = -721.1349
pc20_OU1_test #AICc = -731.7048, lnL = 368.9774, BIC = -724.1393

aic_20 <- setNames(c(pc20_BM1_test$AIC, pc20_BMS_test$AIC, pc20_OU1_test$AIC, pc20_OUM_test$AIC), c("BM1", "BMS", "OU1", "OUM"))
aic_20
aic.w(aic_20) # almost all weight favours OU1 model


#####Switch to mvMORPH function

#do a version of the models using mvMORPH where we compare individual PCs with multivariate models with 1-4 pc axes etc
pc1_OU1_mvmorph <- mvOU(simmap_tree, pca_ou$scores[,1], model = "OU1") #single optimum
pc2_OU1_mvmorph <- mvOU(simmap_tree, pca_ou$scores[,2], model = "OU1") #single optimum
pc1_OUM_mvmorph <- mvOU(simmap_tree, pca_ou$scores[,1], model = "OUM") #multiple optima
pc2_OUM_mvmorph <- mvOU(simmap_tree, pca_ou$scores[,2], model = "OUM") #multiple optima

#mv model
pc1.2_OU1_mvmorph <- mvOU(simmap_tree, pca_ou$scores[,1:2], model = "OU1") #single optimum
pc1.2_OUM_mvmorph <- mvOU(simmap_tree, pca_ou$scores[,1:2], model = "OUM") #multiple optima

#PC1
AIC(pc1_OU1_mvmorph); AIC(pc1_OUM_mvmorph) #not much of a difference but favours OUM

#mv
AIC(pc1.2_OU1_mvmorph); AIC(pc1.2_OUM_mvmorph) #OUM favoured when doing multiple PC axes

#PC2
AIC(pc2_OU1_mvmorph); AIC(pc2_OUM_mvmorph) #more of a difference than for axis 1

#mv model for first 4 pc axes
pc1_to_4_OU1_mvmorph <- mvOU(simmap_tree, pca_ou$scores[,1:4], model = "OU1") #single optimum
pc1_to_4_OUM_mvmorph <- mvOU(simmap_tree, pca_ou$scores[,1:4], model = "OUM") #multiple optima

#PC2
AIC(pc1_to_4_OU1_mvmorph); AIC(pc1_to_4_OUM_mvmorph) #very clearly OUM is better; OU1 result was also unreliable

