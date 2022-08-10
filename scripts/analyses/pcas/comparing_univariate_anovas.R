#intercept model PCAs
library(mvMORPH)
library(phytools)
library(dplyr)
library(ggplot2)
library(OUwie)
library(nortest)
library(nlme)
source("./scripts/plot_phylo_pca_edited_function.R")
#load plot_pca_phylo function

#read in data
data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")
tree_myc <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

model_list <- readRDS("./analysis/pca_analysis/best_intercept_models_for_pcas_92sp.rds")

myc_named <- setNames(data_spectra$myc, rownames(data_spectra$spectra))

simmap_tree <- readRDS("./analysis/pca_analysis/myc_simmap_trees.rds")

i <- 4
#comparing different anova/regressions for pcas
pca_best_model <- mvgls.pca(model_list[[i]], plot = FALSE)  

output_aov_per <- data.frame(matrix(nrow=10, ncol = 7))
colnames(output_aov_per)[1] <- "iteration"
class(output_aov_per[,1]) <- "numeric"
colnames(output_aov_per)[2] <- "PC_axis"
class(output_aov_per[,2]) <- "numeric"
colnames(output_aov_per)[3] <- "F_stat"
class(output_aov_per[,3]) <- "numeric"
colnames(output_aov_per)[4] <- "p_value"
class(output_aov_per[,4]) <- "numeric"
colnames(output_aov_per)[5] <- "df"
class(output_aov_per[,5]) <- "numeric"
colnames(output_aov_per)[6] <- "sum_sq"
class(output_aov_per[,6]) <- "numeric"
colnames(output_aov_per)[7] <- "mean_sq"
class(output_aov_per[,7]) <- "numeric"

#run ouwie models

#get data
j <- 1
data <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named, X=as.numeric(pca_best_model$scores[,j]))
  
#non phylo methods
pc_aov <- aov(X ~ myc, data = data)
summary_aov <- summary(pc_aov) #0.00314 significant #this is non phylogenetic on phyloPCA

pc_lm <- lm(X~myc,data=data)
summary(pc_lm) #0.00314 significant


##actual phylogenetic analysis
simmap_tree
vcv_matrix <- vcv.phylo(simmap_tree[[i]], cor=TRUE)
colnames(data)[1] <- "Species"

#1 in this is lambda
correlation=corPagel(1, pruned_tree_pig, form = ~Species)

#phylo (PCM)
phy_gls <- gls(X~myc,data=data, correlation=corPagel(1,simmap_tree[[i]], form = ~species))
summary(phy_gls) #0.0505 not quite significant

phy_gls_lam <- gls(X~myc,data=data, correlation=corPagel(model_list[[4]]$param,simmap_tree[[i]], form = ~species))
sum_gls <- summary(phy_gls_lam) #0.0505 - same no matter what lambda is 
str(sum_gls)

#lambda
sum_gls[[1]]$corStruct[[1]]
#AIC
sum_gls$AIC
#pvalue
sum_gls$tTable[[8]]



plot(phy_gls_lam$fitted, phy_gls_lam$residuals)

qqnorm(phy_gls_lam$residuals)
plot(fitted(phy_gls_lam), phy_gls_lam$residuals)

abline(0,0, col = "red")

plot(density(phy_gls_lam$residuals))

plot(data$X, phy_gls_lam$residuals)


model_list[[4]]$param

data$X

#calculate phylogenetic signal in pc axis
pc_1 <- setNames(data$X, data$species)
sig1 <- phylosig(simmap_tree[[i]], pc_1, method = "K", test=TRUE, nsim=999) 

j <- 2
data <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named, X=as.numeric(pca_best_model$scores[,j]))

pc_2 <- setNames(data$X, data$species)
sig2 <- phylosig(simmap_tree[[i]], pc_2, method = "K", test=TRUE, nsim=999) 
sig2

simmap_tree[[i]]$tip.label


################
#test for normality of residuals
lillie.test(residuals(phy_gls)) #0.0001215
lillie.test(chol(solve(vcv(simmap_tree[[i]])))%*%residuals(phy_gls)) #0.0001517

#and for non phylo residuals
lillie.test(residuals(pc_lm))  #0.001154
lillie.test(chol(solve(vcv(simmap_tree[[i]])))%*%residuals(pc_lm)) #0.0002027
