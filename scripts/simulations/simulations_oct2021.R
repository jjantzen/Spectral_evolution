#Evolution of leaf spectra
#Script for simulations to test methods for analysis 

#load libraries
library(phytools) #tree simulations
library(hsdar) #spectra simulations
library(parallel) #for simulations of trait data?
library(mvMORPH) #for modeling
library(RPANDA) #for modeling
library(spectrolab)

######################
#1 Generate phylogeny
######################
#List of taxa, seqeuence data, alignment, make trees (here simulate)
#Make a distribution of trees for analysis
#Test with standard error for intraspecific variation

#set seed for reproducibility
set.seed(123)

#use approximately the same number of taxa as dataset - 100

# #different functions for simulating trees
# trees <- rmtree(100, n= 100) 
tree <- pbtree(n=100) #redo with smaller dataset (eg 10x100)

sm_tree <- pbtree(n=10)

# ult_tree <- rcoal(100)

# #make list of ultrametric trees
# ult_trees <- c()
# 
# for (i in 1:25){
#   if (i == 1) {
#     tree <- rcoal(100)
#     ult_trees <- tree
#   } else if (i > 1) {
#     tree <- rcoal(100)
#     ult_trees <- c(ult_trees, tree)
#   }
# }
# 

######################
#2 Evolve parameters for spectra
######################
#need 100 spectra for the 100 species (plus variation within species?)

##simulate traits under different processes

#random spectra - no brownian motion on phylogeny
#one optimum spectra - OU model - specify one optima
#three optima spectra - OU model - specify three optima

#PROSPECT

#parameters that can vary
#cab = varies from 10-400?
#nit = varies from 1-1.5
#car = varies from 10-200
#anth = 0.8-1.2
#cbrown
#cw
#cm = varies from 0.01 to 0.02

#focus on evolving anth and nit (centered around 1)

#first just use single tree
tree

#set regime (two states)
regime1<-as.vector(c(rep("Woody",65),rep("Nonwoody",35))); names(regime1)<-tree$tip.label

#plot states onto tree
tree<-make.simmap(tree, regime1, model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Woody","Nonwoody")

# Plot of the phylogeny for illustration
plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

#evidently strong correlation between woodiness and phylogeny (how do I vary this?)

## Simulate trait evolution according to a bivariate "BMM" model
# Number of traits - why not evolve just one? - just N, or Cw - effect on entire range of wavelengths
ntraits<-2

# Number of simulated (pairs of) traits
nsim<-1

# Rates matrices for the "Woody" and the "Nonwoody" regimes
#this is a 2x2 matrix with covariances for two traits (anth and nit), for two states
sigma_bm<-list(Woody=matrix(c(0.02,0.06,0.06,0.02),2), Nonwoody=matrix(c(0.01,0.04,0.04,0.01),2))

#ancestral states for each trait (centered around 1)
theta_bm<-c(1,1)

#simulate both traits under multiple rate BM - try single rate
sim_data_bm <- mvSIM(tree, model="BMM", nsim=nsim, param=list(ntraits=ntraits,sigma=sigma_bm,theta=theta_bm))

# Rates matrices for the "Woody" and the "Nonwoody" regimes
#this is a 2x2 matrix with covariances for two traits (anth and nit), for two states
sigma_ou<-matrix(c(0.02,0.06,0.06,0.02),2)

#ancestral states for each trait (centered around 1)
theta_ou<-c(1,1)
alpha <- c(2,4) #but should they be the same?

#simulate both traits under single OU model
sim_data_ou <- mvSIM(tree, model="OU1", nsim=nsim, param=list(ntraits=ntraits,sigma=sigma_ou,theta=c(1,1)))#,alpha = alpha

#model spectra using evolved parameters
spec_list_bm <- c()

for (i in 1:100){
  ind_param <- sim_data_bm[i,]
  cab <-  50
  nit <- ind_param[1]
  car <- 100
  anth <- ind_param[2]
  cm <- 0.01
  spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")
  if (i == 1) {
    spec_list_bm <- spec
  } else {
    spec_list_bm <- c(spec_list_bm, spec)
  }
}

#model spectra using evolved parameters
spec_list_ou <- c()

for (i in 1:100){
  ind_param <- sim_data_ou[i,]
cab <-  50
nit <- ind_param[1]
car <- 100
anth <- ind_param[2]
cm <- 0.01
spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")
if (i == 1) {
  spec_list_ou <- spec
} else {
  spec_list_ou <- c(spec_list_ou, spec)
}
}

#plot to visualize
plot(spec_list_ou[[1]],ylim=c(0,0.5),xlim=c(400,2500))
points(400:2500,spec_list_ou[[2]]@spectra@spectra_ma,type="l",col="red")
points(400:2500,spec_list_ou[[3]]@spectra@spectra_ma,type="l",col="blue")
points(400:2500,spec_list_ou[[4]]@spectra@spectra_ma,type="l",col="green")
points(400:2500,spec_list_ou[[5]]@spectra@spectra_ma,type="l",col="yellow")
points(400:2500,spec_list_ou[[6]]@spectra@spectra_ma,type="l",col="brown")

#To summarize:
#two rates for woody/nonwoody spectra
#one optimum (not different for woody/not woody)
#use mvMORPH to see if I get the same result


######################
#3 Model spectra
######################
#using mvMORPH, RPANDA fit_t_pl, phyl.pca_pl for PCA
rownames(sim_data_ou)[1]
#data
spectral_data_ou <- matrix(NA, nrow = 100, ncol = 2101) #make this the right dimensions then plug in
colnames(spectral_data_ou) <- spec_list_ou[[1]]@wavelength
row.names(spectral_data_ou) <- row.names(sim_data_ou)
for (i in 1:length(spec_list_ou)){
  spectral_data_ou[i,]<- spec_list_ou[[i]]@spectra@spectra_ma[1,]
}

spectral_data_bm <- matrix(NA, nrow = 100, ncol = 2101) #make this the right dimensions then plug in
colnames(spectral_data_bm) <- spec_list_bm[[1]]@wavelength
row.names(spectral_data_bm) <- row.names(sim_data_bm)
for (i in 1:length(spec_list_bm)){
  spectral_data_bm[i,]<- spec_list_bm[[i]]@spectra@spectra_ma[1,]
}

spectral_data_bm[1:10,1001:1010]
spectral_data_bm[1:10,1:10] #why inf?

growth_form <- matrix(regime1)

#combine data into one dataframe
data_bm <- list(trait=spectral_data_bm, habit=as.factor(growth_form))
data_ou <- list(trait=spectral_data_ou, habit=as.factor(growth_form))


length(data_bm$habit)

str(data_bm$trait)

dim(data_bm$trait)
length(tree$tip.label)

# Fit the multivariate linear model
fit_1 <- mvgls(trait ~ habit, data=data_bm, tree=tree, model="BM") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_2 <- mvgls(trait ~ habit, data=data_bm, tree=tree, model="OU") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_3 <- mvgls(trait ~ habit, data=data_bm, tree=tree, model="EB") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity

fit_1b <- mvgls(trait ~ 1, data=data_bm, tree=tree, model="BM") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_2b <- mvgls(trait ~ 1, data=data_bm, tree=tree, model="OU") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_3b <- mvgls(trait ~ 1, data=data_bm, tree=tree, model="EB") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity

fit_1 <- readRDS("./data/simulations/fit_1_bm_data.rds")
fit_2 <- readRDS("./data/simulations/fit_2_bm_data.rds")

GIC(fit_1) #-5511260
GIC(fit_2) #-5512857

row.names(data_bm$trait)
tree$tip.label

nrow(model_fr)!=length(tree$tip.label)

nrow(data_bm$trait)

#want each row to be an individual (eg spectra item in list) and each column to be wavelength
#then add dataframe for woody/not (factors)


########
#fit set of models (these aren't high dimensional)
fit1 <- mvgls()

model_fit<-mvOU(tree, ou_sim[[1]],model="OUM")
model_fit<-mvOU(tree, ou_sim[[1]],model="OU1")
model_fit<-mvBM(tree,bm_data[[1]],model="BM1")
model_fit<-mvBM(tree,bm_data[[1]],model="BMM")

#do model test (GIC)

#three datasets simulated under different conditions

#spec_list - random actual spectra
#ou_data - ou model for two states (OUM), 2100 traits but not spectra
#bm_data - bm model for two states (BM1), 2100 traits but not spectra



