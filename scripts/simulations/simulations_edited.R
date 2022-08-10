#Evolution of leaf spectra
#Script for simulations to test methods for analysis 

#load libraries
library(phytools) #tree simulations
library(hsdar) #spectra simulations
#library(parallel) #for simulations of trait data?
library(mvMORPH) #for modeling
library(RPANDA) #for modeling
library(spectrolab)

set.seed(123)

######################
#SMALLER DATASET
######################
#1 Generate phylogeny
######################
# #do not redo simulation of tree - read in saved object

sm_tree <- readRDS("./sim_10_tree.tre")
sm_tree2 <- readRDS("./data/trees/sim_10_tree.tre")

######################
#2 Set regime for tree
######################
regime2<-as.vector(c(rep("Woody",6),rep("Nonwoody",4))); names(regime2)<-sm_tree$tip.label
sm_tree<-make.simmap(sm_tree, regime2, model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Woody","Nonwoody")


######################
#3 Evolve parameters for spectra
######################

#do multiple simulations for distributions

#univariate - single trait evolved
sim_trait_bm_single <- rTraitCont(sm_tree2, model = "BM", sigma = 0.1, root.value = 1) #, alpha = 1, theta = 0, ancestor = FALSE, root.value = 0, ...)

sim_trait_ou_single <- rTraitCont(sm_tree2, model = "OU", alpha = 1, theta = 1.5, sigma = 0.1, root.value = 1) #ancestor = FALSE, root.value = 0, ...)


######################
#4 Model spectra from evolved parameter
######################
#Single BM 
single_spec_list_bm <- c()

for (i in 1:10){
  ind_param <- sim_trait_bm_single[i]
  cab <-  50
  nit <- ind_param
  car <- 100
  anth <- 1
  cm <- 0.01
  spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")
  if (i == 1) {
    single_spec_list_bm <- spec
  } else {
    single_spec_list_bm <- c(single_spec_list_bm, spec)
  }
}

#Single OU
sm_spec_list_ou <- c()

for (i in 1:10){
  ind_param <- sim_trait_ou_single[i]
  cab <-  50
  nit <- ind_param
  car <- 100
  anth <- 1
  cm <- 0.01
  spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")
  if (i == 1) {
    sm_spec_list_ou <- spec
  } else {
    sm_spec_list_ou <- c(sm_spec_list_ou, spec)
  }
}

######################
#5 Prep for modeling spectra
######################

#Do not subset spectra

#str(sim_trait_bm_single)
#Single BM data
spectral_data_bm_single <- matrix(NA, nrow = 10, ncol = 100) #make this the right dimensions then plug in
colnames(spectral_data_bm_single) <- single_spec_list_bm[[1]]@wavelength[1001:1100]
row.names(spectral_data_bm_single) <- names(sim_trait_bm_single)
for (i in 1:length(single_spec_list_bm)){
  spectral_data_bm_single[i,]<- single_spec_list_bm[[i]]@spectra@spectra_ma[1,1001:1100]
}

#Single OU data
spectral_data_ou_sm <- matrix(NA, nrow = 10, ncol = 100) #make this the right dimensions then plug in
colnames(spectral_data_ou_sm) <- sm_spec_list_ou[[1]]@wavelength[1001:1100]
row.names(spectral_data_ou_sm) <- names(sim_trait_ou_single)
for (i in 1:length(sm_spec_list_ou)){
  spectral_data_ou_sm[i,]<- sm_spec_list_ou[[i]]@spectra@spectra_ma[1,1001:1100]
}

#assign predictor from simulated trait data
#What effect does BM of predictor have on results?
growth_form2 <- matrix(regime2)

#combine data into one dataframe
# data_bm_sm <- list(trait=spectral_data_bm_sm, habit=as.factor(growth_form2))
data_bm_single <- list(trait=spectral_data_bm_single, habit=as.factor(growth_form2))
data_ou_sm <- list(trait=spectral_data_ou_sm, habit=as.factor(growth_form2))


######################
#6 Model spectra
######################
fit_bm_sm <- list()
fit_ou_sm <- list()
fit_eb_sm <- list()
fit_lam_sm <- list()

gic_bm_sm <- list()
gic_ou_sm <- list()
gic_eb_sm <- list()
gic_lam_sm <- list()

for (i in 1:10){
  bm <- mvgls(trait ~ habit, data=data_bm_single, tree=sm_tree2, model="BM") 
  fit_bm_sm[[length(fit_bm_sm)+1]] <- list(summary(bm))
  gic_bm_sm[[length(gic_bm_sm)+1]] <- list(GIC(bm))
  ou <- mvgls(trait ~ habit, data=data_bm_single, tree=sm_tree2, model="OU")
  fit_ou_sm[[length(fit_ou_sm)+1]] <- list(summary(ou))
  gic_ou_sm[[length(gic_ou_sm)+1]] <- list(GIC(ou))
  eb <- mvgls(trait ~ habit, data=data_bm_single, tree=sm_tree2, model="EB")
  fit_eb_sm[[length(fit_eb_sm)+1]] <- list(summary(eb))
  gic_eb_sm[[length(gic_eb_sm)+1]] <- list(GIC(eb))
  lam <- mvgls(trait ~ habit, data=data_bm_single, tree=sm_tree2, model="lambda")
  fit_lam_sm[[length(fit_lam_sm)+1]] <- list(summary(lam))
  gic_lam_sm[[length(gic_lam_sm)+1]] <- list(GIC(lam))
}

#Do model comparison (GIC)
capture.output(gic_bm_sm, file = "./gic_bm_sm.txt")
capture.output(gic_ou_sm, file = "./gic_ou_sm.txt")
capture.output(gic_eb_sm, file = "./gic_eb_sm.txt")
capture.output(gic_lam_sm, file = "./gic_lam_sm.txt")

saveRDS(fit_bm_sm, "./fit_bm_sm.rds")
saveRDS(fit_ou_sm, "./fit_ou_sm.rds")
saveRDS(fit_eb_sm, "./fit_eb_sm.rds")
saveRDS(fit_lam_sm, "./fit_lam_sm.rds")

# Fit the multivariate linear model for BM single data (small dataset)
fit_1_single <- mvgls(trait ~ habit, data=data_bm_single, tree=sm_tree2, model="BM") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_2_single <- mvgls(trait ~ habit, data=data_bm_single, tree=sm_tree2, model="OU") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_3_single <- mvgls(trait ~ habit, data=data_bm_single, tree=sm_tree2, model="EB") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_4_single <- mvgls(trait ~ habit, data=data_bm_single, tree=sm_tree2, model="lambda") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity

fit_1b_single <- mvgls(trait ~ 1, data=data_bm_single, tree=sm_tree2, model="BM") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_2b_single <- mvgls(trait ~ 1, data=data_bm_single, tree=sm_tree2, model="OU") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_3b_single <- mvgls(trait ~ 1, data=data_bm_single, tree=sm_tree2, model="EB") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_4b_single <- mvgls(trait ~ 1, data=data_bm_single, tree=sm_tree2, model="lambda") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity

#Do model comparison (GIC)
GIC(fit_1_single) #-20113.58
GIC(fit_2_single) #-20111.58
GIC(fit_3_single) #-20127.32 #lowest
GIC(fit_4_single) #-20111.58

GIC(fit_1b_single) #-22787.11
GIC(fit_2b_single) #-22785.1
GIC(fit_3b_single) #-22790.45 #lowest
GIC(fit_4b_single) #-22785.11

# Fit the multivariate linear model for OU single data (small dataset)
fit_1_ou_sm <- mvgls(trait ~ habit, data=data_ou_sm, tree=sm_tree2, model="BM") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_2_ou_sm <- mvgls(trait ~ habit, data=data_ou_sm, tree=sm_tree2, model="OU") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_3_ou_sm <- mvgls(trait ~ habit, data=data_ou_sm, tree=sm_tree2, model="EB") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_4_ou_sm <- mvgls(trait ~ habit, data=data_ou_sm, tree=sm_tree2, model="lambda") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity

fit_1b_ou_sm <- mvgls(trait ~ 1, data=data_ou_sm, tree=sm_tree2, model="BM") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_2b_ou_sm <- mvgls(trait ~ 1, data=data_ou_sm, tree=sm_tree2, model="OU") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_3b_ou_sm <- mvgls(trait ~ 1, data=data_ou_sm, tree=sm_tree2, model="EB") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity
fit_4b_ou_sm <- mvgls(trait ~ 1, data=data_ou_sm, tree=sm_tree2, model="lambda") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity

#Model comparison (GIC)
GIC(fit_1_ou_sm) #-20133.1
GIC(fit_2_ou_sm) #-20244.84
GIC(fit_3_ou_sm) #-20131.1
GIC(fit_4_ou_sm) #-20278.18 #lowest

GIC(fit_1b_ou_sm) #-22687.81
GIC(fit_2b_ou_sm) #-22895.74 #lowest
GIC(fit_3b_ou_sm) #-22685.81
GIC(fit_4b_ou_sm) #-22890.73

gics_ou <- list(GIC(fit_1_ou_sm)$GIC,GIC(fit_2_ou_sm)$GIC,GIC(fit_3_ou_sm)$GIC, GIC(fit_4_ou_sm)$GIC)
gics_bm <- list(GIC(fit_1_single)$GIC,GIC(fit_2_single)$GIC,GIC(fit_3_single)$GIC, GIC(fit_4_single)$GIC)

capture.output(gics_ou, file = "./gic_output_ou.txt")
capture.output(gics_bm, file = "./gic_output_bm.txt")


saveRDS(fit_1_single, "./fit_1_bm_data.rds")
saveRDS(fit_2_single, "./fit_2_bm_data.rds")
saveRDS(fit_3_single, "./fit_3_bm_data.rds")
saveRDS(fit_4_single, "./fit_4_bm_data.rds")
saveRDS(fit_1b_single, "./fit_1b_bm_data.rds")
saveRDS(fit_2b_single, "./fit_2b_bm_data.rds")
saveRDS(fit_3b_single, "./fit_3b_bm_data.rds")
saveRDS(fit_4b_single, "./fit_4b_bm_data.rds")

saveRDS(fit_1_ou_sm, "./fit_1_ou_data.rds")
saveRDS(fit_2_ou_sm, "./fit_2_ou_data.rds")
saveRDS(fit_3_ou_sm, "./fit_3_ou_data.rds")
saveRDS(fit_4_ou_sm, "./fit_4_ou_data.rds")
saveRDS(fit_1b_ou_sm, "./fit_1b_ou_data.rds")
saveRDS(fit_2b_ou_sm, "./fit_2b_ou_data.rds")
saveRDS(fit_3b_ou_sm, "./fit_3b_ou_data.rds")
saveRDS(fit_4b_ou_sm, "./fit_4b_ou_data.rds")