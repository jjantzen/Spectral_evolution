#Evolution of leaf spectra
#Script for simulations to test methods for analysis 

#load libraries
library(phytools) #tree simulations
library(hsdar) #spectra simulations
#library(parallel) #for simulations of trait data?
library(mvMORPH) #for modeling
library(RPANDA) #for modeling
library(spectrolab)
library(OUwie)

set.seed(123)

######################
#BIG DATASET
######################
#1 Generate phylogeny
######################

tree <- pbtree(n=100) 

saveRDS(tree, "./data/simulations/full/100_taxa_tree_sim.tre")

######################
#2 Set regime for tree
######################

regime <- as.vector(c(rep("Woody",60),rep("Nonwoody",40))); names(regime) <- tree$tip.label
simmap_tree<-make.simmap(tree, regime, model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Woody","Nonwoody")

######################
#Evolve spectra
######################

#do multiple simulations for distributions
#write function with simmap tree as input

#convert to four functions - different models each

simulate_spectra_bm <- function(simmap_tree){
  ######################
  #3 Evolve single parameter for spectra
  ######################
  #simplest models
  sim_trait_bm <- OUwie.sim(simmap_tree, simmap.tree = TRUE, sigma.sq = c(0.1, 0.1), alpha = c(1e-10,1e-10), theta = c(0,0), theta0 = 1)
  
  ######################
  #4 Model spectra from evolved parameter
  ######################
  
  #Single BM 
  spec_list_bm <- c()
  
  for (i in 1:100){
    ind_param <- sim_trait_bm[i,2]
    cab <-  50
    nit <- ind_param
    car <- 100
    anth <- 1
    cm <- 0.01
    spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")
    if (i == 1) {
      spec_list_bm <- spec
    } else {
      spec_list_bm <- c(spec_list_bm, spec)
    }
  }
  
  ######################
  #5 Prep for modeling spectra
  ######################
  
  #Single BM data
  spectral_data_bm <- matrix(NA, nrow = 100, ncol = 2101) #make this the right dimensions then plug in
  colnames(spectral_data_bm) <- spec_list_bm[[1]]@wavelength
  row.names(spectral_data_bm) <- sim_trait_bm$Genus_species
  for (i in 1:length(spec_list_bm)){
    spectral_data_bm[i,]<- spec_list_bm[[i]]@spectra@spectra_ma[1,]
  }
  
  #assign predictor from simulated trait data
  growth_form <- matrix(regime)
  
  #combine data into one dataframe
  data_bm <- list(trait=spectral_data_bm, habit=as.factor(growth_form))
  
  return(data_bm)
}

simulate_spectra_ou <- function(simmap_tree){
  ######################
  #3 Evolve single parameter for spectra
  ######################
  #simplest models
  sim_trait_ou <- OUwie.sim(simmap_tree, simmap.tree = TRUE, sigma.sq = c(0.1, 0.1), alpha = c(1, 1), theta0 = 1, theta = c(1.5, 1.5))
  
  ######################
  #4 Model spectra from evolved parameter
  ######################
  
  #Single OU
  spec_list_ou <- c()
  
  for (i in 1:100){
    ind_param <- sim_trait_ou[i,2]
    cab <-  50
    nit <- ind_param
    car <- 100
    anth <- 1
    cm <- 0.01
    spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")
    if (i == 1) {
      spec_list_ou <- spec
    } else {
      spec_list_ou <- c(spec_list_ou, spec)
    }
  }
  
  ######################
  #5 Prep for modeling spectra
  ######################
  
  #Single OU data
  spectral_data_ou <- matrix(NA, nrow = 100, ncol = 2101) #make this the right dimensions then plug in
  colnames(spectral_data_ou) <- spec_list_ou[[1]]@wavelength
  row.names(spectral_data_ou) <- sim_trait_ou$Genus_species
  for (i in 1:length(spec_list_ou)){
    spectral_data_ou[i,]<- spec_list_ou[[i]]@spectra@spectra_ma[1,]
  }
  
  #assign predictor from simulated trait data
  growth_form <- matrix(regime)
  
  #combine data into one dataframe
  data_ou <- list(trait=spectral_data_ou, habit=as.factor(growth_form))
  
  return(data_ou)
}

simulate_spectra_bmm <- function(simmap_tree){
  ######################
  #3 Evolve single parameter for spectra
  ######################
  
  #do more complex models
  sim_trait_bmm <- OUwie.sim(simmap_tree, simmap.tree = TRUE, sigma.sq = c(0.1, 0.4), alpha = c(1,1), theta = c(0, 0), theta0 = 1)
  
  ######################
  #4 Model spectra from evolved parameter
  ######################
  
  #Multi BM
  spec_list_bmm <- c()
  
  for (i in 1:100){
    ind_param <- sim_trait_bmm[i,2]
    cab <-  50
    nit <- ind_param
    car <- 100
    anth <- 1
    cm <- 0.01
    spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")
    if (i == 1) {
      spec_list_bmm <- spec
    } else {
      spec_list_bmm <- c(spec_list_bmm, spec)
    }
  }
  
  ######################
  #5 Prep for modeling spectra
  ######################
  
  #Multi BM data
  spectral_data_bmm <- matrix(NA, nrow = 100, ncol = 2101) #make this the right dimensions then plug in
  colnames(spectral_data_bmm) <- spec_list_bmm[[1]]@wavelength
  row.names(spectral_data_bmm) <- sim_trait_bmm$Genus_species
  for (i in 1:length(spec_list_bmm)){
    spectral_data_bmm[i,]<- spec_list_bmm[[i]]@spectra@spectra_ma[1,]
  }
  
  #assign predictor from simulated trait data
  growth_form <- matrix(regime)
  
  #combine data into one dataframe
  data_bmm <- list(trait=spectral_data_bmm, habit=as.factor(growth_form))
 
  return(data_bmm)
}

simulate_spectra_oum <- function(simmap_tree){
  ######################
  #3 Evolve single parameter for spectra
  ######################
  #simplest models
 
  #do more complex models
  sim_trait_oum <- OUwie.sim(simmap_tree, simmap.tree = TRUE, sigma.sq = c(0.1, 0.1), alpha = c(1, 1), theta0 = 1, theta = c(1, 2))
  
  ######################
  #4 Model spectra from evolved parameter
  ######################
  
  #Multi OU
  spec_list_oum <- c()
  
  for (i in 1:100){
    ind_param <- sim_trait_oum[i,2]
    cab <-  50
    nit <- ind_param
    car <- 100
    anth <- 1
    cm <- 0.01
    spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")
    if (i == 1) {
      spec_list_oum <- spec
    } else {
      spec_list_oum <- c(spec_list_oum, spec)
    }
  }
  
  ######################
  #5 Prep for modeling spectra
  ######################
  
  #Multi OU data
  spectral_data_oum <- matrix(NA, nrow = 100, ncol = 2101) #make this the right dimensions then plug in
  colnames(spectral_data_oum) <- spec_list_oum[[1]]@wavelength
  row.names(spectral_data_oum) <- sim_trait_oum$Genus_species
  for (i in 1:length(spec_list_oum)){
    spectral_data_oum[i,]<- spec_list_oum[[i]]@spectra@spectra_ma[1,]
  }
  
  #assign predictor from simulated trait data
  growth_form <- matrix(regime)
  
  #combine data into one dataframe
  data_oum <- list(trait=spectral_data_oum, habit=as.factor(growth_form))
  
  return(data_oum)
}

#set dimension of lists
num_ports <- c(100, 2)

#create empty lists and run simulate functions
simulated_bm <- Reduce(function(x, y) replicate(y, x, F), rev(num_ports), NULL)

for (i in 1:100){
  output <- simulate_spectra_bm(simmap_tree)
  simulated_bm[[i]] <- output
}

#create empty lists and run simulate functions
simulated_ou <- Reduce(function(x, y) replicate(y, x, F), rev(num_ports), NULL)

for (i in 1:100){
  output <- simulate_spectra_ou(simmap_tree)
  simulated_ou[[i]] <- output
}

#create empty lists and run simulate functions
simulated_bmm <- Reduce(function(x, y) replicate(y, x, F), rev(num_ports), NULL)

for (i in 1:100){
  output <- simulate_spectra_bmm(simmap_tree)
  simulated_bmm[[i]] <- output
}

#create empty lists and run simulate functions
simulated_oum <- Reduce(function(x, y) replicate(y, x, F), rev(num_ports), NULL)

for (i in 1:100){
  output <- simulate_spectra_oum(simmap_tree)
  simulated_oum[[i]] <- output
}

str(simulated_bm)
#str(output)

saveRDS(simulated_bm, "./data/simulations/full/bm_data_sim.tre")
saveRDS(simulated_ou, "./data/simulations/full/ou_data_sim.tre")
saveRDS(simulated_bmm, "./data/simulations/full/bmm_data_sim.tre")
saveRDS(simulated_oum, "./data/simulations/full/oum_data_sim.tre")

######################
#6 Model spectra
######################

#dataframe version

#BM data
bm_fit_bm <- list()
ou_fit_bm <- list()
eb_fit_bm <- list()
lam_fit_bm <- list()

#create empty dataframe
gic_comparison_bm_df <- as.data.frame(matrix(ncol = 5, nrow = 100))
colnames(gic_comparison_bm_df) <- c("nrep", "BM", "OU", "EB", "LAM")

#iterate over reps
for (i in 1:100){
  gic_comparison_bm_df$nrep[i] <- i
  bm <- mvgls(trait ~ habit, data=simulated_bm[[i]], tree=simmap_tree, model="BM") 
  bm_fit_bm[[length(bm_fit_bm)+1]] <- list(summary(bm))
  gic_comparison_bm_df$BM[i] <- GIC(bm)$GIC
  ou <- mvgls(trait ~ habit, data=simulated_bm[[i]], tree=simmap_tree, model="OU")
  ou_fit_bm[[length(ou_fit_bm)+1]] <- list(summary(ou))
  gic_comparison_bm_df$OU[i] <- GIC(ou)$GIC
  eb <- mvgls(trait ~ habit, data=simulated_bm[[i]], tree=simmap_tree, model="EB")
  eb_fit_bm[[length(eb_fit_bm)+1]] <- list(summary(eb))
  gic_comparison_bm_df$EB[i] <- GIC(eb)$GIC
  lam <- mvgls(trait ~ habit, data=simulated_bm[[i]], tree=simmap_tree, model="lambda")
  lam_fit_bm[[length(lam_fit_bm)+1]] <- list(summary(lam))
  gic_comparison_bm_df$LAM[i] <- GIC(lam)$GIC
}


#save output
saveRDS(bm_fit_bm, "./data/simulations/full/bm_fit_bm.rds")
saveRDS(ou_fit_bm, "./data/simulations/full/ou_fit_bm.rds")
saveRDS(eb_fit_bm, "./data/simulations/full/eb_fit_bm.rds")
saveRDS(lam_fit_bm, "./data/simulations/full/lam_fit_bm.rds")

saveRDS(gic_comparison_bm_df, "./data/simulations/full/gic_comparison_bm_df.rds")

#OU data
bm_fit_ou <- list()
ou_fit_ou <- list()
eb_fit_ou <- list()
lam_fit_ou <- list()

#create empty dataframe
gic_comparison_ou_df <- as.data.frame(matrix(ncol = 5, nrow = 100))
colnames(gic_comparison_ou_df) <- c("nrep", "BM", "OU", "EB", "LAM")


for (i in 1:100){
  gic_comparison_ou_df$nrep[i] <- i
  bm <- mvgls(trait ~ habit, data=simulated_ou[[i]], tree=simmap_tree, model="BM") 
  bm_fit_ou[[length(bm_fit_ou)+1]] <- list(summary(bm))
  gic_comparison_ou_df$BM[i] <- GIC(bm)$GIC
  ou <- mvgls(trait ~ habit, data=simulated_ou[[i]], tree=simmap_tree, model="OU")
  ou_fit_ou[[length(ou_fit_ou)+1]] <- list(summary(ou))
  gic_comparison_ou_df$OU[i] <- GIC(ou)$GIC
  eb <- mvgls(trait ~ habit, data=simulated_ou[[i]], tree=simmap_tree, model="EB")
  eb_fit_ou[[length(eb_fit_ou)+1]] <- list(summary(eb))
  gic_comparison_ou_df$EB[i] <- GIC(eb)$GIC
  lam <- mvgls(trait ~ habit, data=simulated_ou[[i]], tree=simmap_tree, model="lambda")
  lam_fit_ou[[length(lam_fit_ou)+1]] <- list(summary(lam))
  gic_comparison_ou_df$LAM[i] <- GIC(lam)$GIC
}

saveRDS(bm_fit_ou, "./data/simulations/full/bm_fit_ou.rds")
saveRDS(ou_fit_ou, "./data/simulations/full/ou_fit_ou.rds")
saveRDS(eb_fit_ou, "./data/simulations/full/eb_fit_ou.rds")
saveRDS(lam_fit_ou, "./data/simulations/full/lam_fit_ou.rds")

saveRDS(gic_comparison_ou_df, "./data/simulations/full/gic_comparison_ou_df.rds")


#BMM data
bm_fit_bmm <- list()
ou_fit_bmm <- list()
eb_fit_bmm <- list()
lam_fit_bmm <- list()

#create empty dataframe
gic_comparison_bmm_df <- as.data.frame(matrix(ncol = 5, nrow = 100))
colnames(gic_comparison_bmm_df) <- c("nrep", "BM", "OU", "EB", "LAM")


for (i in 1:100){
  gic_comparison_bmm_df$nrep[i] <- i
  bm <- mvgls(trait ~ habit, data=simulated_bmm[[i]], tree=simmap_tree, model="BM") 
  bm_fit_bmm[[length(bm_fit_bmm)+1]] <- list(summary(bm))
  gic_comparison_bmm_df$BM[i] <- GIC(bm)$GIC
  ou <- mvgls(trait ~ habit, data=simulated_bmm[[i]], tree=simmap_tree, model="OU")
  ou_fit_bmm[[length(ou_fit_bmm)+1]] <- list(summary(ou))
  gic_comparison_bmm_df$OU[i] <- GIC(ou)$GIC
  eb <- mvgls(trait ~ habit, data=simulated_bmm[[i]], tree=simmap_tree, model="EB")
  eb_fit_bmm[[length(eb_fit_bmm)+1]] <- list(summary(eb))
  gic_comparison_bmm_df$EB[i] <- GIC(eb)$GIC
  lam <- mvgls(trait ~ habit, data=simulated_bmm[[i]], tree=simmap_tree, model="lambda")
  lam_fit_bmm[[length(lam_fit_bmm)+1]] <- list(summary(lam))
  gic_comparison_bmm_df$LAM[i] <- GIC(lam)$GIC
}

saveRDS(bm_fit_bmm, "./data/simulations/full/bm_fit_bmm.rds")
saveRDS(ou_fit_bmm, "./data/simulations/full/ou_fit_bmm.rds")
saveRDS(eb_fit_bmm, "./data/simulations/full/eb_fit_bmm.rds")
saveRDS(lam_fit_bmm, "./data/simulations/full/lam_fit_bmm.rds")

saveRDS(gic_comparison_bmm_df, "./data/simulations/full/gic_comparison_bmm_df.rds")


#OUM data
bm_fit_oum <- list()
ou_fit_oum <- list()
eb_fit_oum <- list()
lam_fit_oum <- list()

#create empty dataframe
gic_comparison_oum_df <- as.data.frame(matrix(ncol = 5, nrow = 100))
colnames(gic_comparison_oum_df) <- c("nrep", "BM", "OU", "EB", "LAM")


for (i in 1:100){
  gic_comparison_oum_df$nrep[i] <- i
  bm <- mvgls(trait ~ habit, data=simulated_oum[[i]], tree=simmap_tree, model="BM") 
  bm_fit_oum[[length(bm_fit_oum)+1]] <- list(summary(bm))
  gic_comparison_oum_df$BM[i] <- GIC(bm)$GIC
  ou <- mvgls(trait ~ habit, data=simulated_oum[[i]], tree=simmap_tree, model="OU")
  ou_fit_oum[[length(ou_fit_oum)+1]] <- list(summary(ou))
  gic_comparison_oum_df$OU[i] <- GIC(ou)$GIC
  eb <- mvgls(trait ~ habit, data=simulated_oum[[i]], tree=simmap_tree, model="EB")
  eb_fit_oum[[length(eb_fit_oum)+1]] <- list(summary(eb))
  gic_comparison_oum_df$EB[i] <- GIC(eb)$GIC
  lam <- mvgls(trait ~ habit, data=simulated_oum[[i]], tree=simmap_tree, model="lambda")
  lam_fit_oum[[length(lam_fit_oum)+1]] <- list(summary(lam))
  gic_comparison_oum_df$LAM[i] <- GIC(lam)$GIC
}

saveRDS(bm_fit_oum, "./data/simulations/full/bm_fit_oum.rds")
saveRDS(ou_fit_oum, "./data/simulations/full/ou_fit_oum.rds")
saveRDS(eb_fit_oum, "./data/simulations/full/eb_fit_oum.rds")
saveRDS(lam_fit_oum, "./data/simulations/full/lam_fit_oum.rds")

saveRDS(gic_comparison_oum_df, "./data/simulations/full/gic_comparison_oum_df.rds")



# #Do model comparison (GIC)
# 
# highest <- c()
# 
# for (i in 1:nrow(dataframe)){
#   highest[i] <- colnames(dataframe)[which(max(dataframe[i,]))]
# }
# 
