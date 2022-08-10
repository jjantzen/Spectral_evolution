#summarizing simulation
library(ape)
library(phytools)
library(spectrolab)
library(stringr)
library(dplyr)
library(mvMORPH)

#read in objects
fit_2_ou <- readRDS("./data/simulations/rds/fit_ou_full_ou.rds")
fit_2_bm <- readRDS("./data/simulations/rds/fit_ou_full_bm.rds")
fit_4_ou <- readRDS("./data/simulations/rds/fit_lam_full_ou.rds")
fit_4_bm <- readRDS("./data/simulations/rds/fit_lam_full_bm.rds")

fit_3_ou <- readRDS("./data/simulations/rds/fit_eb_full_ou.rds")
fit_3_bm <- readRDS("./data/simulations/rds/fit_eb_full_bm.rds")
fit_1_ou <- readRDS("./data/simulations/rds/fit_bm_full_ou.rds")
fit_1_bm <- readRDS("./data/simulations/rds/fit_bm_full_bm.rds")

ou_spectra <- readRDS("./data/simulations/rds/ou_sim_spectra.rds")
ou_spectra_df <- readRDS("./data/simulations/rds/ou_sim_spectra_df.rds")
bm_spectra <- readRDS("./data/simulations/rds/bm_sim_spectra.rds")
bm_spectra_df <- readRDS("./data/simulations/rds/bm_sim_spectra_df.rds")


GIC(fit_1_bm)
GIC(fit_1_ou) 
GIC(fit_2_bm)
GIC(fit_2_ou)

GIC(fit_3_bm) 
GIC(fit_3_ou)
GIC(fit_4_bm)
GIC(fit_4_ou)

gic_list_bm <- c()
for (i in 1:length(fit_bm)){
  gic <- GIC(fit_bm[[i]][[1]])$GIC
  gic_list_bm <- c(gic_list_bm, gic)
}

gic_list_ou <- c()
for (i in 1:length(fit_ou)){
  gic <- GIC(fit_ou[[i]][[1]])$GIC
  gic_list_ou <- c(gic_list_ou, gic)
}

gic_list_eb <- c()
for (i in 1:length(fit_eb)){
  gic <- GIC(fit_eb[[i]][[1]])$GIC
  gic_list_eb <- c(gic_list_eb, gic)
}

gic_list_lam <- c()
for (i in 1:length(fit_lam)){
  gic <- GIC(fit_lam[[i]][[1]])$GIC
  gic_list_lam <- c(gic_list_lam, gic)
}

gic_list_lam
#ou model is best for the BM-generated data


GIC(fit_bm) #best
GIC(fit_ou)
GIC(fit_eb)
GIC(fit_lam)


#do manova
aov_1_ac <- manova.gls(fit_2_ou, nperm=99, test="Wilks", verbose=TRUE)


regime2<-as.vector(c(rep("Woody",6),rep("Nonwoody",4))); names(regime2)<-sm_tree2$tip.label
sm_tree2<-make.simmap(sm_tree2, regime2, model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Woody","Nonwoody")

regime2

#univariate - single trait evolved
sim_trait_bm_single <- rTraitCont(sm_tree2, model = "BM", sigma = 0.1, root.value = 1) #, alpha = 1, theta = 0, ancestor = FALSE, root.value = 0, ...)

sim_trait_ou_single <- rTraitCont(sm_tree2, model = "OU", alpha = 1, theta = 1.5, sigma = 0.1, root.value = 1) #ancestor = FALSE, root.value = 0, ...)


#str(ou_spectra)

jpeg("./output/simulated_data_parameter1_bm.jpg", width = 10, height = 10, res = 300, units = "in")
phenogram(sm_tree2, sim_trait_bm_single)
title("BM-evolved parameter")
dev.off()



jpeg("./output/simulated_data_parameter1_ou.jpg", width = 10, height = 10, res = 300, units = "in")
phenogram(sm_tree2, sim_trait_ou_single)
title("OU-evolved parameter")
dev.off()


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
spectral_data_bm_single <- matrix(NA, nrow = 10, ncol = 2101) #make this the right dimensions then plug in
colnames(spectral_data_bm_single) <- single_spec_list_bm[[1]]@wavelength
row.names(spectral_data_bm_single) <- names(sim_trait_bm_single)
for (i in 1:length(single_spec_list_bm)){
  spectral_data_bm_single[i,]<- single_spec_list_bm[[i]]@spectra@spectra_ma[1,]
}

#Single OU data
spectral_data_ou_sm <- matrix(NA, nrow = 10, ncol = 2101) #make this the right dimensions then plug in
colnames(spectral_data_ou_sm) <- sm_spec_list_ou[[1]]@wavelength
row.names(spectral_data_ou_sm) <- names(sim_trait_ou_single)
for (i in 1:length(sm_spec_list_ou)){
  spectral_data_ou_sm[i,]<- sm_spec_list_ou[[i]]@spectra@spectra_ma[1,]
}


#plot
jpeg("./output/sim_spectra_bm_and_ou.jpg", height = 8, width = 8, units = "in", res = 400)
plot(single_spec_list_bm[[1]],ylim=c(0,0.6),xlim=c(400,2500), type = "l", col = "brown")
points(400:2500,single_spec_list_bm[[2]]@spectra@spectra_ma,type="l",col="brown")
points(400:2500,single_spec_list_bm[[3]]@spectra@spectra_ma,type="l",col="brown")
points(400:2500,single_spec_list_bm[[4]]@spectra@spectra_ma,type="l",col="brown")
points(400:2500,single_spec_list_bm[[5]]@spectra@spectra_ma,type="l",col="brown")
points(400:2500,single_spec_list_bm[[6]]@spectra@spectra_ma,type="l",col="brown")
points(400:2500,single_spec_list_bm[[7]]@spectra@spectra_ma,type="l",col="brown")
points(400:2500,single_spec_list_bm[[8]]@spectra@spectra_ma,type="l",col="brown")
points(400:2500,single_spec_list_bm[[9]]@spectra@spectra_ma,type="l",col="brown")
points(400:2500,single_spec_list_bm[[10]]@spectra@spectra_ma,type="l",col="brown")

points(400:2500,sm_spec_list_ou[[2]]@spectra@spectra_ma,type="l",col="green")
points(400:2500,sm_spec_list_ou[[2]]@spectra@spectra_ma,type="l",col="green")
points(400:2500,sm_spec_list_ou[[3]]@spectra@spectra_ma,type="l",col="green")
points(400:2500,sm_spec_list_ou[[4]]@spectra@spectra_ma,type="l",col="green")
points(400:2500,sm_spec_list_ou[[5]]@spectra@spectra_ma,type="l",col="green")
points(400:2500,sm_spec_list_ou[[6]]@spectra@spectra_ma,type="l",col="green")
points(400:2500,sm_spec_list_ou[[7]]@spectra@spectra_ma,type="l",col="green")
points(400:2500,sm_spec_list_ou[[8]]@spectra@spectra_ma,type="l",col="green")
points(400:2500,sm_spec_list_ou[[9]]@spectra@spectra_ma,type="l",col="green")
points(400:2500,sm_spec_list_ou[[10]]@spectra@spectra_ma,type="l",col="green")

dev.off()


#create empty dataframe
gic_comparison_bm_df <- as.data.frame(matrix(ncol = 5, nrow = 100))
colnames(gic_comparison_bm_df) <- c("nrep", "BM", "OU", "EB", "LAM")

length(fit_2_ou)

#iterate over reps
bm_fit_bm <- list()
#BM
for (i in 1:100){
  gic_comparison_bm_df$nrep[i] <- i
  #bm <- mvgls(trait ~ habit, data=simulated_bm[[i]], tree=simmap_tree, model="BM") 
  bm <- fit_2_bm
  bm_fit_bm[[length(bm_fit_bm)+1]] <- list(summary(bm))
  gic_comparison_bm_df$BM[i] <- GIC(bm)$GIC
}
i <- 1
for (i in 1:100){
  gic_comparison_bm_df$nrep[i] <- i
  bm <- bm_fit_bm[[i]]
  gic_comparison_bm_df$BM[i] <- GIC(bm)$GIC
}


