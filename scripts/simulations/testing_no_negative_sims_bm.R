library(OUwie)
library(hsdar)

#read tree
simmap_tree <- readRDS("./analysis/from_cluster/full/100_taxa_simmap_tree.rds")

#sim data
sim_trait_bm <- OUwie.sim(simmap_tree, simmap.tree = TRUE, sigma.sq = c(0.1, 0.1), alpha = c(1e-10,1e-10), theta = c(0,0), theta0 = 2)

#make sure no values are negative
if (all(sim_trait_bm$X >= 0) ){
  sim_trait_bm$X <- sim_trait_bm$X
} else {
  sim_trait_bm$X[which(sim_trait_bm$X < 0)] <- 0 
}

#make spectra
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

spectral_data_bm <- matrix(NA, nrow = 100, ncol = 2101) #make this the right dimensions then plug in
colnames(spectral_data_bm) <- spec_list_bm[[1]]@wavelength
row.names(spectral_data_bm) <- sim_trait_bm$Genus_species
for (i in 1:length(spec_list_bm)){
  spectral_data_bm[i,]<- spec_list_bm[[i]]@spectra@spectra_ma[1,]
}


#check spectra for negatives
spectral_data_bm
str(spectral_data_bm)

if (any(spectral_data_bm < 0)){
  spectral_data_bm[which(spectral_data_bm < 0)] <- 0 
} else {
  spectral_data_bm <- spectral_data_bm
}


