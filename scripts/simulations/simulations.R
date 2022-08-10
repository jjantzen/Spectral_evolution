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
#Test with polytomies for intraspecific variation
#Test with standard error for intraspecific variation

#use approximately the same number of taxa as dataset - 100

# #different functions for simulating trees
# trees <- rmtree(100, n= 100) 
tree <- pbtree(n=100) 
# ult_tree <- rcoal(100)
# 
# plot(ult_tree, edge.width = 2)
plot(tree)

#make list of ultrametric trees
ult_trees <- c()

for (i in 1:25){
  if (i == 1) {
    tree <- rcoal(100)
    ult_trees <- tree
  } else if (i > 1) {
    tree <- rcoal(100)
    ult_trees <- c(ult_trees, tree)
  }
}

#tips to add polytomy to
tip_labels_trees <- ult_trees_2[[1]]$tip.label
# tip_labels_to_polytomy1 <- tip_labels_trees[1:10]
# tip_labels_to_polytomy2 <- tip_labels_trees[20:30]
# tip_labels_to_polytomy3 <- tip_labels_trees[40:45]

tip_labels_to_polytomy <- tip_labels_trees[c(1:10, 20:30, 40:45)]

#add polytomy to certain tips
poly_trees <- c()

for (i in 1:length(ult_trees_2)){
  test_tree <- ult_trees_2[[i]]
  #for ten species, have five individuals per species
  for (j in 1:10){
    node <- which(test_tree$tip.label == tip_labels_to_polytomy[j])
    for (k in 1:5){
      test_tree <- bind.tip(test_tree, tip.label = paste0(tip_labels_to_polytomy[i], "_", k), where = node, edge.length = 0)    
    }
  }
  #for ten species, have ten individuals per species
  for (j in 11:20){
    node <- which(test_tree$tip.label == tip_labels_to_polytomy[j])
    for (k in 1:10){
      test_tree <- bind.tip(test_tree, tip.label = paste0(tip_labels_to_polytomy[i],"_", k), where = node, edge.length = 0)    
    }
  }
  #for five species, have 20 individuals per species
  for (j in 21:25){
    node <- which(test_tree$tip.label == tip_labels_to_polytomy[j])
    for (k in 1:20){
      test_tree <- bind.tip(test_tree, tip.label = paste0(tip_labels_to_polytomy[i],"_", k), where = node, edge.length = 0)    
    }
  }
  if (i == 1){
    poly_trees <- test_tree
    } else {
    poly_trees <- c(poly_trees, test_tree)
    }
}

#examine trees to make sure it worked
plot(poly_trees[[1]])

#check names - there are n+1 individuals because added tips not changed original tip (eg t13, t13_1 etc)
poly_trees[[2]]$tip.label

######################
#2 Prepare spectra
######################
#In correct format with intraspecific variation and without (means)

#need 100 spectra for the 100 species plus variation within species

##simulate traits under different processes

#random spectra - no brownian motion on phylogeny
#one optimum spectra - OU model - specify one optima
#three optima spectra - OU model - specify three optima

#PROSPECT

#parameters that need to vary
#cab = varies from 10-400?
#nit = varies from 1-1.5
#car = varies from 10-200
#anth = 0.8-1.2
#cbrown
#cw
#cm = varies from 0.01 to 0.02

#first just use single tree
tree

#set regime (two states)
regime1<-as.vector(c(rep("Woody",60),rep("Nonwoody",40))); names(regime1)<-tree$tip.label
regime2<-as.vector(c(rep("Woody",50),rep("Nonwoody",50))); names(regime2)<-tree$tip.label

#plot states onto tree
tree<-make.simmap(tree, regime2, model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Woody","Nonwoody")

# Plot of the phylogeny for illustration
plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

## Simulate trait evolution according to a bivariate "BMM" model
# Number of traits
ntraits<-2
ntraits<-5

# Number of simulated (pairs of) traits
nsim<-1

#set lower limit (only works with fit_t_pl)
low = 1e-5  
#I need to figure out relationship between matrices and actual rates across the tree

# Rates matrices for the "Woody" and the "Nonwoody" regimes
sigma<-list(Woody=matrix(c(0.05,0.02,0.01,0.02,0.05),5, 5), Nonwoody=matrix(c(0.5,0.2,0.1,0.2,0.5),5, 5))

#Two parameter sigmas
sigma<-list(Woody=matrix(c(0.05,0.02,0.02,0.05),2), Nonwoody=matrix(c(0.1,0.2,0.2,0.1),2))

#ancestral states for each traits
theta<-c(1,100,1,100,1)

all_param_list <- mvSIM(tree, model="BMM", nsim=nsim, param=list(ntraits=ntraits,sigma=sigma,theta=theta))

#two paramter simulation

theta<-c(1,1)
two_param_list <- mvSIM(tree, model="BMM", nsim=nsim, param=list(ntraits=ntraits,sigma=sigma,theta=theta))

#evolve parameters (just focus on one)
theta = c(1,1)
anth_list <- mvSIM(tree, model="BMM", nsim=nsim, param=list(ntraits=ntraits,sigma=sigma,theta=theta))

theta<-c(100,100)
cab_list <- mvSIM(tree, model="BMM", nsim=nsim, param=list(ntraits=ntraits,sigma=sigma,theta=theta))

theta<-c(1,1)
nit_list <- mvSIM(tree, model="BMM", nsim=nsim, param=list(ntraits=ntraits,sigma=sigma,theta=theta))

theta<-c(100,100)
car_list <- mvSIM(tree, model="BMM", nsim=nsim, param=list(ntraits=ntraits,sigma=sigma,theta=theta))

theta<-c(1,1)
cm_list <- mvSIM(tree, model="BMM", nsim=nsim, param=list(ntraits=ntraits,sigma=sigma,theta=theta))

#random sample use either sample(low:high, num nums, replace = T/F) or runif(num nums, low, high)
spec_list <- c()

for (i in 1:100){
  ind_param <- all_param_list[i,]
  cab <-  ind_param[1]
  nit <- ind_param[2]
  car <- ind_param[3]
  anth <- ind_param[4]
  cm <- ind_param[5]
  spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")
  if (i == 1) {
    spec_list <- spec
  } else {
    spec_list <- c(spec_list, spec)
  }
}

for (i in 1:100){
  ind_param <- two_param_list[i,]
  cab <-  cab_list[i]
  nit <- nit_list[i]
  car <- car_list[i]
  anth <- anth_list[i]
  cm <- cm_list[i]
  spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")
  if (i == 1) {
    spec_list <- spec
  } else {
    spec_list <- c(spec_list, spec)
  }
}

for (i in 1:100){
  cab <-  cab_list[i]
  nit <- nit_list[i]
  car <- car_list[i]
  anth <- anth_list[i]
  cm <- cm_list[i]
  spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")
  if (i == 1) {
    spec_list <- spec
  } else {
    spec_list <- c(spec_list, spec)
  }
}
spec <- PROSPECT(N = nit, Cab = cab, Car = car, Anth = anth, Cm = cm, transmittance = FALSE, parameterList = NULL, version = "D")

plot(spec,ylim=c(0,0.6),xlim=c(400,2500))


# plot(spec_list[[1]],ylim=c(0,2),xlim=c(400,2500))
# points(400:2500,spec_list[[2]]@spectra@spectra_ma,type="l",col="red")
# points(400:2500,spec_list[[3]]@spectra@spectra_ma,type="l",col="blue")
# points(400:2500,spec_list[[4]]@spectra@spectra_ma,type="l",col="green")

spec_list[[1]]

class(spec_list[[1]])

spec_list[[1]]@spectra@spectra_ma

#mvMORPH

#under OU model with default parameters with two states on tree

#have tree from 
tree <- pbtree(n=100)
set.seed(14)

#this puts discrete character onto tree
state <- as.vector(c(rep("Woody", 60), rep("Herb", 40)))
names(state) <- tree$tip.label

tree <- make.simmap(tree, state, model = "ER", nsim =1)
tree

#colour/plot tree
col = c("brown", "green")
names(col) <- c("Woody", "Herb")
plotSimmap(tree, col, fsize = 0.6, node.numbers = FALSE, lwd = 3, pts= FALSE)

#then simulate trait data onto tree (MV trait)
set.seed(101)

#different parameters - need to adjust matrices to size of traits (2100) - how to do this? Otherwise, uses default values
alpha<-matrix(c(1.1,-0.9,-0.9,1.1),2) #strength of pull?
sigma<-matrix(c(0.35,0.06,0.06,0.35),2) #rate of evolution?
theta<-c(5.5,5.1,1.2,1.4) #optimum?
ntraits<-2
#this is where my data differ from tutorial - need 2100 traits not 2
ou_data<-mvSIM(tree, model="OUM", nsim=2100, param=list(ntraits=ntraits, sigma=sigma, alpha=alpha, theta=theta, names_traits=c("limb.length","limb.width")))
str(ou_data)


#make new parameters that make sense with spectral data

#two optima of 2101 traits, so concat two spectra into one theta list
length(spec_list[[1]]@spectra@spectra_ma[1,])
theta <- c(spec_list[[1]]@spectra@spectra_ma[1,],spec_list[[2]]@spectra@spectra_ma[1,])

#length of spectra is number of traits
ntraits<-2101
ntraits<-100

#the problem is that the examples have the same number of regimes and number of traits simulated
#what if there are different number of traits and regimes?

#make alpha (strength of selectiOn) the same for all traits?
#alpha <- 
#make sigma (rate of evolution) the same for all traits?
#sigma<-list(Woody=matrix(c(2,0.5,0.5,1),2), Herb=matrix(c(5,3,3,4),2))

#This isn't right
ou_sim<-mvSIM(tree, model="OUM", nsim=1, param=list(ntraits=ntraits, theta=theta))#alpha=alpha, sigma=sigma, , param=list(ntraits=ntraits, theta=theta), names_traits=c("limb.length","limb.width"), sigma=sigma, alpha=alpha, 

model_fit<-mvOU(tree,ou_sim[[1]],model="OUM")

#could figure out how to specify different parameters beside defaults but matrix size needs to match number of traits


#set up simulation
ntraits <- 2

# Number of simulated (pairs of) traits
nsim <- 2100

# Rates matrices for the "Forest" and the "Savannah" regimes
sigma<-list(Woody=matrix(c(2,0.5,0.5,1),2), Herb=matrix(c(5,3,3,4),2))

# ancestral states for each traits
theta<-c(0,0)

# Simulate - need to figure out sigma matrix that matches in dimensions to ntraits
bm_data<-mvSIM(tree,nsim=nsim, model="BMM",param=list(ntraits=ntraits,sigma=sigma,theta=theta))

model_fit<-mvBM(tree,bm_data[[1]],model="BM1")
model_fit2<-mvBM(tree,bm_data[[1]],model="BMM")

#three datasets simulated under different conditions

#spec_list - random actual spectra
#ou_data - ou model for two states (OUM), 2100 traits but not spectra
#bm_data - bm model for two states (BM1), 2100 traits but not spectra

######################
#3 Model spectral evolution
######################
#RPANDA fit_t_pl
#phyl.pca_pl for PCA

#Based on simulations from Clavel 2019


# ----- Packages dependencies
library(MASS)
library(mvMORPH)    # >=1.0.9
library(RPANDA)     # https://github.com/hmorlon/PANDA/tree/PenalizedMLPhylo
library(phylocurve) # last gitHub update
library(glassoFast) # https://github.com/JClavel/glassoFast

# ----- Parameters simulations
nsim = 1                        # Number of simulations
n <- 100                        # Number of species
ndim <- c(100)                  # Number of traits #c(2,5,10,25,31,50)

tol=NULL                        # lower limit for the regularization parameter
up = 1                          # upper bound for the parameter search (to adjust with the model)
low = 1e-5                      # lower bound for the parameter search (to adjust with the model)
scale=1                         # Height of the tree root

# ----- Models parameters
model="OU"                       # Model simulated ("BM", "OU", "EB", or "lambda")
REML=FALSE                       # REML or ML likelihood

# Parameters values for simulations
modelValue<-switch(model,
                   "EB"={ c(-log(2)/(scale/c(0.5,1,3,5)))},
                   "OU"={ c((log(2)/(scale/c(0.5,1,3,5))))},
                   "lambda"={ c(0.2,0.4,0.6,0.8)})

# Methods for phylocurves functions (used as benchmark for the ML and Pairwise likelihood approaches
if(REML){
  methodPW <- "Pairwise REML"
  methodML <- "Full REML"
}else{
  methodPW <- "Pairwise ML"
  methodML <- "Full ML"
}

theta <- c(spec_list[[1]]@spectra@spectra_ma[1,],spec_list[[2]]@spectra@spectra_ma[1,])

theta <- spec_list[[1]]@spectra@spectra_ma[1,][1:100]

# ----- Simulations
simulations <- lapply(ndim, function(p){
  lapply(modelValue, function(valueSimul){
    lapply(1:nsim, function(x){
      
      # Phylogenetic tree
      tree <- pbtree(n=n, scale=scale)
      state <- as.vector(c(rep("Woody", 60), rep("Herb", 40)))
      names(state) <- tree$tip.label
      
      tree <- make.simmap(tree, state, model = "ER", nsim =1)
      # trait
      # Covariance matrix and traits
      Sigma <- Posdef(p)
      X <- mvSIM(tree, model="OUM", nsim=1, param=list(sigma=Sigma, theta=theta))
      
      # # Nominal LOOCV (best approximation)
      # # Archetype Ridge with H&L algorithm
      # mod1 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="Variance", REML=REML, method="RidgeArch")
      # 
      # # Alternative Ridge Nominal LOOCV
      # mod2 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="Variance", REML=REML, method="RidgeAlt")
      # 
      # # Alternative Ridge (null target) Nominal LOOCV
      # mod3 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="null", REML=REML, method="RidgeAlt")
      # 
      # # LASSO nominal LOOCV
      # mod4 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="null", REML=REML, method="LASSO")
      # 
      # # Approximations to LOOCV (faster but less accurate)
      # if(p<n){
      #   # Alternative Ridge 1st Taylor approximation
      #   mod5 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="Variance", REML=REML, method="RidgeAltapprox")
      #   
      #   # LASSO approximate LOOCV
      #   mod6 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="null", REML=REML, method="LASSOapprox")
      # }else{
      #   mod5 <- mod6 <- 0
      # }
      # 
      # # Pairwise likelihood - phylocurve (Goolsby 2016)
      # mod7 <- evo.model(tree = tree, Y = X ,method =  methodPW, model=model, plot.LL.surface = FALSE, bounds=c(low,up))
      # 
      # # ML estimate (works with p<n)
      # if(p<n){
      #   mod8 <- evo.model(tree = tree, Y = X ,method =  methodML, model=model, plot.LL.surface = FALSE, bounds=c(low,up))
      # }else{
      #   mod8 <- 0
      # }
      results <- list(Sigma=Sigma, X=X, tree=tree)
      return(results)
    })
  })
})

simulations[[1]]
str(simulations)

saveRDS(simulations, "./data/simulations/first_sim.rds")

