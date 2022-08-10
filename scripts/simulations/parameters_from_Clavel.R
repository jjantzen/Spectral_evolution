## ---- Simulations for model comparisons
source("utils.r")

# Simulations have to be conducted on the model parameters as well as the number of dimensions
# Index for simulations (arguments sent by Condor)
args <- commandArgs(TRUE)
indice <- as.numeric(args[1])


# ----- Packages dependencies
library(MASS)
library(mvMORPH)    # >=1.0.9
library(RPANDA)     # https://github.com/hmorlon/PANDA/tree/PenalizedMLPhylo
library(phylocurve) # last gitHub update
library(glassoFast) # https://github.com/JClavel/glassoFast

# ----- Random number generator
set.seed(indice+1)

# ----- Parameters simulations
nsim = 100                      # Number of simulations
n <- 100                        # Number of species
ndim <- c(2,5,10,25,31,50)      # Number of traits

tol=NULL                        # lower limit for the regularization parameter
up = 1                          # upper bound for the parameter search (to adjust with the model)
low = 1e-5                      # lower bound for the parameter search (to adjust with the model)
scale=1                         # Height of the tree root

# ----- Models parameters
model="lambda"                   # Model simulated ("BM", "OU", "EB", or "lambda")
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

# ----- Simulations
simulations <- lapply(ndim, function(p){
  lapply(modelValue, function(valueSimul){
    lapply(1:nsim, function(x){
      
      # Phylogenetic tree
      tree <- pbtree(n=n, scale=scale)
      
      # trait
      # Covariance matrix and traits
      Sigma <- Posdef(p)
      X <- mvSIM(phyblTrans(tree, valueSimul, model=model), model="BM1", nsim=1, param=list(sigma=Sigma))
      
      # Nominal LOOCV (best approximation)
      # Archetype Ridge with H&L algorithm
      mod1 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="Variance", REML=REML, method="RidgeArch")
      
      # Alternative Ridge Nominal LOOCV
      mod2 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="Variance", REML=REML, method="RidgeAlt")
      
      # Alternative Ridge (null target) Nominal LOOCV
      mod3 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="null", REML=REML, method="RidgeAlt")
      
      # LASSO nominal LOOCV
      mod4 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="null", REML=REML, method="LASSO")
      
      # Approximations to LOOCV (faster but less accurate)
      if(p<n){
        # Alternative Ridge 1st Taylor approximation
        mod5 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="Variance", REML=REML, method="RidgeAltapprox")
        
        # LASSO approximate LOOCV
        mod6 <- fit_t_pl(X, tree=tree, model=model, tol=tol, up=up, low=low, targM="null", REML=REML, method="LASSOapprox")
      }else{
        mod5 <- mod6 <- 0
      }
      
      # Pairwise likelihood - phylocurve (Goolsby 2016)
      mod7 <- evo.model(tree = tree, Y = X ,method =  methodPW, model=model, plot.LL.surface = FALSE, bounds=c(low,up))
      
      # ML estimate (works with p<n)
      if(p<n){
        mod8 <- evo.model(tree = tree, Y = X ,method =  methodML, model=model, plot.LL.surface = FALSE, bounds=c(low,up))
      }else{
        mod8 <- 0
      }
      results <- list(mod1=mod1, mod2=mod2, mod3=mod3, mod4=mod4, mod5=mod5, mod6=mod6, mod7=mod7, mod8=mod8,  Sigma=Sigma, X=X, tree=tree)
      return(results)
    })
  })
})



# # save each files with a different name
save_file <- paste("Parameters_",model,"_",indice,".rdata", sep="")
save.image(save_file)
