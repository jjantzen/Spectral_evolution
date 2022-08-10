## ---- Required libraries
library(mvMORPH) # >=1.0.9
library(corpcor)
library(glassoFast) # on gitHub for now https://github.com/JClavel/glassoFast

## ---- Function to compute the Quadratic loss
# mod.par = model parameter
# tuning = regularization parameter (0 for ML)
# Y = traits dataset
# R = (population) covariance matrix
# tree = phylogenetic tree (ape object of class "phylo")
# penalty = estimator of the covariance ("RidgeArchVar","RidgeAltVar", "RidgeAltNull", "LASSO" and "ML")
# phy.mod = evolutionary model (tree transformation). "lambda" is Pagel lambda, "OU" is Ornstein-Uhlenbeck, "BM" is Brownian-motion, "EB" is Early-Burst
# REML = (TRUE or FALSE), the method used to compute the estimate of the covariance

Qloss <- function(mod.par=0, tuning=0, Y, R, tree, penalty=c("RidgeAltNull","RidgeAltVar","RidgeArchVar","LASSO","ML"), phy.mod=c("lambda","BM","EB","OU"), REML=FALSE){
  require(corpcor)
  
  # Select the penalty
  penalty <- match.arg(penalty)
  # Select the model
  phy.mod <- match.arg(phy.mod)
  
  n <- nrow(Y)
  p <- ncol(Y)
  if(REML) n <- n-1
  
  tr <- phyblTrans(tree, mod.par, phy.mod)             
  Yk <-  apply(Y,2,function(i) pic(i,tr))   
  S <- crossprod(Yk)/n                      
  
  if(any(!is.finite(S))) return(NA)
  
  switch(penalty,
         "RidgeAltNull"={
           Target <- matrix(0,p,p)
           P <- .makePenaltyFun(S,tuning,Target,"null")
         },
         "RidgeAltVar"={
           Target <- diag(1/diag(S))
           P <- .makePenaltyFun(S,tuning,Target,"Variance")
         },
         "RidgeArchVar"={
           Target <- diag(diag(S))
           P <- (1-tuning)*S + tuning*Target
         },
         "LASSO"={
           LASS <- glassoFast(S,tuning)
           P <- LASS$wi
         },
         "ML"={
           P <- S
         }
  )
  

    if(any(!is.finite(P))){
      loss=NA 
    }else{
      if(penalty!="LASSO") P <- pseudoinverse(P) 
      loss <- sum(((P %*% R - diag(ncol(P))))^2)
    }
  
  return(loss)
}



## ---- Tree branch length transformation
# phy = ape object of class phylo
# mod.par = parameter of the model
# model = evolutionary model. Evolutionary model (tree transformation). "lambda" is Pagel lambda, "OU" is Ornstein-Uhlenbeck, "BM" is Brownian-motion, "EB" is Early-Burst

phyblTrans <- function(phy, mod.par, model="lambda"){
  if(mod.par!=0 | model=="lambda"){
    switch(model,
           "lambda"={
             if(mod.par==1) return(phy)
             rootOrig <- max(branching.times(phy))
             tips <- match(c(1:Ntip(phy)), phy$edge[,2])
             phy$edge.length <- phy$edge.length * mod.par
             phy$edge.length[tips] <- phy$edge.length[tips] + (rootOrig * (1-mod.par))
             
           },
           "OU"={
             if(mod.par<=.Machine$double.eps) return(phy)
             times <- branching.times(phy)
             names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
             Tmax<-times[1]
             phy2<-phy
             
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age <- times[which(names(times) == phy$edge[i, 1])]
               t1 <- max(times) - age
               t2 <- t1+bl
               phy2$edge.length[i] <- (1/(2*mod.par))*exp(-2*mod.par * (Tmax-t2)) * (1 - exp(-2 * mod.par * t2)) - 
                 (1/(2*mod.par))*exp(-2*mod.par * (Tmax-t1)) * (1 - exp(-2 * mod.par * t1))
             }
             phy <- phy2
           },
           "EB"={
             if(abs(mod.par)<=.Machine$double.eps) return(phy)
             
             times <- branching.times(phy);
             
             names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
             
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age = times[which(names(times) == phy$edge[i, 1])]
               t1 = max(times) - age
               t2 = t1+bl
               phy$edge.length[i] = (exp(mod.par*t2)-exp(mod.par*t1))/(mod.par)
             }    
           },
           "BM"={
             phy = phy
           })
  }
  return(phy)
}


# Alternative penalty of van Wieringen & Peeters 2016 - Computational Statistics and Data Analysis
# see also Witten & Tibshirani 2009 - J. R. Statist. Proc. B
.makePenaltyFun <- function(S,lambda,target,targM){
  
  switch(targM,
         "Variance"={
           D <- (S - lambda * target)
           Alt <- D/2 + .sqM((D %*% D)/4 + lambda * diag(nrow(S)))
         },
         "unitVariance"={
           eig  <- eigen(S, symmetric = TRUE)
           Q <- eig$vectors
           d <- eig$values - lambda*target[1]
           D <- diag(sqrt(lambda + d^2/4) + d/2)
           Alt <- Q %*% D %*% t(Q)
         },
         "null"={
           eig  <- eigen(S, symmetric = TRUE)
           Q <- eig$vectors
           d <- eig$values
           D <- diag(sqrt(lambda + d^2/4) + d/2)
           Alt <- Q %*% D %*% t(Q)
         }
  )
  
  return(Alt)
}

# Matrix square root
.sqM <- function(x){
  if(!all(is.finite(x))) return(Inf)
  eig <- eigen(x)
  sqrtM <- eig$vectors %*% diag(sqrt(eig$values)) %*% solve(eig$vectors)
  return(sqrtM)
}

# Vec operator
.vec <- function(x) as.numeric(x)

# Trace operator
.tr <- function(x) sum(diag(x))

# From Uyeda et al. 2015 - Systematic Biology
Posdef <- function (n, ev = rexp(n, 1/100)) {
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

