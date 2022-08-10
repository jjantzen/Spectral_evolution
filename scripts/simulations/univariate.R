
#Simulate univariate trait under BM
X<-fastBM(anoletree,nsim=10)

#multiple rates - need variance covariance matrix
X2 <- sim.corrs(tree, vcv, anc=NULL, internal=FALSE)

