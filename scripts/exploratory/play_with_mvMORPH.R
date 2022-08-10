#playing with mvMORPH
#library(datelife)
library(ape)
library(mvMORPH)
#load tree and spectral data

tree_sp <- read.tree("./data/tol_tree.tre")

spectra <- readRDS("./data/spectra_by_species.rds")

meta(spectra)

#need to add branch lenghts
tree_brlen <- compute.brlen(tree_sp, method = "Grafen", power = 0.5)

class(tree_brlen)
plot(tree_brlen)
spectra

names <- strsplit(tree_brlen$tip.label, "_")
for (i in 1:length(tree_brlen$tip.label)){
  tree_brlen$tip.label[i] <- paste0(names[[i]][1],"_", names[[i]][2])
}




i <- 1
#try fitting a model?
fit1 <- mvgls()


#######their data
set.seed(14)
par(mfrow=c(1,3))
tree <- pbtree(n=10)
# Set a different regime to the monophyletic group on node "12"
tree2 <- paintSubTree(tree, node=12, state="group_2", anc.state="group_1")
plot(tree2);nodelabels(,cex=0.8)
# We can set the character state to the stem branch leading to the subtree
tree3 <- paintSubTree(tree, node=12, state="group_2", anc.state="group_1", stem=TRUE)
plot(tree3)
# Finally we can also set a different regime to the branch leading to the tip "t10"
branch_1 <- match("t10",tree3$tip.label) # alternatively: which(tree$tip.label=="t10")
tree4 <- paintBranches(tree3, branch_1, state="group_1")
# set also a change to a third state along the branch
# leading to "t2" using the "stem" argument
branch_2 <- match("t2",tree3$tip.label)
tree4 <- paintSubTree(tree4, branch_2, state="group_3", stem=0.5)
plot(tree4)


set.seed(14)
# Generating a random tree with 50 species
tree<-pbtree(n=50)
# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label
# Making the simmap tree with mapped states
tree<-make.simmap(tree, sta , model="ER", nsim=1)
# Number of simulated datasets
nsim<-1
# Rates matrices for the "Forest" and the "Savannah" regimes
# Note: use lists for multiple rates (matrices or scalars)
sigma<-list(Forest=matrix(c(2,0.5,0.5,1),2), Savannah=matrix(c(5,3,3,4),2))
# ancestral states for each traits
theta<-c(0,0)
# Simulate the "BMM" model
simul_1<-mvSIM(tree, nsim=nsim, model="BMM", param=list(sigma=sigma, theta=theta))
head(simul_1)

fit <- mvBM(tree, simul_1)

simul_2 <- simulate(fit, tree=tree, nsim=100)

bootstrap <- lapply(simul_2, function(x) mvBM(tree, x, echo=F, diagnostic=F))

log_distribution <- sapply(bootstrap, logLik)
hist(log_distribution, main="Log-likelihood distribution")
abline(v=fit$LogLik, lty=2, lwd=5, col="red")
