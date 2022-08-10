#mvMORPH tutorial from Clavel and Morlon, 2019 (Phylogenetic Regressions for Multivariate Comparative Data)

#Johanna Jantzen
#January 13, 2021

##################################

#Multivariate GLS

#functions:
#mvgls - linear model fit; counterpart of gls in "nlme"
  # uses ML (method = "LL") when p<n-q; can allow for computing multivariate counterpart to R^2
  # uses PL (method = "PL-LOOCV") for any p
#manova.gls - multivariate hypotheses tests; counterpart of anova "nlme"
  # currently (1.1.1) only uses "RidgeArch" penalty for PL approach
  # options:
    # test 
      #Pillai: performs better when heterogeneity of covariance matrices in factorial designs
      #Wilks: connection with likelihood ratio tests (LRTs) and approximations of F-statistic
      #Lawley-Hotteling:
      #Roy: more powerful when mean vectors are collinear, less powerful when mean vectors are more diffuse
    # type (I, II, III)
      #Type I (Sequential) - additive, the same for all possible ordering of explanatory variables; all variation in model accounted for
          #But importance of factors depend on order in which they are entered in model - may obtain different results with different ordering of predictors
          #Useful for assessing effect of adjusting for given variable
      #Type II (Adjusted) - do not depend on ordering of predictors; obey the principle of marginality - respect the hierarchy of predictors level (main effects and interactions)
          #Significance of tests for main effects does not account for higher-order effects terms (interactions)
          #Tests for higher-orders account for all lower-orders terms
      #Type III (Adjusted) - similar to Type II but include higher-level terms (interactions) when making adjustments
          #Tests for an effect are performed after adjusting for all other model terms - violate the marginality principle
          #Controversial test - type II tests are more powerful when interaction is nonsignificant and interpreting main effects in presence of significant interactions is of little interest
          #Do not depend on the differences in sample sizes for factorial designs (?contrasts)
      #When multiple predictors but not interaction terms - need to choose between sequential and adjusted tests
          #If there are low covariations between predictors, adjusted/non-adjusted won't make a big difference
      #When multiple regression model includes interaction terms (higher order terms), look at results of highest-order term first
          #All types of tests have same result for highest-order term - if significant, indicate we can reject null hypothesis that main effects are additive because they are interdependent
          #Lower order terms will vary between tests but if significant interactions (higher order), less informative/little interest b/c significant interaction
    # nperm (default 999)
      #Permutations to approximate the distribution of the mv statistic under null hypothesis
      #Plot permuted distribution using plot() on manova.gls output
      #Can parallelize to increase speed for higher nperm
    # nbcores (default 1); only works on Mac and Linux
    # permutation (approx or full)
      #If not using BM model (eg OU, lambda, or EB), can compute stats differently
      #Approximate strategy fixes evolutionary model to parameters estimated by mvgls function and only optimize level of regularization/penalization on each permuted dataset
      #Exact strategy ("full") optimizes both evolutionary mdoel and regularization/penalization parameters on each permuted dataset
      #No noticeable differences in simulations but approx faster and easier
    # L (contrasts coding matrix, default NULL)
    # rhs (right hand side equation of GL hypothesis test; default = 0)

##################################

#First example: one-way MANOVA with model fit and hypothesis testing

# load the package and simulate some fake data
library(mvMORPH)
set.seed(1, sample.kind = "Rounding") # for reproducibility with older R versions
n = 36; p=5
tree = pbtree(n=n)
R = crossprod(matrix(runif(p*p),p,p)) # random covariance matrix
theta = rep(0,p) # the mean of the p traits
trait <- mvSIM(tree, model="BM1", param=list(sigma=R, theta=theta)) # phenotypic traits
grp <- rep(c("forest","open"), each=n/2)
pheno <- trait + rep(c(0,10), each=n/2) # introduce some differences (~10) on "open"
size <- rTraitCont(tree) # an hypothetical "body-size" for the species

#Need trait dataset (eg spectra) and grouping categories (other variables eg habit)
#Grouping variable needs to be as a factor for MANOVA design

#First create a list object with the required data
data = list(trait=pheno, habitat=as.factor(grp), size=size)

# Fit the multivariate linear model
fit <- mvgls(trait ~ habitat, data=data, tree=tree, model="lambda") # here, model is Pagel's lambda for phylo structure, default parameters eg RidgeArch and target matrix proportional to identity

# Print the model fit and other functions for retrieving information
print(fit)
summary(fit)
residuals(fit)
vcov(fit)
GIC(fit)
coef(fit)
fitted(fit)
fit$coefficients

#Key output: lambda value (1) is equal to expected under BM
#Key output: regularization parameter small (0.0056) so no real need for strong regularization of covariance matrix here
#Key output: Coefficients (first one (ie intercept) is first factor level of grouping variable), can by eye determine differences between two groups (rows)

# The overall MANOVA test:
aov <- manova.gls(fit, nperm=999, test="Wilks", verbose=TRUE) # testing null model of no difference between two groups

# display the test statistic
plot(aov)

# Fit with parallel computing on 4 cores for 9999 permutations
aov2 <- manova.gls(fit, nperm=9999, nbcores=4L, test="Wilks", verbose=TRUE) # higher number of permutations so smoother distribution of null

# display the test statistic
plot(aov2)

#Key output: p value - determines significance (if < 0.5 etc reject null)
#Key output: estimated (Wilks) statistic, can compare with null distribution using plot
#Higher nperm, more smooth distribution of null


##################################

#Second example: multivariate multiple regression tests (multiple predictors)

#Fit the MANCOVA model
fit2 <- mvgls(trait ~ size + habitat, data=data, tree=tree, model="lambda") # with second predictor (continuous covariate) = Analysis of Covariance

fit3 <- mvgls(trait ~ habitat + size, data=data, tree=tree, model="lambda") # with second predictor (continuous covariate) = Analysis of Covariance

# The overall MANOVA test: (type I = sequential = tests effect of each predictor in order after accounting for effect of previous ones)
aov3 <- manova.gls(fit2, nperm=999, test="Wilks", verbose=TRUE) # size is not significant because was simulated independent of phylogeny/trait

aov3b <- manova.gls(fit3, nperm=999, test="Wilks", verbose=TRUE) # size is not significant because was simulated independent of phylogeny/trait

# The overall (type II) MANOVA test: (type II = adjusted = marginal contribution of each factor)
aov4 <- manova.gls(fit2, nperm=999, test="Wilks", type="II", verbose=TRUE)

aov4b <- manova.gls(fit3, nperm=999, test="Wilks", type="II", verbose=TRUE) # you get the same result regardless of order of variables in model

# Fit the MANCOVA model with interaction term
fit4 <- mvgls(trait ~ size + habitat + size*habitat, data=data, tree=tree, model="lambda") # interaction term using *; again, order shouldn't matter since type II test later

#The overall (type II) MANOVA test:
aov5 <- manova.gls(fit4, nperm=999, test="Wilks", type="II", verbose=TRUE) # no effect of interaction because simulated independently


##################################

#Third example: General Linear Hypothesis Testing (L argument)
#In cases where testing k linear combinations of parameters B (beta)

grp2 = rep(1:3, each=n/3)
pheno2 = trait + rep(c(0,0,10), each=n/3) # introduce some differences on the third group
data2 = list(pheno=pheno2, grp=as.factor(grp2))

# fit the model
fit5 <- mvgls(pheno~grp, data=data2, tree=tree, model="lambda")

# The overall MANOVA test:
aov6 <- manova.gls(fit5, verbose=TRUE) # significant overall because one of groups different, but not tested which here

#To test more specific hypotheses, include L (linear combinations)
#Where difference between the two options is zero, and the sum of all factors is zero???

L1 <- matrix(c(0,1,-1), ncol=3) # Contrasts (vector or matrix) - this one is asking second and third are the same?; first one is not involved in hypothesis so set to zero

# Test the first contrast:
aov6b <- manova.gls(fit5, L=L1, nperm=999, verbose=TRUE) # significant, so b2 and b3 are not the same

#different hypothesis - where testing sameness of b1 and b2
L2 <- matrix(c(0,1,0), ncol=3) # Contrasts (vector or matrix)

# Test the second contrast:
aov6c <- manova.gls(fit5, L=L2, nperm=999, verbose=TRUE) # not significant, so the first and second are the same

#implied that because 1=2 and 2 != 3, 1 != 3

#significance means difference here because testing deviation from sameness

# fit the model without "intercept" - different approach than the above, where testing all groups explicitly and not compared to baseline of first factor
fit5bis <- mvgls(pheno ~ grp + 0, data=data2, tree=tree, model="lambda")

#new linear (intuitive??) relationship testing for difference in one and two?
L2bis <- matrix(c(1,-1,0), ncol=3) # Contrasts (vector or matrix)

# Test the second contrast:
aov6d <- manova.gls(fit5bis, L=L2bis, nperm=999, verbose=TRUE) # should be same (p value may vary) as above second contrast

#Be careful about hypothesis that is being tested especially when using non intercept models

##################################

#Fourth example: Testing hypothesis on values of estimated parameters (rhs argument)
#For testing specific values for parameters or differences (rather than just different/not)

L3 <- matrix(c(0,0,1), ncol=3) # define coding - this is the most confusing part, how to define the L matrix
rhs <- 2 # is the difference between the first and the third group equal to 2? (in numerical terms)

# Test the first contrast:
aov7 <- manova.gls(fit5, L=L3, rhs=rhs, nperm=999, verbose=TRUE) # significant because different from 2

# indeed it is rejected (recall that beta3=10)... try with rhs=10 instead:
aov7b <- manova.gls(fit5, L=L3, rhs=10, nperm=999, verbose=TRUE) # not different from 10, as expected, because it was 10

#I don't think I'd want to do this, but can test rhs without L matrix?
aov7c <- manova.gls(fit5, rhs=10, nperm=999, verbose=TRUE) # testing without the L matrix, assumes matrix contains only that value - not sure how this works


##################################

#Fifth example: Bat example from paper (Monteiro and Nogueira 2001 data)

#Load the data
data <- get(data(phyllostomid))

# head(data$mandible) # the dataset is high-dimensional - i.e. p>(n-q)
#n = 6, p = 36*2?

str(data)
#list of different types of data, each as lists of factors or numbers
#first, mandible data which is the trait data (ie spectral data)
#next, grouping variables (diet, class of eating), eg different categories (habit, biome, clade etc)
#last, tree, with tip labels matching names in data lists

# Fit the mandible data to the four diet and feeding modes schemes
fit_grp1 <- mvgls(mandible~grp1, data=data, data$tree, model="lambda", method="PL-LOOCV")
fit_grp2 <- mvgls(mandible~grp2, data=data, data$tree, model="lambda", method="PL-LOOCV")
fit_grp3 <- mvgls(mandible~grp3, data=data, data$tree, model="lambda", method="PL-LOOCV")
fit_grp4 <- mvgls(mandible~grp4, data=data, data$tree, model="lambda", method="PL-LOOCV")

# Then we perform the MANOVAs as in the paper
nbcores=4L # to speed up the calculations, but doesn't work on my computer

man1 <- manova.gls(fit_grp1, test="Wilks", nbcores=nbcores, verbose=TRUE) # Grouping 1
man2 <- manova.gls(fit_grp2, test="Wilks", nbcores=nbcores, verbose=TRUE) # Grouping 2
man3 <- manova.gls(fit_grp3, test="Wilks", nbcores=nbcores, verbose=TRUE) # Grouping 3
man4 <- manova.gls(fit_grp4, test="Wilks", nbcores=nbcores, verbose=TRUE) # Grouping 4

plot(man1)

# First define the set of contrasts for diet (factors considered in alphabetical order); here use intercept to allow estimation for each group instead of from a baseline
L <- rbind(c(0,0,0,1,-1), # Nectarivory vs Sanguivory
           c(0.5,-1,0.5,0,0), # Frugivory vs. Animalivory
           c(1,0,-1,0,0)) # Insectivory vs. Carnivory

#fit the model with the intercept (+0) - for grp4 which is the diet grouping with 5 levels
fit_glh <- mvgls(mandible~grp4+0, data=data, data$tree, model="lambda", method="PL-LOOCV")

#then conduct three separate tests for each row in the L matrix - tests for specific levels within grouping
L1 <- L[1,,drop=FALSE] # Nectarivory vs Sanguivory
manL1 <- manova.gls(fit_glh, test="Wilks", L=L1, verbose=TRUE)

L2 <- L[2,,drop=FALSE] # Frugivory vs. Animalivory
manL2 <- manova.gls(fit_glh, test="Wilks", L=L2, verbose=TRUE)

L3 <- L[3,,drop=FALSE] # Insectivory vs. Carnivory
manL3 <- manova.gls(fit_glh, test="Wilks", L=L3, verbose=TRUE)

plot(manL3)

#other functions to get more information from fitted model

#add more information to these lines (what are they doing?)
fit_grp1 <- mvgls(mandible~grp1, data=data, data$tree, model="lambda", method="PL-LOOCV")
fit_grpOU <- mvgls(mandible~grp1, data=data, data$tree, model="OU", method="PL-LOOCV")
result <- mvOU(tree, data$trait)

coef(fit_grp1)
fitted(fit_grp1)
halflife(result)

mv.Precalc(tree, nb.traits = 34)

signal <- mvgls(mandible~1, data=data, data$tree, model="lambda", penalty="RidgeArch")
summary(signal)
residuals(signal)
vcov(signal)
summary(mvgls(mandible~1, data=data, data$tree, model="BM", penalty="RidgeArch", method="H&L"))

#me testing stuff
ou_fit <- mvgls(mandible~1, data=data, data$tree, model="OU", penalty="RidgeArch", method="H&L")
ou_fit2 <- mvgls(mandible~grp1, data=data, data$tree, model="OU", penalty="RidgeArch", method="H&L")
residuals(ou_fit)
ou_fit$param
summary(ou_fit2)
unique(data$grp1)
ou_manova <- manova.gls(ou_fit, test="Wilks", verbose=TRUE)
ou_manova$pvalue
data$grp1

ggplot()+
  #geom_histogram(aes(x = data$grp1, y = data$mandible[,1]))
  #geom_histogram(aes(x = data$mandible[,1], colour = data$grp1))
  geom_boxplot(aes(x = data$grp1, y = data$mandible[,3]))


head(data)
#try to recreate their pca figure
fit_pca <- mvgls(Y~1, tree=tree, model="lambda", method="LOO", penalty="RidgeAlt")
bats_pca <- mvgls(mandible~1, data=data, data$tree, model="lambda", method="PL-LOOCV") # not sure what the ~1 means but it gives the same results as theirs
mvgls.pca(bats_pca, plot=TRUE, axes = c(1,3))

##################################

#Web example: Phylogenetic PCA

set.seed(1)
n <- 32 # number of species
p <- 30 # number of traits

tree <- pbtree(n=n) # phylogenetic tree
R <- crossprod(matrix(runif(p*p),p))  # a random symmetric matrix (covariance)

# simulate a dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R))

#using other exampe of mvgls output to do pca
data$grp3
morph_fit_lamb <- mvgls(data$mandible ~ 1, tree=data$tree, model="lambda", method = "PL-LOOCV")
morph_reg_lamb <- mvgls(data$mandible ~ data$grp3, tree=data$tree, model="lambda", method = "PL-LOOCV") 
panda_fit_lamb <- fit_t_pl(data$mandible, data$tree, model = "lambda", method = "RidgeAlt", SE = TRUE)

# The conventional phylogenetic PCA
phylo_pca <- mvgls(Y~1, tree=tree, model="BM", method="LL") # true likelihood for smaller datasets?
mvgls.pca(phylo_pca, plot=TRUE) 

pca_lamb <- mvgls.pca(morph_reg_lamb, plot=TRUE) #for regression model
mvgls.pca(morph_fit_lamb, plot=TRUE) #for regression model to intercept (~1) - no difference from regular regression
mvgls.pca(panda_fit_lamb, plot=TRUE) #for panda model of evolution of mv trait - doesn't work bc wrong class of object

head(pca_lamb$scores) #but how to interpret this? Maybe do regression of primary axes instead???

# fit a multivariate Pagel lambda model with Penalized likelihood
fit_pca <- mvgls(Y~1, tree=tree, model="lambda", method="LOO", penalty="RidgeAlt") # I will do this approach because of high dimensionality

# Perform a regularized phylogenetic PCA using the model fit (Pagel lambda model)
pca_results <- mvgls.pca(fit_pca, plot=TRUE) 

# retrieve the scores
head(pca_results$scores)

#try doing regression using primary PC axes
fit_pca_reg <- mvgls(pca_lamb$scores ~ data$grp3, tree = data$tree, model="lambda", method = "PL-LOOCV")
fit_pca_reg_reduced <- mvgls(pca_lamb$scores[,1:3] ~ data$grp3, tree = data$tree, model="lambda", method = "PL-LOOCV")

#compare to when using all data
morph_reg_lamb #trad version
fit_pca_reg #pca version with all axes (same as trad)
fit_pca_reg_reduced #pca version with fewer axes (different from others)

##################################

#Assess model fit using GIC for models computed with penalized likelihood - want the biggest negative number
#Testing GIC for model comparison

#Making up example (multiple predictors not been fully tested so proceed with caution with multiple)
gic_fit_lam <- mvgls(mandible~1, data=data, data$tree, model="lambda", method="PL-LOOCV")
lamb_reg <- mvgls(mandible~grp2, data=data, data$tree, model="lambda", method="PL-LOOCV")

gic_fit_bm <- mvgls(mandible~1, data=data, data$tree, model="BM", method="PL-LOOCV")
gic_fit_ou <- mvgls(mandible~1, data=data, data$tree, model="OU", method="PL-LOOCV")
gic_fit_eb <- mvgls(mandible~1, data=data, data$tree, model="EB", method="PL-LOOCV")
GIC(gic_fit_lam)
GIC(gic_fit_bm)
GIC(gic_fit_ou)
GIC(gic_fit_eb)

gic_fit_ou$opt
gic_fit_eb$param
gic_fit_bm$model


gic_fit_lam$param
gic_fit_eb$param

#AIC and LRT for those models under ML (not PL)

##################################

#Modeling evolution outside of regression (not mvgls)

#For non high-dimensional data, you can use mvBM, mvEB, and mvOU to model evolution (and with multiple rates/optima of selective regimes mapped onto tree)
#But the PL method (method = "PL-LOOCV") doesn't seem to work with those functions
#How to do the equivalent of mvOU with multiple optima/selective regimes with high dimensional data?
#Can do mvgls with (trait ~ 1) but this doesn't give the same results (values or parameters) as the mvOU and mvBM functions
#Regressions (x ~ y) versus just fitting the data (X~1) (against the intercept?)

#Just fitting BM model to MV data - allows for comparison of different models

# try it with the original tutorial data
set.seed(1, sample.kind = "Rounding") # for reproducibility with older R versions
n = 36; p=5
tree = pbtree(n=n)
R = crossprod(matrix(runif(p*p),p,p)) # random covariance matrix
theta = rep(0,p) # the mean of the p traits
trait <- mvSIM(tree, model="BM1", param=list(sigma=R, theta=theta)) # phenotypic traits
grp <- rep(c("forest","open"), each=n/2)
pheno <- trait + rep(c(0,10), each=n/2) # introduce some differences (~10) on "open"
size <- rTraitCont(tree) # an hypothetical "body-size" for the species
data = list(trait=pheno, habitat=as.factor(grp), size=size)

# make models
bm_1 <- mvBM(tree, data$trait, model = "BM1", method = "pic")
bm_2 <- mvBM(tree, data$trait, model = "BM1", method = "pic", param=list(constraint="diagonal"))
bm_3 <- mvBM(tree, data$trait, model = "BM1", method = "pic", param=list(constraint="equal"))

# compare models
results <- list(bm_1,bm_2,bm_3)
aicw(results, aicc=TRUE)

# test simulating ancestral states on tree

set.seed(14)
# Generating a random tree with 50 species
tree<-pbtree(n=50)
# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label
# Making the simmap tree with mapped states
tree<-make.simmap(tree, sta , model="ER", nsim=1)

# Rates matrices for the "Forest" and the "Savannah" regimes
# Note: use lists for multiple rates (matrices or scalars)
sigma<-list(Forest=matrix(c(2,0.5,0.5,1),2), Savannah=matrix(c(5,3,3,4),2))
# ancestral states for each traits
theta<-c(0,0)

# Simulate the "BMM" model
simul_1<-mvSIM(tree, nsim=1, model="BMM", param=list(sigma=sigma, theta=theta))

mvBM_sim <- mvBM(tree, simul_1) # it's not very close but kind of approximate to what was simulated

#These functions likely won't work with my data because they're too highly dimensional - but I could try with reduced PCA axes???

##################################

#Comparing mvgls and mvOU and mvBM functions with higher dimensional data

#From "How to use mvMORPH" vignette

#Testing multiple means (BM)
set.seed(1)
tree<-pbtree(n=50)

# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label #categorical variable for different states to map

# Making the simmap tree with mapped states
tree<-make.simmap(tree,sta , model="ER", nsim=1) #need to map states onto tree for multiple rates etc

# Simulate the traits
sigma<-matrix(c(0.1,0.05,0.05,0.1),2)
theta<-c(0,2,0,2)
data<-mvSIM(tree, param=list(sigma=sigma, theta=theta, smean=FALSE), model="BM1", nsim=1)

# Model simulated data (using mvgls and mvBM respectively) - these are very inconsistent, why?
fit_1 <- mvgls(data ~ 1, tree=tree, model="OU") 
fit_2 <- mvgls(data ~ 1, tree=tree, model="BM")

fit_3 <- mvBM(tree, data, model="BM1")
fit_4 <- mvBM(tree, data, model="BMM")

summary(fit_2)
fit_2$param

#Simulate OU data
set.seed(100)
tree<-rtree(100)

# Simulate the traits
alpha<-matrix(c(0.2,0.05,0.05,0.1),2)
sigma<-matrix(c(0.1,0.05,0.05,0.1),2)
theta<-c(0,2,0,1.3)
data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha, 
                             theta=theta, root=TRUE), model="OU1", nsim=1)


fit_5 <- mvgls(data ~ 1, tree=tree, model="OU") 
fit_6 <- mvOU(tree, data, model="OU1", method = "PL-LOOCV", param=list(root=TRUE))
tree

#I may need to see if I can get the code where they did the models for the 2019 paper from the authors
#They seem to do a similar thing to the mvOU models but with highly dimensional data

#They used a different package (RPANDA) for the 2019 paper
library(RPANDA)

#load the bat data
data <- get(data(phyllostomid))

#I should get what they did? Except they used different data (not availble)
#Account for intraspecific variation - option in mvgls function with "error = TRUE"
#Think about type of penalty (ridgearch vs lasso etc)
#I still think this is more limited than the mvOU method because of the inability to specify the different selective regimes by mapping another trait?
panda_fit <- fit_t_pl(data$mandible, data$tree, model = "OU", method = "RidgeAlt", SE = TRUE)
panda_fit_bm <- fit_t_pl(data$mandible, data$tree, model = "BM", method = "RidgeAlt", SE = TRUE)
panda_fit_eb <- fit_t_pl(data$mandible, data$tree, model = "EB", method = "RidgeAlt", SE = TRUE)
panda_fit_lamb <- fit_t_pl(data$mandible, data$tree, model = "lambda", method = "RidgeAlt", SE = TRUE)

#These are just models of the spectra (not regressions) - appropriate for testing best fit model of evolution to high dimensional traits

#This works to calculate the actual phylosignal of the data but only for a single trait (not mv)
phylosig_lamba <- phylosig(data$tree, data$mandible, method = "lambda", test = TRUE)
phylosig_lamba

#Compare GIC
results <- list(GIC(panda_fit_bm),GIC(panda_fit_lamb),GIC(panda_fit), GIC(panda_fit_eb))
results

panda_fit$model.par
panda_fit$model
panda_fit_lamb$model.par
panda_fit_bm$model.par
panda_fit_eb$model.par

#This is a regression against itself (intercept) which gives different results (why?)

morph_fit_lamb <- mvgls(data$mandible ~ 1, tree=data$tree, model="lambda", method = "PL-LOOCV") 
morph_fit_OU <- mvgls(data$mandible ~ 1, tree=data$tree, model="OU", method = "PL-LOOCV") 
morph_fit_BM <- mvgls(data$mandible ~ 1, tree=data$tree, model="BM", method = "PL-LOOCV") 
morph_fit_EB <- mvgls(data$mandible ~ 1, tree=data$tree, model="EB", method = "PL-LOOCV") 

GIC(morph_fit_lamb)
GIC(morph_fit_OU)
GIC(morph_fit_BM)
GIC(morph_fit_EB)

#Or try Bootstrap Information Criterion, or EIC
eb_eic <- EIC(morph_fit_EB)
ou_eic <- EIC(morph_fit_OU)
lamb_eic <- EIC(morph_fit_lamb)
bm_eic <- EIC(morph_fit_BM)


#How to go from PC back to trait itself? Like in Clavel et al 2019
#Reconstruct ancestral states using ancestral?
recon_trait <- ancestral(panda_fit_lamb)

##################################

#Thinking about normality of data and/or residuals: http://blog.phytools.org/2013/02/a-comment-on-distribution-of-residuals.html
#We do not expect the residuals or data to be normal, but we do expect the residual error CONTROLLING FOR BOTH MAIN EFFECTS IN MODEL AND FOR PHYLOGENY to be normal

# is x2 normal? (should fail)
lillie.test(x2)

#this is the model
fit<-gls(y~x1+x2,data=data.frame(x1,x2,y), correlation=corBrownian(1,tree))

# are the residuals normal? (should fail)
lillie.test(residuals(fit))

# are the residuals controlling for phylogeny normal?
# (should pass)
lillie.test(chol(solve(vcv(tree)))%*%residuals(fit))


##################################

#Assumptions about the tree

#Ultrametric

#Looking for phylo signal in residuals of regression model

###############
#From Book
#Tree uncertainty, check out De Villemereuil et al 2012 for regressions over range of trees
#Account for phylogenetic uncertainty

#Details of estimation process can have impacts on reliability of results (Revell 2010)
  #Heterogeneities of the underlying process across the clade investigated (Garland and Ives 2000)
  #Errors in phylogenetic tree (Diaz-Uriarte and Garland 1998)

#Assumptions about distribution of residuals
#Distribution of predictor(s)
#Number of predictors compared to sample size
#Absence of collinearity and influential cases

#Questions of model validity, stability, and reliability

#Assumption of correct and complete data and predictors
#Missing data occur at random (phylogenetically speaking)

#Is the model appropriate for the question? (conceptually)
  #Which predictors? - needs to be reasoned not stats
  #make decisions about interactions a priori too (disregarding data but with reasoning)
  #give priority to the right model (including right terms) over size of data limitations

#Technical validity of the model
  #two main effects of two-way interaction must be in model (considering interactions)
  #three-way interaction: must include three two-way interactions and three main effects
  #squared covariate: must include unsquared covariate

#Scaling
  #z-transform = mean of zero and SD of 1 (after other transformations which may be needed)
  #not needed but can facilitate interpretation because coefficients represent average change in response per Standard deviation(SD)
    #therefore directly comparable between covariates
  #can also enhance model interpretability when interactions included
    #without transforming: interaction between two covariates: coefficients derived for main effects indicate their effect at the respective other covariate having a value of zero
    #with z transformation: coefficients indicate effect of two covariates at average of respective other covariate (more reasonable)
  #create plots easier (ignore all other terms in model when plotting implies assuming them to be at their average)
  #disadvantage: coefficients for different data not comparable because SD will vary between datasets
    #so report original SD of covariates

###Prior to fitting model
 
#sample size N and number of predictors k
  #N should be >>> k normally
  #power decreases otherwise
  #sample size = think about expected effect sizes ---- in phylo analyses, N should be ~ and > 10*k
  #otherwise, do power analysis based on simulations of expected/minimum effect size using phylogenetic data close to ones analyzed
  #how to exclude predictors if needed:
    #based on reasoning of which would be least likely of relevance for response under question
    #do PC or factor anaysis and use derived scores rather than original covariates
    #exclude predictors based on collinearity
  #for phylo analyses: limitation in taxa, so better to do right model than oversimplified to meet assumed limit
  #model too complex may appear unstable

#Distribution of quantitative predictors/covariates
  #No assumptions about distribution of covariates
  #But inspect distribution using histogram/qq-plot
  #See influential cases
  #Want roughly symmetrical distribution (normal or uniform) to avoid problematic cases
    #can use log or square root transformation

#Categorical predictors (factors)
  #Dummy code them - then modeled as normal predictors
  #But inspect frequency of levels - don't want any too rare - would be unduly influential
  #May be better to drop than to add additional factors into model
  #with multiple factors, check number of times combinations of levels occur
  #with interactions, only include when each combo of levels occurs at least twice in data

#Collinearity
  #absence of strong collinearity among predictors important for validity of results
  #ie not redundant or not providing new information
  #check for collinearity: VIF (variance inflation factors)
    #When R^2 too large, VIF gets large; want a value of 1 for VIF
    #no set threshold for what's too big
  #if have collinearity, can omit, or can combine predictors using PC eg
  #sometimes unavoidable - if included to account for predictor, fine to interpret other key test predictors
          # - assess model stability to determine if collinearity affects conclusion (eg with and without collinear variables)


###After fitting model

#Distribution of residuals
  #Normality of residuals - critical assumption
    #inspect distribution using qq-plot (or histogram)
    #if skewed, consider transformation of predictor and/or response
    #if outliers, are there predictors missing? or evolutionary singularity?
  #Homogeneity of residuals (homoskedasticity)
    #variation in residuals should be the same regardless of the constellation of values of predictors
    #inspect by plotting residuals of model (y) agasint fitted values (x) - want to see no pattern
      #same scatter around zero across entire range of fitted values
      #can transform it if there is heteroskedasticity
    #can be from misspecified model eg non-linear effect of covariate, impact of predictor varies across phylogeny, or main effect or interaction missing
    #if make changes after seeing this, need to account for a posteriori change when inferences made

#Absence of influential cases
  #essentially model stability - tested by removing cases and evaluating impact
  #no simple answer if there are influential cases
  #how to do this in R?

###Drawing conclusions

#Full-Null Model Comparison
  #to protect from false positives
  #null model is model comprising only the intercept (~1); compare with full model and see p value for impact of predictors
  #predictors are considered significant if overall model significant?
  #test predictor vs control predictor - eg those important for hypothesis vs those just for controlling for effects

#Inference about individual terms
  #Factors with more than two levels: dummy variable for each level of the factor except reference level
    #p values given for arbitrary comparison to reference level
    #do full-null comparison with full model and model without factor (F-test); if test is significant, then factor significant impact on response
  #Conditional interpretation of interaction dependent on other main effect it is interacting with

#Decisions to be made
  #Which main effects, interactions, nonlinear terms to include
  #potential transformations of predictors/response and subsequent z-transformations of covariates
  #checks of model validity and stability
  #report these all in paper
  #outline reasoning when formulating model, clearly formulate model analysed, describe preparatory steps, describe evaluation steps (assumptions, stability)
  #

  

#############
#From chapter 15 on OU models

#uncertainty in parameters - parametric bootstrapping
  #assumption that model is true, what distribution of parameter estimates recovered if evolution were rerun under that model?
#low and high alpha (strength of selection) - impacts ability to infer processes
  #high alpha - erases past processes
  #low alpha - final trait values depend less on theta (optima)
#state of root can be problematic
  #if little info of past due to strong attraction, biased towards zero (value)
  #even separately estimating root can cause errors - little info about root state
#OU doesn't always approach BM if alpha is zero (depends on software)

#meaning of parameters
  #theta - units of trait under investigation
  #sigma squared - units of trait units squared over branch length units
  #alpha - units of reciprocal branch length units (eg my^-1)
  #phylogenetic halflife - time units
  #need to report these units
  #are trees rescaled? - rescale or rethink re interpretation

#think about complexity of model relative to number of taxa
#think about using different starting points for parameter optimization
  #or rerun exact same analyses and compare results (inconsistency in parameter estmiates - issue)
#if nested models (restriction of other) - more restricted model likelihood must be less than general model (but not AIC score)
#run same analysis in two programs (eg mvMORPH and RPANDA)
#check if parameters at min or max of preset bounds
#consider using model averaging if similar support for different models
#map regimes (ie discrete characters eg herbaceous/woody) onto phylogeny - then model evolution of other trait (eg spectra)
#see this chapter for example with woody/herbaceous plants



##############
##Chapter 7 - Uncertainties Due to Within-Species Variation in Comparative Studies: Measurement Errors adn Statistical Weights
#Laszlo Zsolt Garamszegi

#Sources of variation ie measurement error s.l.
  #Instrument-related errors
  #Observer errors - did different people consistently sample different taxa?
  #Individual physiological or behavioural differences - over time/contexts
  #Sex or age-specific effects
  #Population differences due to phylogeography/gene flow
#Methods to deal with error deal with it in the aggregate (hard to separate out the difrerent components)

#Statistical treatment of correlating measurement errors - see Ives et al 2007 and Hansen and Bartoszek 2012

#Consequences
  #If random measurement error (eg measurement itself) and univariate - loss of power and failure to reject false null hypothesis
  #If two or more variables of interest - can affect parameter estimtes too (as well as precision/significance)
    #Assume at least predictors have been measured without error (regression)
    #If variation within species (subjects), will bias toward zero (underestimate true parameters)
      #Underestimation of R^2 and standardized regression coefficients (if variation in response variable) AND
      #Underestimation of unstandardized regression estimates (if variation in predictor variable)
      #Even more complex in non-linear regressions
    #Attenuation bias - downward bias in parameter estimates by within-species variability
    
#Covariances between multiple traits - between-subject and within-subject
  #


#Testing phylosignal aspect of fitted models
library(phytools)
library(mvMORPH)

set.seed(1)
n <- 32 # number of species
p <- 50 # number of traits (p>n)

tree <- pbtree(n=n, scale=1) # phylogenetic tree
R <- crossprod(matrix(runif(p*p), ncol=p)) # a random covariance matrix
# simulate a BM dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R, theta=rep(0,p))) 
data=list(Y=Y)

signal <- mvgls(Y~1, data=data, tree, model="lambda", penalty="RidgeArch")
summary(signal)


fit1 <- mvgls(Y~1, data=data, tree, model="BM", penalty="RidgeArch")
fit2 <- mvgls(Y~1, data=data, tree, model="OU", penalty="RidgeArch")
fit3 <- mvgls(Y~1, data=data, tree, model="EB", penalty="RidgeArch")

GIC(fit1); GIC(fit2); GIC(fit3) # BM have the lowest GIC value

summary(fit1)

library(geomorph)
