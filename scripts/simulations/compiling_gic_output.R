#compiling output from cluster (gic comparisons)
library(dplyr)

#BM data

#load
gic_comparison_bm_df_bm <- readRDS("./analysis/from_cluster/full/gic_comparison_bm_df_bm.rds")
gic_comparison_bm_df_lam <- readRDS("./analysis/from_cluster/full/gic_comparison_bm_df_lam.rds")
gic_comparison_bm_df_ou1 <- readRDS("./analysis/from_cluster/full/ou_gic_comparison_bm_df_1.rds")
gic_comparison_bm_df_ou2 <- readRDS("./analysis/from_cluster/full/ou_gic_comparison_bm_df_2.rds")

#combine
bm_gics <- gic_comparison_bm_df_bm

bm_gics$LAM <- gic_comparison_bm_df_lam$LAM
bm_gics$OU <- gic_comparison_bm_df_ou1$OU
bm_gics$OU[51:100] <- gic_comparison_bm_df_ou2$OU[51:100]

#get smalles number - put column name into new column
count_reps_bm <- bm_gics %>% 
  mutate(best_model = apply(.[,2:5], 1, function(x) names(x)[which.min(x)])) %>% 
  #dplyr::select(best_model) %>% 
  group_by(best_model) %>% 
  summarize(best_model = best_model, count = n()) %>% 
  unique()


#OU data

#load
gic_comparison_ou_df_bm <- readRDS("./analysis/from_cluster/full/gic_comparison_ou_df_bm.rds")
gic_comparison_ou_df_lam <- readRDS("./analysis/from_cluster/full/gic_comparison_ou_df_lam.rds")
gic_comparison_ou_df_ou <- readRDS("./analysis/from_cluster/full/gic_comparison_ou_df_ou.rds")
gic_comparison_ou_df_eb <- readRDS("./analysis/from_cluster/full/gic_comparison_ou_df_eb.rds")

#combine
ou_gics <- gic_comparison_ou_df_bm
ou_gics$EB <- gic_comparison_ou_df_eb$EB
ou_gics$LAM <- gic_comparison_ou_df_lam$LAM
ou_gics$OU <- gic_comparison_ou_df_ou$OU

ou_gics

#get smalles number - put column name into new column
count_reps_ou <- ou_gics %>% 
  mutate(best_model = apply(.[,2:5], 1, function(x) names(x)[which.min(x)])) %>% 
  #dplyr::select(best_model) %>% 
  group_by(best_model) %>% 
  summarize(best_model = best_model, count = n()) %>% 
  unique()

count_reps_ou_nolam <- ou_gics %>% 
  mutate(best_model = apply(.[,2:4], 1, function(x) names(x)[which.min(x)])) %>% 
  #dplyr::select(best_model) %>% 
  group_by(best_model) %>% 
  summarize(best_model = best_model, count = n()) %>% 
  unique()

#BM redo data

#load
gic_comparison_bm_df_bm_redo <- readRDS("./analysis/from_cluster/full/gic_comparison_bm_df_bm_redo.rds")
gic_comparison_bm_df_lam_redo <- readRDS("./analysis/from_cluster/full/gic_comparison_bm_df_lam_redo.rds")
gic_comparison_bm_df_ou_redo <- readRDS("./analysis/from_cluster/full/gic_comparison_bm_df_ou_redo.rds")
gic_comparison_bm_df_eb_redo <- readRDS("./analysis/from_cluster/full/gic_comparison_bm_df_eb_redo.rds")

#combine
bm_gics_redo <- gic_comparison_bm_df_bm_redo

bm_gics_redo$LAM <- gic_comparison_bm_df_lam_redo$LAM
bm_gics_redo$OU <- gic_comparison_bm_df_ou_redo$OU
bm_gics_redo$EB <- gic_comparison_bm_df_eb_redo$EB

#get smalles number - put column name into new column
count_reps_bm_redo <- bm_gics_redo %>% 
  mutate(best_model = apply(.[,2:5], 1, function(x) names(x)[which.min(x)])) %>% 
  #dplyr::select(best_model) %>% 
  group_by(best_model) %>% 
  summarize(best_model = best_model, count = n()) %>% 
  unique()

count_reps_bm_redo_nolam <- bm_gics_redo %>% 
  mutate(best_model = apply(.[,2:4], 1, function(x) names(x)[which.min(x)])) %>% 
  #dplyr::select(best_model) %>% 
  group_by(best_model) %>% 
  summarize(best_model = best_model, count = n()) %>% 
  unique()

#OUM data

#load
gic_comparison_oum_df_bm <- readRDS("./analysis/from_cluster/full/gic_comparison_oum_df_bm.rds")
gic_comparison_oum_df_lam <- readRDS("./analysis/from_cluster/full/gic_comparison_oum_df_lam.rds")
gic_comparison_oum_df_ou <- readRDS("./analysis/from_cluster/full/gic_comparison_oum_df_ou.rds")
gic_comparison_oum_df_eb <- readRDS("./analysis/from_cluster/full/gic_comparison_oum_df_eb.rds")

#combine
oum_gics <- gic_comparison_oum_df_bm

oum_gics$LAM <- gic_comparison_oum_df_lam$LAM
oum_gics$OU <- gic_comparison_oum_df_ou$OU
oum_gics$EB <- gic_comparison_oum_df_eb$EB

#get smalles number - put column name into new column
count_reps_oum <- oum_gics %>% 
  mutate(best_model = apply(.[,2:5], 1, function(x) names(x)[which.min(x)])) %>% 
  #dplyr::select(best_model) %>% 
  group_by(best_model) %>% 
  summarize(best_model = best_model, count = n()) %>% 
  unique()

count_reps_oum_nolam <- oum_gics %>% 
  mutate(best_model = apply(.[,2:4], 1, function(x) names(x)[which.min(x)])) %>% 
  #dplyr::select(best_model) %>% 
  group_by(best_model) %>% 
  summarize(best_model = best_model, count = n()) %>% 
  unique()

#BMM data

#load
gic_comparison_bmm_df_bm <- readRDS("./analysis/from_cluster/full/gic_comparison_bmm_df_bm.rds")
gic_comparison_bmm_df_lam <- readRDS("./analysis/from_cluster/full/gic_comparison_bmm_df_lam.rds")
gic_comparison_bmm_df_ou <- readRDS("./analysis/from_cluster/full/gic_comparison_bmm_df_ou.rds")
gic_comparison_bmm_df_eb <- readRDS("./analysis/from_cluster/full/gic_comparison_bmm_df_eb.rds")

#combine
bmm_gics <- gic_comparison_bmm_df_bm

bmm_gics$LAM <- gic_comparison_bmm_df_lam$LAM
bmm_gics$OU <- gic_comparison_bmm_df_ou$OU
bmm_gics$EB <- gic_comparison_bmm_df_eb$EB

#get smalles number - put column name into new column
count_reps_bmm <- bmm_gics %>% 
  mutate(best_model = apply(.[,2:5], 1, function(x) names(x)[which.min(x)])) %>% 
  #dplyr::select(best_model) %>% 
  group_by(best_model) %>% 
  summarize(best_model = best_model, count = n()) %>% 
  unique()

#look at output
count_reps_ou
count_reps_bm

#incomplete
count_reps_oum
count_reps_bm_redo

#looking for parameters in ang_conifer output models
str(fit_1_ac_error)
str(fit_1_ac_error$sigma)
str(fit_1_ac_error$sigma$S)
sigma_ang_conifer <- fit_1_ac_error$sigma$S
sigma_ang_conifer[,1:10]

identical(fit_1_ac_error$sigma$Pinv, fit_1_ac_error$sigma$S)
all.equal(fit_1_ac_error$sigma$Pinv, fit_1_ac_error$sigma$S)


#sum the diagonals of the matrix? to get rate (sigma squared)

#from Shi et al 2021 code
#variance, sigma^2 in diag
#covariance, sigma in off-diag
rates <- vector("list", 1)
rates[[1]] <- sum(diag(one_partition_model$sigma)) #[PC, PC] is diagonal, so always is the rate of trait in question

#try with Pi then try with S
rate_sigma_squared_Pi <- sum(diag(fit_1_ac_error$sigma$Pinv))
rate_sigma_squared_S <- sum(diag(fit_1_ac_error$sigma$S))

#they are the same (in this case)

#question for Simon:
#sigma of mvgls object - includes Pi (covariance), P (precision) matrices, and S (sample estimates)
#What are the differences between those?
#Am I calculating sigma squared correctly this way (sum of diagonals of covariance matrix)?
#There is also corrSt but that seems to be the transformed tree


#compare rates for different models
BM_error_sigma_squared_Pi <- sum(diag(fit_1_ac_error$sigma$Pinv))
BM_orig_sigma_squared_Pi <- sum(diag(fit_1_ac$sigma$Pinv))
OU_error_sigma_squared_Pi <- sum(diag(fit_2_ac_error$sigma$Pinv))
OU_orig_sigma_squared_Pi <- sum(diag(fit_2_ac$sigma$Pinv))

BM_error_sigma_squared_Pi #0.05888094
BM_orig_sigma_squared_Pi #130.9017
OU_error_sigma_squared_Pi #0.02143446
OU_orig_sigma_squared_Pi #243.513

str(fit_1_ac)

#trying to figure out mvMORPH paper sigma (from Uyeda et al 2015)
p <- 20
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

Sigma_mv <- Posdef(p)
n
Z <- matrix(ncol=p, rnorm(p^2))

#seems like sigma matrix is a normally distributed random matrix with dimensions of the number of traits

#same for both
BM_self_sigma_squared_Pi <- sum(diag(fit_BM_self$sigma$Pinv))
BM_self_sigma_squared_S <- sum(diag(fit_BM_self$sigma$S))

str(fit_BM_self$sigma)
