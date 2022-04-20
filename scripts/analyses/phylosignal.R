#Phylosignal of spectra

library(geomorph)
library(mvMORPH)
library(phytools)

#calculate Kmulti using physignal function - for procrustes shape variables


spectra_matrix <- readRDS("./data/tidy/spectra_not_reordered_to_tree.rds")

str(spectra_matrix)

data_spectra <- list(trait=spectra_matrix)

new_trees <- readRDS("./data/tidy/new_trees_matched_spectra.rds")

#rescaled_trees <- readRDS("./data/tidy/rescaled_trees_matched_spectra.rds")

tree1 <- new_trees[[1]]
#tree1 <- rescaled_trees[[1]]

class(tree1)

dimnames(spectra_matrix)
tree1$tip.label

#A is matrix or 3D array
#phy is tree of class phylo
#physignal(A, phy, iter = 9, seed = NULL, print.progress = TRUE) #999

#calculates p value with significance

#full spectra
physignal(spectra_matrix, tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.0768 pvalue 0.003, ES 2.3865


#arbitrary ranges and individual wavelengths
physignal(spectra_matrix[,1:1000], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.0546 pvalue 0.075 ES 1.4163
physignal(spectra_matrix[,1:500], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.044 pvalue 0.15
physignal(spectra_matrix[,501:1000], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.0625 pvalue 0.067
physignal(spectra_matrix[,1001:2001], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.1503 pvalue 0.001
physignal(spectra_matrix[,1000:1100], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.2401 pvalue 0.001
physignal(spectra_matrix[,1000:1050], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.2614 pvalue 0.001
physignal(spectra_matrix[,1000:1200], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.1968 pvalue 0.001
physignal(spectra_matrix[,1150:1200], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.154 pvalue 0.001
physignal(spectra_matrix[,1050:1100], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.2215 pvalue 0.001
physignal(spectra_matrix[,1100:1150], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.188 pvalue 0.001
physignal(spectra_matrix[,1001:1500], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.1513 pvalue 0.001
physignal(spectra_matrix[,1501:2001], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.1485 pvalue 0.001

physignal(spectra_matrix[,1000], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.2528 pvalue 0.001 ES 4.5993
physignal(spectra_matrix[,1], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.0047 pvalue 0.969 ES -1.8952
physignal(spectra_matrix[,2], tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.0047 pvalue 0.971 ES -1.8947

#comparison with other method of calculating K and lambda
phylosig(tree1, spectra_matrix[,1000], method="lambda", test=TRUE, nsim=999) #0.877844 p value v small
phylosig(tree1, spectra_matrix[,1000], method="K", test=TRUE, nsim=999) #0.252761 pvalue 0.001001
phylosig(tree1, spectra_matrix[,5], method="K", test=TRUE, nsim=999) #0.00469549 pvalue 0.977978
phylosig(tree1, spectra_matrix[,5], method="lambda", test=TRUE, nsim=999) #6.03548e-05 pvalue 1
phylosig(tree1, spectra_matrix[,500], method="K", test=TRUE, nsim=999) #0.0489246 pvalue 0.168168
phylosig(tree1, spectra_matrix[,500], method="lambda", test=TRUE, nsim=999) #0.701759 pvalue 3.39795e-07
phylosig(tree1, spectra_matrix[,1500], method="K", test=TRUE, nsim=999) #0.157186 pvalue 0.001001
phylosig(tree1, spectra_matrix[,1500], method="lambda", test=TRUE, nsim=999) #0.824058 pvalue small
phylosig(tree1, spectra_matrix[,2000], method="K", test=TRUE, nsim=999) #0.148842 pvalue 0.001001
phylosig(tree1, spectra_matrix[,2000], method="lambda", test=TRUE, nsim=999) #0.82143 pvalue small
phylosig(tree1, spectra_matrix[,100], method="K", test=TRUE, nsim=999) #0.0146972 pvalue 0.746747
phylosig(tree1, spectra_matrix[,100], method="lambda", test=TRUE, nsim=999) #0.644987 pvalue small
phylosig(tree1, spectra_matrix[,750], method="K", test=TRUE, nsim=999) #0.0611665 pvalue 0.0860861
phylosig(tree1, spectra_matrix[,750], method="lambda", test=TRUE, nsim=999) #0.688226 pvalue small

####Getting phylogenetic coordinates and profiles
# #phy sig profile
# psig_full <- physignal(spectra_matrix, tree1, iter = 999, seed = NULL, print.progress = TRUE) #0.0768 pvalue 0.003, ES 2.3865
# summary(psig_full)
# jpeg("./output/kmult_distribution_full_spectra.jpg")
# plot(psig_full)
# dev.off()
# psig_full$K.by.p
# 
# plot(psig_full$PACA, phylo = TRUE)
# 
# psig_subset <- physignal(spectra_matrix[,1001:1500], tree1, iter = 999, seed = NULL, print.progress = TRUE)
# 
# plot(psig_subset)
# plot(psig_subset$PACA, phylo = TRUE)
# psig_subset$K.by.p


#for main regions of spectrum

#visible 400-699 nm
physignal(spectra_matrix[,1:300], tree1, iter = 999, seed = NULL, print.progress = TRUE) #K = 0.0332, p = 0.361

#NIR 700-1399 nm
physignal(spectra_matrix[,301:1000], tree1, iter = 999, seed = NULL, print.progress = TRUE) #K = 0.0572, p = 0.092

#SWIR 1400-2400 nm
swir <- physignal(spectra_matrix[,1001:2001], tree1, iter = 999, seed = NULL, print.progress = TRUE) #K = 0.1503, p = 0.001

colnames(spectra_matrix[,1001:2001])
# jpeg("./output/kmult_distribution_swir.jpg")
# plot(swir)
# dev.off()

#try it for PCs - redo with latest models

fit_2_ac_error <- readRDS("./data/ang_conifer/rds/fit_2_ac_error.rds")

pca_ang_conifer <- mvgls.pca(fit_2_ac_error, plot=TRUE)

head(pca_ang_conifer$scores)[,1:10]

tree1
dimnames(pca_ang_conifer$scores)[[1]][which(dimnames(pca_ang_conifer$scores)[[1]] %in% tree1$tip.label == FALSE)]

pc1_edited <- pca_ang_conifer$scores[-which(dimnames(pca_ang_conifer$scores)[[1]] == "Picea rubens"),]
tree1_edited <- tree1
tree1_edited <- drop.tip(tree1_edited, tree1_edited$tip.label[which(tree1_edited$tip.label %in% dimnames(pca_ang_conifer$scores)[[1]] == FALSE)])

physignal(pc1_edited[,9], tree1_edited, iter = 999, seed = NULL, print.progress = TRUE) #0.1147 pvalue 0.003, ES 2.3865

phylosig(tree1_edited, pc1_edited[,1], method="K", test=TRUE, nsim=999) #0.066653 p value 0.108108
phylosig(tree1_edited, pc1_edited[,9], method="K", test=TRUE, nsim=999) #0.114683 p value 0.002002

phylosig(tree1_edited, pc1_edited[,1], method="lambda", test=TRUE, nsim=999) #0.495778 p value 0.0099809
phylosig(tree1_edited, pc1_edited[,9], method="lambda", test=TRUE, nsim=999) #0.322394 p value 0.0404377

#get lambda for entire model
fit_lambda_intercept <- readRDS("./data/physignal/fit_lambda_intercept.rds")
fit_bm_intercept <- readRDS("./data/physignal/fit_bm_intercept.rds")
fit_ou_intercept <- readRDS("./data/physignal/fit_ou_intercept.rds")
fit_eb_intercept <- readRDS("./data/physignal/fit_eb_intercept.rds")

fit_ou_intercept_rescaled <- readRDS("./data/physignal/fit_ou_intercept_rescaled.rds")
fit_eb_intercept_rescaled <- readRDS("./data/physignal/fit_eb_intercept_rescaled.rds")


GIC(fit_eb_intercept_rescaled) #GIC: -2607702 | Log-likelihood 1305790 (biggest negative)
GIC(fit_ou_intercept_rescaled) #GIC: -2607570 | Log-likelihood 1305724 


GIC(fit_lambda_intercept) #GIC: -2540204 | Log-likelihood 1271836
GIC(fit_bm_intercept) #GIC: -2388897 | Log-likelihood 1195795
GIC(fit_ou_intercept) #GIC: -2677467 | Log-likelihood 1340925 (biggest negative so best model)
GIC(fit_eb_intercept) #GIC: -2646101 | Log-likelihood 1325125

fit_eb_intercept_rescaled$param #-5.545177 #a pretty strong slow down???

fit_lambda_intercept$param #lambda = 0.1999411
fit_bm_intercept$param #NA
alpha <- fit_ou_intercept$param #alpha = 0.02008654 (non-scaled)
fit_ou_intercept$sigma #gives matrices of 2001 rows

alpha_rs <- fit_ou_intercept_rescaled$param #alpha = 5.545177 (rescaled)

#interpretation like Ives and Garland 2010 from Cooper et al 2015?
-log(alpha_rs) #-1.712929 so higher rather than lower (negative is higher, positive is lower (BM))

fit_lambda_intercept$mserr #estimated standard error is 115.3676

#####
#THIS DOESN'T QUITE SEEM RIGHT - since halflife is in branch lengths, should use halflife, not time_to_opt as phy halflife
#doesn't match Artuso paper

#phylogenetic halflife
halflife <- log(2)/alpha #in units of branch lengths 

halflife_rs <- log(2)/alpha_rs #0.125 in units of branch lengths (so scaled to 1)



#get total tree length (crown age of clade)
tree_max <- max(node.depth.edgelength(tree1))

#time
time_to_opt <- halflife*tree_max

#halflife times total tree length (replicating Artuso paper)
halflife_art <- log(2)/1.1*11.5


#Rescale tre to height of 1, then interpret log(alpha): -log(alpha) = 4 is very low (BM), while -log(alpha)=-4 is very high
#or do phylogenetic halflife (see Hansen and Bartoszek, 2012 or Slater 2015) - from Cooper et al 2015






#For checking matching of tree and spectra taxa
# tree1$tip.label %in% dimnames(spectra_matrix)[[1]]
# 
# tree1$tip.label[duplicated(tree1$tip.label)]
