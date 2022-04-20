#comparing pcas from different models

#read in models

myc_model_trial_BM <- readRDS("./analysis/testing_models/testing_myc/myc_model_trial_BM.rds")
myc_model_trial_OU <- readRDS("./analysis/testing_models/testing_myc/myc_model_trial_OU.rds")
myc_model_trial_lambda <- readRDS("./analysis/testing_models/testing_myc/myc_model_trial_lambda.rds")

#do pcas on both
BM_pca <- mvgls.pca(myc_model_trial_BM, plot = FALSE)

OU_pca <- mvgls.pca(myc_model_trial_OU, plot = FALSE)

lam_pca <- mvgls.pca(myc_model_trial_lambda, plot = FALSE)

#compare pc 1 from each (are they)
col.group <- data_spectra$myc
col.group <- gsub("AM", "orange", col.group)
col.group <- gsub("EM", "blue", col.group)

#plot pc1 vs pc4
jpeg("./output/PCAs/myc_phylo_pc1_4_ou.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(OU_pca, myc_model_trial_OU, axes = c(1,4), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
dev.off()

#plot pc1 vs pc4
jpeg("./output/PCAs/myc_phylo_pc1_4_lambda.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(lam_pca, myc_model_trial_lambda, axes = c(1,4), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
dev.off()

#plot pc1 vs pc4
jpeg("./output/PCAs/myc_phylo_pc1_4_bm.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(BM_pca, myc_model_trial_BM, axes = c(1,4), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
dev.off()

#plot pc1 vs pc4
jpeg("./output/PCAs/myc_phylo_pc2_3_ou.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(OU_pca, myc_model_trial_OU, axes = c(2,3), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
dev.off()

#plot pc1 vs pc4
jpeg("./output/PCAs/myc_phylo_pc2_3_lambda.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(lam_pca, myc_model_trial_lambda, axes = c(2,3), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
dev.off()

#plot pc1 vs pc4
jpeg("./output/PCAs/myc_phylo_pc2_3_bm.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(BM_pca, myc_model_trial_BM, axes = c(2,3), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
dev.off()

identical(BM_pca$scores, OU_pca$scores)
all.equal(BM_pca$scores, OU_pca$scores)

identical(BM_pca$scores, lam_pca$scores)
all.equal(BM_pca$scores, lam_pca$scores)

#pcas from different models are different