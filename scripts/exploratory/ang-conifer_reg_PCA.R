#Do regular PCA of spectra

#libraries to load
library(ggplot2)
library(ggbiplot)
library(mvMORPH)
library(nlme)

#load data

spectra <- readRDS("./data/ang_conifer/rds/data_spectra_ang_conifer.rds")
str(spectra)
spec.pca <- prcomp(spectra$trait, center = TRUE, scale. = TRUE)


pdf("./output/ang_conifer_nonphylo_pc1_pc2.pdf", width = 10, height = 10)
ggbiplot(spec.pca, ellipse=TRUE,  groups=spectra$lineage, var.axes=FALSE)#labels=rownames(spectra$trait), 
dev.off()

jpeg("./output/ang_conifer_nonphylo_pc1_pc2.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot(spec.pca, ellipse=TRUE,  groups=spectra$lineage, var.axes=FALSE, size = 10) +
  geom_point(aes(colour=spectra$lineage), size = 8)+
  scale_color_manual(name = "Lineage", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

fit_2_ac_error <- readRDS("./data/ang_conifer/rds/fit_2_ac_error.rds")
pca_ang_conifer <- mvgls.pca(fit_2_ac_error, plot=FALSE)

head(pca_ang_conifer$scores[1:100])
plot(pca_ang_conifer$scores, col = c("blue", "green"))#levels(fit_2_ac_error$variables$X[,2]))

str(fit_2_ac_error$variables$X)

str(pca_ang_conifer)
str(spec.pca)
pca_ang_conifer$scores
fit_2_ac_error
plot(pca_ang_conifer$scores[,1], pca_ang_conifer$scores[,2])


rownames(pca_ang_conifer$scores)

lineage_data <- as.data.frame(fit_2_ac_error$variables$X[,2], row.names = rownames(fit_2_ac_error$variables$X))
colnames(lineage_data) <- "data"
# 
# reorder_idx <- match(rownames(pca_ang_conifer$scores),rownames(lineage_data))
# lineage_data_reordered <- lineage_data[reorder_idx]  

str(fit_2_ac_error)

pch.group <- c(rep(17, times=sum(lineage_data$data == "1")), rep(16, times=sum(lineage_data$data == "0")))
col.group <- c(rep("darkgreen", times=sum(lineage_data$data == "1")), rep("mediumpurple1", times=sum(lineage_data$data == "0")))

#xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = "")

pdf("./output/ang_conifer_pc1_pc2_colours.pdf", width = 10, height = 10)
#ggbiplot(pca_ang_conifer$scores, ellipse=TRUE,  groups=spectra$lineage, var.axes=FALSE)
mvgls.pca_ed(fit_2_ac_error, plot=TRUE, pch = pch.group, col=col.group)#col = "black",
dev.off()

jpeg("./output/ang_conifer_pc1_pc2_colours.jpg", res = 600, width = 10, height = 10, units = "in")
#ggbiplot(pca_ang_conifer$scores, ellipse=TRUE,  groups=spectra$lineage, var.axes=FALSE)
mvgls.pca_ed(fit_2_ac_error, plot=TRUE, pch = pch.group, col=col.group)
dev.off()

#regular lm model
data_spectra <- readRDS("./data/ang_conifer/rds/data_spectra_ang_conifer.rds")


fit_1_ac <- gls(trait ~ 1, data = data_spectra)
summary(fit_1_ac)
fit_1_lm <- readRDS("./data/ang_conifer/rds/fit_1_lm_ac.rds")
res_lm_1 <- residuals(fit_1_ac)

summary(fit_1_lm)

res_lm <- residuals(fit_1_lm$residuals)

str(fit_1_lm$model)
fit_2_ac_error$model



full_tree <- readRDS("./data/ang_conifer/rds/pruned_ang_con_tree.rds")

jpeg("./output/tree_figure.jpg", height = 20, width = 8, units = "in", res = 600)
plot(ult_tree)
dev.off()

sig_resid <- phylosig(ult_tree, res_lm_1, method = "lambda", test = TRUE, nsim = 999)

res_phylo <- residuals(fit_2_ac_error)
str(res_phylo)
str(res_lm)

sig_resid_phylo <- phylosig(ult_tree, res_phylo, method = "lambda", test = TRUE, nsim = 999)

lm_manova <- manova(trait ~ lineage, data = data_spectra)

lm_aov <- aov(trait ~ lineage, data = data_spectra)

summary(lm_manova)

summary(lm_aov)



#create polytomy phylogeny with zero branch lengths
testing <- mvgls(trait ~ lineage, data=data_spectra, model="OU") 

#redo analyses with ultrametric tree
ult_tree <- force.ultrametric(full_tree, method=c("nnls"))

di2multi(full_tree)
