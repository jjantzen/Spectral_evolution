#Do PCA of spectra - phylo and reg

#libraries to load
library(ggplot2)
library(ggbiplot)
library(mvMORPH)
library(nlme)
library(vegan)

#load model
model <- readRDS("./data/physignal/fit_ou_intercept.rds")

#load spectra
spectra <- readRDS("./data/tidy/spectra_not_reordered_to_tree.rds")

#load trait data
trait_data <- readRDS("./data/tidy/new_combo_data_matched_spectra.rds")

#make data object
data_spectra <- list(trait=spectra, woody = trait_data$Woody, order = trait_data$Order)
str(data_spectra)


##################
#calculate regular PCA
spectra.pca <- prcomp(data_spectra$trait, center = TRUE, scale. = TRUE)

#plot PCAs - regular
pdf("./output/PCAs/intercept_nonphylo_pc1_pc2.pdf", width = 10, height = 10)
ggbiplot(spectra.pca, ellipse=TRUE,  groups=data_spectra$Woody, var.axes=FALSE)#labels=rownames(spectra$trait), 
dev.off()

jpeg("./output/PCAs/intercept_nonphylo_pc1_pc2.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = 1:2, ellipse=TRUE,  groups=data_spectra$woody, var.axes=FALSE, size = 10, labels = dimnames(data_spectra$trait)[[1]]) +
  geom_point(aes(colour=data_spectra$woody), size = 2)+
  scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

jpeg("./output/PCAs/intercept_nonphylo_pc3_pc4.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = 3:4, ellipse=TRUE,  groups=data_spectra$woody, var.axes=FALSE, size = 10, labels = dimnames(data_spectra$trait)[[1]]) +
  geom_point(aes(colour=data_spectra$woody), size = 2)+
  scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait),
dev.off()

jpeg("./output/PCAs/intercept_nonphylo_pc5_pc6.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca,choices = 5:6, ellipse=TRUE,  groups=data_spectra$woody, var.axes=FALSE, size = 10, labels = dimnames(data_spectra$trait)[[1]]) +
  geom_point(aes(colour=data_spectra$woody), size = 2)+
  scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait),
dev.off()

jpeg("./output/PCAs/intercept_nonphylo_pc7_pc8.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = 7:8, ellipse=TRUE,  groups=data_spectra$woody, var.axes=FALSE, size = 10, labels = dimnames(data_spectra$trait)[[1]]) +
  geom_point(aes(colour=data_spectra$woody), size = 2)+
  scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait),
dev.off()

jpeg("./output/PCAs/intercept_nonphylo_pc1_pc3.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = c(1,3), ellipse=TRUE,  groups=data_spectra$woody, var.axes=FALSE, size = 10, labels = dimnames(data_spectra$trait)[[1]]) +
  geom_point(aes(colour=data_spectra$woody), size = 2)+
  scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait),
dev.off()

jpeg("./output/PCAs/intercept_nonphylo_pc2_pc4.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = c(2,4), ellipse=TRUE,  groups=data_spectra$woody, var.axes=FALSE, size = 10, labels = dimnames(data_spectra$trait)[[1]]) +
  geom_point(aes(colour=data_spectra$woody), size = 2)+
  scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait),
dev.off()

jpeg("./output/PCAs/intercept_nonphylo_pc3_pc5.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = c(3,5), ellipse=TRUE,  groups=data_spectra$woody, var.axes=FALSE, size = 10, labels = dimnames(data_spectra$trait)[[1]]) +
  geom_point(aes(colour=data_spectra$woody), size = 2)+
  scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait),
dev.off()


##################
#calculate phylo pca for model
pca_phylo <- mvgls.pca(model, plot=FALSE)

# pdf("./output/PCAs/intercept_pc1_pc2_default.pdf", width = 10, height = 10)
# mvgls.pca_ed(model, plot=TRUE)#col = "black",, pch = pch.group, col=col.group
# dev.off()

# head(pca_phylo$scores[1:100])
# plot(pca_phylo$scores, col = c("blue", "green"))#levels(fit_2_ac_error$variables$X[,2]))
# plot(pca_phylo$scores[,1], pca_phylo$scores[,3])

#get grouping variable as dataframe
woodiness <- as.data.frame(trait_data$Woody, row.names = trait_data$Species)

#assign numbers for shape by group
pch.group <- as.numeric(woodiness$`trait_data$Woody`)
pch.group <- gsub("1", "15", pch.group)
pch.group <- gsub("2", "19", pch.group)

#assign colours to group
col.group <- as.numeric(woodiness$`trait_data$Woody`)
col.group <- gsub("1", "mediumpurple1", col.group)
col.group <- gsub("2", "darkgreen", col.group)

# rownames(model$variables$X)
# rownames(woodiness)
# col.group <- data_spectra$woody
# levels(col.group) <- c("darkgreen", "mediumpurple1")

#plotting function - figure out ellipses
plot_pca_phylo <- function(pca_object, model, axes, groups,  cols, labels) {
  #getting values for plotting
  U <- pca_object$vectors
  resids <- model$residuals
  S <- resids %*% U
  #make as dataframe for plotting and name columns
  df_S <- as.data.frame(S[,axes[1:2]])
  names(df_S) <- c("xvar", "yvar")
  #assign groups based on separate trait for colouring
  #groups <- levels(factor(col))
  df_S$groups <- groups
  #get labels for axes
  tot <- sum(pca_object$values)
  valX <- round(pca_object$values[axes[1]] * 100/tot, digits = 2)
  valY <- round(pca_object$values[axes[2]] * 100/tot, digits = 2)
  xlabel <- paste("PC", axes[1], " (", valX, " %)", sep = "")
  ylabel <- paste("PC", axes[2], " (", valY, " %)", sep = "")
  #get parameters for ellipses
  ellipse.prob <- 0.68
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  ell <- ddply(df_S, "groups", function(x) {
    sigma <- var(cbind(df_S$xvar, df_S$yvar))
    mu <- c(mean(df_S$xvar), mean(df_S$yvar))
    ed <- sqrt(qchisq(ellipse.prob, df = 2))
    data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = "+"), groups = df_S$groups)
  })
  names(ell)[1:2] <- c("xvar", "yvar")
  str(ell)
  #labels for points
  df_S$labels <- labels
  ggplot(data = df_S, aes(x = xvar, y = yvar))+
    geom_point(aes(colour = groups), shape = 19, size = 2)+
    geom_text_repel(label = labels, max.overlaps = 3, size = 3)+
    geom_path(data = ell[which(ell$groups == 0),])#, aes(colour = groups, group = groups))+ #, inherit.aes = FALSE
    geom_hline(aes(yintercept = 0))+
    geom_vline(aes(xintercept = 0))+
    labs(y = ylabel, x = xlabel)+
    scale_colour_manual(name = "Woodiness", labels = c("Non-woody", "Woody"), values = rev(levels(factor(cols))))+
    theme_bw()
}

plot(ell$xvar)
plot(circle)

pca_object <- pca_phylo

plot_pca_phylo(pca_object, model, axes, groups,  cols, labels)

#different method for ellipses
ordiellipse(pca_object,groups,conf=0.65)


#do pca first using modified function
pca_phylo <- mvgls.pca_ed(model) #, axes = c(1,2)

#testing out ellipeses
plot_pca_phylo(pca_phylo, model, axes = c(1,2), cols = col.group, groups = trait_data$Woody, labels = trait_data$Species)


str(pca_phylo)



#plot for sets of axes
jpeg("./output/PCAs/intercept_phylo_pc1_pc2.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(pca_phylo, model, axes = c(1,2), col = col.group, groups = data_spectra$woody)
dev.off()

jpeg("./output/PCAs/intercept_phylo_pc3_pc4.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(pca_phylo, model, axes = c(3,4), col.group)
dev.off()

jpeg("./output/PCAs/intercept_phylo_pc5_pc6.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(pca_phylo, model, axes = c(5,6), col.group)
dev.off()

jpeg("./output/PCAs/intercept_phylo_pc7_pc8.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(pca_phylo, model, axes = c(7,8), col.group)
dev.off()

jpeg("./output/PCAs/intercept_phylo_pc1_pc3.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(pca_phylo, model, axes = c(1,3), col.group)
dev.off()

jpeg("./output/PCAs/intercept_phylo_pc2_pc4.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(pca_phylo, model, axes = c(2,4), col.group)
dev.off()

jpeg("./output/PCAs/intercept_phylo_pc3_pc5.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(pca_phylo, model, axes = c(3,5), col.group)
dev.off()

jpeg("./output/PCAs/intercept_phylo_pc1_pc2.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(pca_phylo, model, axes = c(1,2), col.group)
dev.off()



#ggbiplot(pca_ang_conifer$scores, ellipse=TRUE,  groups=spectra$lineage, var.axes=FALSE)
#pch = pch.group,, plot=TRUE,  col = col.group



#non-function version of plotting pca output object
ggplot(data = df_S, aes(x = V1, y = V2))+
  geom_point(aes(colour = factor(col.group)), shape = 19, size = 2)+
  scale_colour_manual(name = "Woodiness", labels = c("Woody", "Non-woody"), values = c("darkgreen", "mediumpurple1"))+
  geom_hline(aes(yintercept = 0))+
  geom_vline(aes(xintercept = 0))+
  labs(y = ylabel, x = xlabel)+
  geom_text_repel(label = rownames(df_S), max.overlaps = 3, size = 2)+
  theme_bw()


###################


pdf("./output/ang_conifer_pc1_pc2.pdf", width = 10, height = 10)
mvgls.pca(fit_2_ac_error, plot=TRUE)
dev.off()

pdf("./output/ang_conifer_pc3_pc4.pdf", width = 10, height = 10)
mvgls.pca(fit_2_ac_error, plot=TRUE, axes= c(3,4))
dev.off()

fit_eb_intercept$param
