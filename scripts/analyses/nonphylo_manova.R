#Non-phylo Manova

library(dplyr)
library(ggplot2)
library(factoextra)

#first check assumptions

#1. mshapiro.test() - to test for normal distribution within groups - multivariate normality
#2. homogeneity of variances across predictors (just use var?)
#3. linearity between pairs of dependent variables, pairs of covariates, and all dependent-covariate pairs

#import data

myc_data <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

new_spectra <- readRDS("./data/tidy/new_spectra_matched_trees.rds")

#check data

#dplyr::sample_n(new_spectra, 10)

res.man <- manova(spectra[,c(1:7)] ~ myc, data = myc_data)

summary(res.man)

summary.aov(res.man)

#any more than 7 wavelengths gives error about residuals have rank 7 < 8
#this is because of lack of independence between wavelengths
# #doesn't work 
# library(compositions)
# mv_out <- manova(ilr(clo(spectra)) ~ myc, data = myc_data)
# summary(mv_out)

#alternatively, first do non-phylo pca and then manovas on pcas

#calculate regular PCA
spectra.pca <- prcomp(myc_data$spectra, center = TRUE, scale. = TRUE)

#plot PCAs - regular
jpeg("./output/PCAs/nonphylo_myc_pc1_pc2.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = 1:2, ellipse=TRUE,  groups=myc_data$myc, var.axes=FALSE, size = 10, labels = dimnames(myc_data$spectra)[[1]]) +
  geom_point(aes(colour=myc_data$myc), size = 2)+
  scale_color_manual(name = "Mycorrhizal type", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

jpeg("./output/PCAs/nonphylo_myc_pc3_pc4.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = 3:4, ellipse=TRUE,  groups=myc_data$myc, var.axes=FALSE, size = 10, labels = dimnames(myc_data$spectra)[[1]]) +
  geom_point(aes(colour=myc_data$myc), size = 2)+
  scale_color_manual(name = "Mycorrhizal type", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

jpeg("./output/PCAs/nonphylo_myc_pc5_pc6.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = 5:6, ellipse=TRUE,  groups=myc_data$myc, var.axes=FALSE, size = 10, labels = dimnames(myc_data$spectra)[[1]]) +
  geom_point(aes(colour=myc_data$myc), size = 2)+
  scale_color_manual(name = "Mycorrhizal type", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

#manova on pc axes
str(spectra.pca)

pca_data_frame <- list(pcas = spectra.pca$x, myc = myc_data$myc)

pca_data_frame$pcas[,c(1:6)]

#test multivariate comparison of AM vs EM for top 20 pc axes
res.man <- manova(pcas[,c(1:20)] ~ myc, data = pca_data_frame)

summary(res.man) #pvalue - 0.0004482
summary.aov(res.man)

sum_res <- summary(res.man)
colnames(sum_res$stats)
#get rank of pcs which are most sig to least sig

#get output into new dataframe
df_output <- data.frame(matrix(nrow=20, ncol = 4))
colnames(df_output)[1] <- "iteration"
class(df_output[,1]) <- "numeric"
colnames(df_output)[2] <- "pvalue"
class(df_output[,2]) <- "numeric"
colnames(df_output)[3] <- "df_residuals"
class(df_output[,3]) <- "numeric"
colnames(df_output)[4] <- "F_value"
class(df_output[,4]) <- "numeric"
  
for (i in 1:length(summary.aov(res.man))){
  pvalue <- summary.aov(res.man)[[i]][["Pr(>F)"]][1]
  df_resids <- summary.aov(res.man)[[i]][["Df"]][2]
  F_value <- summary.aov(res.man)[[i]][["F value"]][1]
  df_output$iteration[i] <- i
  df_output$pvalue[i] <- pvalue
  df_output$df_residuals[i] <- df_resids
  df_output$F_value[i] <- F_value
}

#sort by magnitude of pvalue
sorted_pc_axes_pvalues <- df_output %>% 
  arrange(pvalue)

saveRDS(sorted_pc_axes_pvalues, "./analysis/nonphylo_manova/summary_manova_pc_results.rds")

#just looking at top 15 axes
res.man_top_15 <- manova(pcas[,c(1:15)] ~ myc, data = pca_data_frame)
summary(res.man_top_15)

#in general, very significant relationship - in nonphylo method

#most significant pcs - 1,3,12,5

#view pcas for these axes

#plot PCAs - regular
jpeg("./output/PCAs/nonphylo_myc_pc1_pc3_sig.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = c(1,3), ellipse=TRUE,  groups=myc_data$myc, var.axes=FALSE, size = 10, labels = dimnames(myc_data$spectra)[[1]]) +
  geom_point(aes(colour=myc_data$myc), size = 2)+
  scale_color_manual(name = "Mycorrhizal type", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

jpeg("./output/PCAs/nonphylo_myc_pc5_pc12_sig.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = c(5,12), ellipse=TRUE,  groups=myc_data$myc, var.axes=FALSE, size = 10, labels = dimnames(myc_data$spectra)[[1]]) +
  geom_point(aes(colour=myc_data$myc), size = 2)+
  scale_color_manual(name = "Mycorrhizal type", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

#plot loadings

#simple scree plot using factoextra
fviz_eig(spectra.pca)

summ <- summary(spectra.pca)

summ$importance

spectra.pca$rotation

vectors_df <- spectra.pca$rotation %>% 
  as.data.frame() 

vectors_df$wavelength <- c(400:2400)
vectors_for_plotting <- vectors_df %>% tidyr::gather(., pc_axis, value, -wavelength)
vectors_for_plotting$pc_axis <- gsub("PC", "", vectors_for_plotting$pc_axis)
class(vectors_for_plotting$pc_axis) <- "numeric"

pc_labs <- paste("PC", c(1:15))
names(pc_labs) <- c(1:15)

jpeg(paste0("./output/PCAs/nonphylo_myc_loadings_top15.jpg"), res = 400, width = 15, height = 15, units = "in")
ggplot(vectors_for_plotting[which(vectors_for_plotting$pc_axis %in% c(1:15)),], aes(x=wavelength,y=value)) + 
  geom_bar(stat = "identity")+
  facet_wrap(vars(pc_axis), ncol = 3, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
  theme(strip.text.x = element_text(size = 15))+
  xlab("Wavelength (nm)")+
  ylab("Vector unit")+
  theme_bw()
dev.off()

pc_labs <- c("PC1", "PC3", "PC5", "PC12")
names(pc_labs) <- c(1,3,5,12)

jpeg(paste0("./output/PCAs/nonphylo_myc_loadings_sig_pcs.jpg"), res = 400, width = 10, height = 8, units = "in")
ggplot(vectors_for_plotting[which(vectors_for_plotting$pc_axis %in% c(1,3,5,12)),], aes(x=wavelength,y=value)) + 
  geom_bar(stat = "identity")+
  facet_wrap(vars(pc_axis), ncol = 2, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
  theme(strip.text.x = element_text(size = 15))+
  xlab("Wavelength (nm)")+
  ylab("Vector unit")+
  theme_bw()
dev.off()


######ADONIS######
library(vegan)

adonis(dune ~ A1, data = dune.env,strata=dune.env$Management) 


adonis(dune ~ A1, data = dune.env)
