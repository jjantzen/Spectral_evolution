#looking at variation in spectra across sites

#load libraries
library(ggplot2)
library(ggbiplot)
library(dplyr)
library(spectrolab)
library(stringr)

#read in spectral data (not means)

spectra <- readRDS("./data/tidy/new_spectra_matched_trees.rds")
str(spectra)

#reshape data to keep spectra and species and sites
spectra_reshaped <- spectra
spectra_reshaped$meta <- spectra_reshaped$meta[,c(23,28,89,94)]

#calculate regular PCA for full dataset
spectra.pca <- prcomp(spectra_reshaped$value, center = TRUE, scale. = TRUE)

#plot coloured by site
jpeg("./output/PCAs/sample_based_nonphylo_pc1_pc2.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = 1:2, ellipse=TRUE,  groups=spectra_reshaped$meta$site, var.axes=FALSE, size = 10, labels = dimnames(spectra_reshaped$meta$site)[[1]]) +
  geom_point(aes(colour=spectra_reshaped$meta$site), size = 2)+
  #scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

#trim down to one species with lots of samples
list_sites <- as.character(unique(spectra_reshaped$meta$site))
list_sites

spectra_reshaped$meta$location <- as.character(spectra_reshaped$meta$site)

spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "CF_MSB")] <- "MSB"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "MSB_Chemin.L-d-Bou")] <- "MSB"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "MSB_Chem.L-d-Bou")] <- "MSB"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "MSB_Tourbiere")] <- "MSB"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "MSB_Mil.Humide-f")] <- "MSB"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "MSB_Lac-Boul.-f")] <- "MSB"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "MSB_lac_seigneur")] <- "MSB"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "MBP_20N")] <- "MBP"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "MBP_open")] <- "MBP"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "MBP_MSc_RBR")] <- "MBP"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "MBP_5N")] <- "MBP"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "MBP_10N")] <- "MBP"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "Oka_Plage")] <- "Oka"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "Oka_Humide")] <- "Oka"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "Oka_Plage-f")] <- "Oka"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "Oka_parking-f")] <- "Oka"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "SBLUdeM-f")] <- "SBLUdeM"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "SBLUdeM")] <- "SBLUdeM"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "SBLUdeM")] <- "SBLUdeM"
spectra_reshaped$meta$location[which(spectra_reshaped$meta$location == "CGOP_1")] <- "CGOP"

unique(spectra_reshaped$meta$location)

#plot coloured by location
jpeg("./output/PCAs/sample_based_location_nonphylo_pc1_pc2.jpg", res = 300, width = 20, height = 20, units = "in")
ggbiplot_edited(spectra.pca, choices = 1:2, ellipse=TRUE,  groups=spectra_reshaped$meta$location, var.axes=FALSE, size = 10, labels = dimnames(spectra_reshaped$meta$location)[[1]]) +
  geom_point(aes(colour=spectra_reshaped$meta$location), size = 2)+
  #scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

jpeg("./output/PCAs/sample_based_species_nonphylo_pc1_pc2.jpg", res = 300, width = 20, height = 20, units = "in")
ggbiplot_edited(spectra.pca, choices = 1:2, ellipse=TRUE,  groups=spectra_reshaped$meta$species_names, var.axes=FALSE, size = 10, labels = dimnames(spectra_reshaped$meta$species_names)[[1]]) +
  geom_point(aes(colour=spectra_reshaped$meta$species_names), size = 2)+
  #scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()


#plot by geographic region not by site
spectra_reshaped$meta$region <- as.character(spectra_reshaped$meta$site)
spectra_reshaped$meta$region[which(spectra_reshaped$meta$region != "CGOP_1")] <- "East_coast"
spectra_reshaped$meta$region[which(spectra_reshaped$meta$region == "CGOP_1")] <- "West_coast"

#plot by genus not by species
spectra_reshaped$meta$genus <- word(as.character(spectra_reshaped$meta$species_names), 1)

#plot new ones
jpeg("./output/PCAs/sample_based_regions_nonphylo_pc1_pc2.jpg", res = 300, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = 1:2, ellipse=TRUE,  groups=spectra_reshaped$meta$region, var.axes=FALSE, size = 10, labels = dimnames(spectra_reshaped$meta$region)[[1]]) +
  geom_point(aes(colour=spectra_reshaped$meta$region), size = 2)+
  #scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

jpeg("./output/PCAs/sample_based_genus_nonphylo_pc1_pc2.jpg", res = 300, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca, choices = 1:2, ellipse=TRUE,  groups=spectra_reshaped$meta$genus, var.axes=FALSE, size = 10, labels = dimnames(spectra_reshaped$meta$genus)[[1]]) +
  geom_point(aes(colour=spectra_reshaped$meta$genus), size = 2)+
  #scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

unique(as.character(spectra_reshaped$meta$region))

#choose top species
spectra_reshaped$meta$species_names

# spectra_matrix <- as.matrix(spectra_df)
# str(spectra_matrix) #not how I got the other format
spectra_df <- as.data.frame(spectra_reshaped)

species_counts <- spectra_df %>% 
  group_by(species_names) %>%
  dplyr::summarise(number = n()) %>% 
  dplyr::arrange(desc(number))

#use Populus tremuloides, Betula populifolia, and Acer rubrum
keepers <- c(species_counts$species_names[c(1:3)])

#subset dataframe to species - doesn't work
#spectra_populus <- spectra_reshaped[which(spectra_reshaped$meta$species_names %in% keepers[1])][1:4]


colnames(spectra_df)

spectra_df_populus <- spectra_df[which(spectra_df$species_names %in% keepers[1]),]
populus_spec <- as_spectra(spectra_df_populus, name_idx = 1, meta_idxs = c(2:8))

spectra_df_betula <- spectra_df[which(spectra_df$species_names %in% keepers[2]),]
betula_spec <- as_spectra(spectra_df_betula, name_idx = 1, meta_idxs = c(2:8))

spectra_df_acer <- spectra_df[which(spectra_df$species_names %in% keepers[3]),]
acer_spec <- as_spectra(spectra_df_acer, name_idx = 1, meta_idxs = c(2:8))

str(populus_spec)

#do pca for smaller set of taxa
spectra.pca_populus <- prcomp(populus_spec, center = TRUE, scale. = TRUE)

spectra.pca_betula <- prcomp(betula_spec, center = TRUE, scale. = TRUE)

spectra.pca_acer <- prcomp(acer_spec, center = TRUE, scale. = TRUE)

#plot coloured by location
jpeg("./output/PCAs/populus_by_site_nonphylo_pc1_pc2.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca_populus, choices = 1:2, ellipse=TRUE,  groups=populus_spec$meta$location, var.axes=FALSE, size = 10, labels = dimnames(populus_spec$meta$location)[[1]]) +
  geom_point(aes(colour=populus_spec$meta$location), size = 2)+
  #scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

jpeg("./output/PCAs/betula_by_site_nonphylo_pc1_pc2.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca_betula, choices = 1:2, ellipse=TRUE,  groups=betula_spec$meta$location, var.axes=FALSE, size = 10, labels = dimnames(betula_spec$meta$location)[[1]]) +
  geom_point(aes(colour=betula_spec$meta$location), size = 2)+
  #scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()

jpeg("./output/PCAs/acer_by_site_nonphylo_pc1_pc2.jpg", res = 600, width = 10, height = 10, units = "in")
ggbiplot_edited(spectra.pca_acer, choices = 1:2, ellipse=TRUE,  groups=acer_spec$meta$location, var.axes=FALSE, size = 10, labels = dimnames(acer_spec$meta$location)[[1]]) +
  geom_point(aes(colour=acer_spec$meta$location), size = 2)+
  #scale_color_manual(name = "Woody", values = c("mediumpurple1", "darkgreen")) +
  theme_bw()#labels=rownames(spectra$trait), 
dev.off()


#### do phylo pcas for site
#load model
model <- readRDS("./data/physignal/fit_ou_intercept.rds")

#load spectra
spectra <- readRDS("./data/tidy/spectra_not_reordered_to_tree.rds")

#load trait data
trait_data <- readRDS("./data/tidy/new_combo_data_matched_spectra.rds")

#make data object
data_spectra <- list(trait=spectra, woody = trait_data$Woody, order = trait_data$Order)

#calculate phylo pca for model
pca_phylo <- mvgls.pca(model, plot=FALSE)

#get grouping variable as dataframe
species_by_region <- spectra_reshaped$meta %>% 
  #dplyr::select(-sample_name) %>% 
  dplyr::distinct(species_names, region) %>% 
  dplyr::arrange(species_names)

region <- as.data.frame(species_by_region$region, row.names = species_by_region$species_names)

#assign numbers for shape by group
pch.group <- as.numeric(region$`species_by_region$region`)
pch.group <- gsub("1", "15", pch.group)
pch.group <- gsub("2", "19", pch.group)

#assign colours to group
col.group <- as.numeric(region$`species_by_region$region`)
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
    #geom_path(data = ell[which(ell$groups == 0),])#, aes(colour = groups, group = groups))+ #, inherit.aes = FALSE
    geom_hline(aes(yintercept = 0))+
    geom_vline(aes(xintercept = 0))+
    labs(y = ylabel, x = xlabel)+
    scale_colour_manual(name = "Region", labels = c("East coast", "West coast"), values = rev(levels(factor(cols))))+
    theme_bw()
}

#plot for sets of axes
jpeg("./output/PCAs/intercept_phylo_pc1_pc2_by_region.jpg", res = 600, width = 10, height = 10, units = "in")
plot_pca_phylo(pca_phylo, model, axes = c(1,2), cols = col.group, groups = region$`species_by_region$region`, labels = rownames(region))
dev.off()


