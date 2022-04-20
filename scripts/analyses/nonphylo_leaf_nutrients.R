#leaf nutrients non-phylo analysis

library(spectrolab)
library(dplyr)
library(phytools)
library(nlme)

#read data
spectral_data <- readRDS("./data/for_analysis/final_spectra_matched_trees.rds")

predictor_data <- readRDS("./data/for_analysis/final_data.rds")

myc_data <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

consensus_tree <- readRDS("./data/for_analysis/consensus_tree_full.rds")

#prune trees to match data - not needed yet (unless trim both to myc dataset)
exclude <- consensus_tree$tip.label[-which(consensus_tree$tip.label %in% myc_data$species)]
pruned_consensus <- drop.tip(consensus_tree, tip = exclude)

#prep data
metadata <- meta(spectral_data)
colnames(metadata)
trait_data <- metadata[,c(28, 86, 99, 106, 120, 123, 133, 161, 122)]
colnames(trait_data)

#convert to numeric from factor via character
trait_data[,c(3:9)] <- lapply(trait_data[,c(3:9)], function(x) as.numeric(as.character(x)))

#get mean of spectra
trait_means <- trait_data[,c(1,3:9)] %>% 
  #dplyr::select(-c(sample_name)) %>% 
  group_by(species_names) %>% 
  summarize_all(mean, na.rm = TRUE)

#add myc state to trait data
myc_data_df <- as.data.frame(x = myc_data$myc)
myc_data_df <- cbind(myc_data$species, myc_data_df, myc_data$lp, myc_data$gf)
colnames(myc_data_df) <- c("species_names", "myc", "lp", "gf")

trait_data_combo <- merge(trait_data, myc_data_df, by = "species_names", all = TRUE)

#figure out which data missing
trait_means$species_names[which(is.na(trait_means$P_area))]

###########For Nmass
#drop pinus rigida
pruned_tree_2 <- drop.tip(pruned_consensus, tip = "Pinus rigida")

trait_means_2 <- trait_means[which(trait_means$species_names %in% pruned_tree_2$tip.label),]

myc_data_df_2 <- myc_data_df[-which(myc_data_df$species_names == "Pinus rigida"),]

trait_data_combo_2 <- merge(trait_means_2, myc_data_df_2, by = "species_names", all = TRUE)

#non-phylo model: linear regression of trait and mycorrhizal type
model_Nmass <- lm(Nmass ~ myc, data = trait_data_combo_2)
Nmass_output <- summary(model_Nmass)

#phylo model
vcv_matrix <- vcv.phylo(pruned_tree_2, cor=TRUE)
colnames(trait_data_combo_2)[1] <- "Species"

phylo_Nmass <- gls(Nmass ~ myc, data=trait_data_combo_2, correlation=corPagel(1, pruned_tree_2, form = ~Species))

plot(phylo_Nmass)    # residual plot
summary(phylo_Nmass) # coefficients, SE's, fit statistics
confint(phylo_Nmass) # confidence intervals for coefficients

###########Lignin mass
#drop Alnus incana
pruned_tree_2 <- drop.tip(pruned_consensus, tip = "Alnus incana")

trait_means_2 <- trait_means[which(trait_means$species_names %in% pruned_tree_2$tip.label),]

myc_data_df_2 <- myc_data_df[-which(myc_data_df$species_names == "Alnus incana"),]

trait_data_combo_2 <- merge(trait_means_2, myc_data_df_2, by = "species_names", all = TRUE)

#non-phylo model: linear regression of trait and mycorrhizal type
model_ligninmass <- lm(lignin_mass ~ myc, data = trait_data_combo_2) 
ligninmass_output <- summary(model_ligninmass) #sig for both

#phylo model
vcv_matrix <- vcv.phylo(pruned_tree_2, cor=TRUE)

colnames(trait_data_combo_2)[1] <- "Species"

phylo_ligninmass <- gls(lignin_mass ~ myc, data=trait_data_combo_2, correlation=corPagel(1, pruned_tree_2, form = ~Species))

plot(phylo_ligninmass)    # residual plot
summary(phylo_ligninmass) # coefficients, SE's, fit statistics
confint(phylo_ligninmass) # confidence intervals for coefficients

#test correlation of mass vs area and lma

cor(trait_data_combo_2$Nmass, trait_data_combo_2$LMA) #-0.604046
cor(trait_data_combo_2$Narea, trait_data_combo_2$LMA) #0.8631537
