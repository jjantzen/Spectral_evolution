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
trait_data <- metadata[,c(28, 86, 97, 99, 106, 120, 122, 123, 133, 161)]
colnames(trait_data)
str(trait_data)
#convert to numeric from factor via character
trait_data[,c(3:10)] <- lapply(trait_data[,c(3:10)], function(x) as.numeric(as.character(x)))

#get mean of spectra
trait_means <- trait_data[,c(1,3:10)] %>% 
  #dplyr::select(-c(sample_name)) %>% 
  group_by(species_names) %>% 
  summarize_all(mean, na.rm = TRUE)

#add myc state to trait data
myc_data_df <- as.data.frame(x = myc_data$myc)
myc_data_df <- cbind(myc_data$species, myc_data_df, myc_data$lp, myc_data$gf)
colnames(myc_data_df) <- c("species_names", "myc", "lp", "gf")

trait_data_combo <- merge(trait_data, myc_data_df, by = "species_names", all = TRUE)

#figure out which data missing
trait_means$species_names[which(is.na(trait_means$EWT))]

###########For Nmass
#drop pinus rigida
pruned_tree_2 <- drop.tip(pruned_consensus, tip = "Pinus rigida")

trait_means_2 <- trait_means[which(trait_means$species_names %in% pruned_tree_2$tip.label),]

myc_data_df_2 <- myc_data_df[-which(myc_data_df$species_names == "Pinus rigida"),]

trait_data_combo_2 <- merge(trait_means_2, myc_data_df_2, by = "species_names", all = TRUE)

#non-phylo model: linear regression of trait and mycorrhizal type
model_Nmass <- lm(Nmass ~ myc, data = trait_data_combo_2)
Nmass_output <- summary(model_Nmass)

Nmass_aov <- aov(Nmass ~ myc, data = trait_data_combo_2)
summary_aov <- summary(Nmass_aov)

model_Narea <- lm(Narea ~ myc, data = trait_data_combo_2)
Narea_output <- summary(model_Narea)

Nmass_aov <- aov(Nmass ~ myc, data = trait_data_combo_2)
summary(Nmass_aov)

#phylo model
length(pruned_tree_2$tip.label)
vcv_matrix <- vcv.phylo(pruned_tree_2, cor=TRUE)
colnames(trait_data_combo_2)[1] <- "Species"

phylo_Nmass <- gls(Nmass ~ myc, data=trait_data_combo_2, correlation=corPagel(1, pruned_tree_2, form = ~Species))

phylo_Narea <- gls(Narea ~ myc, data=trait_data_combo_2, correlation=corPagel(1, pruned_tree_2, form = ~Species))

plot(phylo_Nmass)    # residual plot
summary(phylo_Nmass) # coefficients, SE's, fit statistics
confint(phylo_Nmass) # confidence intervals for coefficients

plot(phylo_Narea)    # residual plot
summary(phylo_Narea) # coefficients, SE's, fit statistics
confint(phylo_Narea) # confidence intervals for coefficients

###########Lignin mass
#drop Alnus incana
pruned_tree_3 <- drop.tip(pruned_consensus, tip = "Alnus incana")

trait_means_3 <- trait_means[which(trait_means$species_names %in% pruned_tree_3$tip.label),]

myc_data_df_3 <- myc_data_df[-which(myc_data_df$species_names == "Alnus incana"),]

trait_data_combo_3 <- merge(trait_means_3, myc_data_df_3, by = "species_names", all = TRUE)

#non-phylo model: linear regression of trait and mycorrhizal type
model_ligninmass <- lm(lignin_mass ~ myc, data = trait_data_combo_3) 
ligninmass_output <- summary(model_ligninmass) #sig for both

model_ligninarea <- lm(lignin_area ~ myc, data = trait_data_combo_3) 
ligninarea_output <- summary(model_ligninarea) #sig for both

#phylo model
length(pruned_tree_3$tip.label)
vcv_matrix <- vcv.phylo(pruned_tree_3, cor=TRUE)

colnames(trait_data_combo_3)[1] <- "Species"

phylo_ligninmass <- gls(lignin_mass ~ myc, data=trait_data_combo_3, correlation=corPagel(1, pruned_tree_3, form = ~Species))

plot(phylo_ligninmass)    # residual plot
summary(phylo_ligninmass) # coefficients, SE's, fit statistics
confint(phylo_ligninmass) # confidence intervals for coefficients

phylo_ligninarea <- gls(lignin_area ~ myc, data=trait_data_combo_3, correlation=corPagel(1, pruned_tree_3, form = ~Species))

plot(phylo_ligninarea)    # residual plot
summary(phylo_ligninarea) # coefficients, SE's, fit statistics
confint(phylo_ligninarea) # confidence intervals for coefficients


#for LMA and EWT, all taxa

trait_means_all <- trait_means[which(trait_means$species_names %in% pruned_consensus$tip.label),]
trait_data_combo_all <- merge(trait_means_all, myc_data_df, by = "species_names", all = TRUE)

#LMA
#non-phylo model: linear regression of trait and mycorrhizal type
model_LMA <- lm(LMA ~ myc, data = trait_data_combo_all) 
LMA_output <- summary(model_LMA) #sig for both
plot(model_LMA)
#phylo model

LMA_aov <- aov(LMA ~ myc, data = trait_data_combo_all)
summary_aov <- summary(LMA_aov)


vcv_matrix <- vcv.phylo(pruned_consensus, cor=TRUE)

colnames(trait_data_combo_all)[1] <- "Species"

phylo_LMA <- gls(LMA ~ myc, data=trait_data_combo_all, correlation=corPagel(1, pruned_consensus, form = ~Species))


plot(phylo_LMA)    # residual plot
summary(phylo_LMA) # coefficients, SE's, fit statistics
confint(phylo_LMA) # confidence intervals for coefficients

#EWT
model_EWT <- lm(EWT ~ myc, data = trait_data_combo_all) 
EWT_output <- summary(model_EWT) #sig for both

EWT_aov <- aov(EWT ~ myc, data = trait_data_combo_all)
summary_aov <- summary(EWT_aov)


#phylo model
vcv_matrix <- vcv.phylo(pruned_consensus, cor=TRUE)

colnames(trait_data_combo_all)[1] <- "Species"

phylo_EWT <- gls(EWT ~ myc, data=trait_data_combo_all, correlation=corPagel(1, pruned_consensus, form = ~Species))


plot(phylo_EWT)    # residual plot
summary(phylo_EWT) # coefficients, SE's, fit statistics
confint(phylo_EWT) # confidence intervals for coefficients

#############pigments
colnames(metadata)
trait_data_pig <- metadata[,c(28, 86, 107, 108, 109, 135, 137, 139)]
str(trait_data_pig)

#convert to numeric from factor via character
trait_data_pig[,c(3:8)] <- lapply(trait_data_pig[,c(3:8)], function(x) as.numeric(as.character(x)))

#get mean of spectra
trait_means_pig <- trait_data_pig %>% 
  dplyr::select(-c(sample_name)) %>% 
  group_by(species_names) %>% 
  summarize_all(mean, na.rm = TRUE)

#add myc state to trait data
myc_data_df <- as.data.frame(x = myc_data$myc)
myc_data_df <- cbind(myc_data$species, myc_data_df, myc_data$lp, myc_data$gf)
colnames(myc_data_df) <- c("species_names", "myc", "lp", "gf")

#trait_data_combo_pig <- merge(trait_data_pig, myc_data_df, by = "species_names", all = TRUE)

#figure out which data missing
trait_means_pig$species_names[which(is.na(trait_means_pig$chlA_area))]

#Rhus typhina missing from pigments
pruned_tree_pig <- drop.tip(pruned_consensus, tip = "Rhus typhina")

trait_means_pig <- trait_means_pig[which(trait_means_pig$species_names %in% pruned_tree_pig$tip.label),]

myc_data_df_pig <- myc_data_df[-which(myc_data_df$species_names == "Rhus typhina"),]

trait_data_combo_pig <- merge(trait_means_pig, myc_data_df_pig, by = "species_names", all = TRUE)

#non-phylo model: linear regression of trait and mycorrhizal type
model_chla_mass <- lm(chlA_mass ~ myc, data = trait_data_combo_pig) 
chla_mass_output <- summary(model_chla_mass) #sig for both

model_chla_area <- lm(chlA_area ~ myc, data = trait_data_combo_pig) 
chla_area_output <- summary(model_chla_area) #sig for both

model_chlb_mass <- lm(chlB_mass ~ myc, data = trait_data_combo_pig) 
chlb_mass_output <- summary(model_chlb_mass) #sig for both

model_chlb_area <- lm(chlB_area ~ myc, data = trait_data_combo_pig) 
chlb_area_output <- summary(model_chlb_area) #sig for both

model_car_mass <- lm(car_mass ~ myc, data = trait_data_combo_pig) 
car_mass_output <- summary(model_car_mass) #sig for both

model_car_area <- lm(car_area ~ myc, data = trait_data_combo_pig) 
car_area_output <- summary(model_car_area) #sig for both

#summary: pigments by mass significant, pigments by area not (nonphylo)

#phylo model
length(pruned_tree_pig$tip.label)
vcv_matrix <- vcv.phylo(pruned_tree_pig, cor=TRUE)
colnames(trait_data_combo_pig)[1] <- "Species"

phylo_chla_mass <- gls(chlA_mass ~ myc, data=trait_data_combo_pig, correlation=corPagel(1, pruned_tree_pig, form = ~Species))

phylo_chla_area <- gls(chlA_area ~ myc, data=trait_data_combo_pig, correlation=corPagel(1, pruned_tree_pig, form = ~Species))

phylo_chlb_mass <- gls(chlB_mass ~ myc, data=trait_data_combo_pig, correlation=corPagel(1, pruned_tree_pig, form = ~Species))

phylo_chlb_area <- gls(chlB_area ~ myc, data=trait_data_combo_pig, correlation=corPagel(1, pruned_tree_pig, form = ~Species))

phylo_car_mass <- gls(car_mass ~ myc, data=trait_data_combo_pig, correlation=corPagel(1, pruned_tree_pig, form = ~Species))

phylo_car_area <- gls(car_area ~ myc, data=trait_data_combo_pig, correlation=corPagel(1, pruned_tree_pig, form = ~Species))

summary(phylo_chla_mass)
summary(phylo_chla_area)
summary(phylo_chlb_mass)
summary(phylo_chlb_area)
summary(phylo_car_mass)
summary(phylo_car_area)

#phylo not significant




#####################################
#test correlation of mass vs area and lma

cor(trait_data_combo_2$Nmass, trait_data_combo_2$LMA) #-0.604046
cor(trait_data_combo_2$Narea, trait_data_combo_2$LMA) #0.8631537
