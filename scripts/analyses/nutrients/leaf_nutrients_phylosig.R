#leaf nutrients phylogenetic signal

summarized_best_models_100trees_mycdataset
colnames(summarized_best_models_100trees_mycdataset)
best_models_myc_models_92sp

library(mvMORPH)
library(phytools)
library(dplyr)
library(spectrolab)
library(OUwie)

#read data
spectral_data <- readRDS("./data/for_analysis/final_spectra_matched_trees.rds")

predictor_data <- readRDS("./data/for_analysis/final_data.rds")

myc_data <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

trees_big <- readRDS("./data/for_analysis/final_trees_matched_spectra.rds")

consensus_tree <- readRDS("./data/for_analysis/consensus_tree_full.rds")

#prep data
metadata <- meta(spectral_data)
colnames(metadata)
trait_data <- metadata[,c(28, 86, 95:164)]
colnames(trait_data)

#convert to numeric from factor via character
trait_data[,c(3:72)] <- lapply(trait_data[,c(3:72)], function(x) as.numeric(as.character(x)))

#get mean of spectra
trait_means <- trait_data[,c(1,3:72)] %>% 
  #dplyr::select(-c(sample_name)) %>% 
  group_by(species_names) %>% 
  summarize_all(mean, na.rm = TRUE)

unique(trait_means$species_names)

#prune trees to match data - not needed yet (unless trim both to myc dataset)
exclude <- consensus_tree$tip.label[-which(consensus_tree$tip.label %in% myc_data$species)]
pruned_consensus <- drop.tip(consensus_tree, tip = exclude)

trees_object <- lapply(trees_big,drop.tip,tip=exclude)

trait_means_keep <- trait_means[-which(trait_means$species_names %in% exclude),]

#function to calculate phylosig over 100 trees
physig_iterations <- function(trait_data, trees_object, name){ #name is name of leaf trait
  df_output <- data.frame(matrix(nrow=100, ncol = 3))
  colnames(df_output)[1] <- "iteration"
  class(df_output[,1]) <- "numeric"
  colnames(df_output)[2] <- "K"
  class(df_output[,2]) <- "numeric"
  colnames(df_output)[3] <- "pvalue"
  class(df_output[,3]) <- "numeric"
  
  #do physignal calc
  for (i in 1:length(trees_object)){
    sig <- phylosig(trees_object[[i]], trait_data, method = "K", test=TRUE, nsim=999) 
    df_output$iteration[i] <- i
    df_output$K[i] <- sig$K
    df_output$pvalue[i] <- sig$P
  }
  #get averages
  # df_output$iteration[101] <- NA
  # df_output$K[101] <-  mean(df_output$K[c(1:100)])
  # df_output$pvalue[101] <-  mean(df_output$pvalue[1:100])
  # df_output$iteration[102] <- NA
  # df_output$K[102] <-  max(df_output$K[c(1:100)])
  # df_output$pvalue[102] <-  max(df_output$pvalue[1:100])
  # 
  saveRDS(df_output, paste0("./analysis/trait_phylosig/", name, "_K_100trees_leaf_nutrients.rds"))
  write.csv(df_output, paste0("./analysis/trait_phylosig/", name, "_K_100trees_leaf_nutrients.csv"))
  return(df_output)
}

nrow(trait_means_keep[which(is.nan(trait_means_keep$P_mass) == TRUE),])

unique(trait_data$species_names[which(!is.na(trait_data$P_area))]) #24 species have phosphorus in original dataset

#missing Nmass and Narea from Pinus rigida

#missing Pmass and Parea from 73 species

#missing lignin and cellulose from Alnus incana

trait_means_keep[which(is.nan(trait_means_keep$P_mass) == TRUE),c(1, 27, 68)]

colnames(trait_means_keep)



#save output
SLA <- setNames(trait_means_keep$SLA, trait_means_keep$species_names)
SLA_K <- physig_iterations(SLA, trees_object, "SLA")
sd(SLA_K$pvalue)
nrow(SLA_K[which((SLA_K$pvalue < 0.05) == TRUE),])

Cmass <- setNames(trait_means_keep$Cmass, trait_means_keep$species_names)
Cmass_K <- physig_iterations(Cmass, trees_object, "Cmass")
sd(Cmass_K$pvalue)
nrow(Cmass_K[which((Cmass_K$pvalue < 0.05) == TRUE),])

Nmass <- setNames(trait_means_keep$Nmass, trait_means_keep$species_names)
Nmass_K <- physig_iterations(Nmass, trees_object, "Nmass")
sd(Nmass_K$pvalue)
nrow(Nmass_K[which((Nmass_K$pvalue < 0.05) == TRUE),])

lignin_mass <- setNames(trait_means_keep$lignin_mass, trait_means_keep$species_names)
lignin_mass_K <- physig_iterations(lignin_mass, trees_object, "lignin_mass")
sd(lignin_mass_K$pvalue)
nrow(lignin_mass_K[which((lignin_mass_K$pvalue < 0.05) == TRUE),])

chlA_mass <- setNames(trait_means_keep$chlA_mass, trait_means_keep$species_names)
chlA_mass_K <- physig_iterations(chlA_mass, trees_object, "chlA_mass")
sd(chlA_mass_K$pvalue)
nrow(chlA_mass_K[which((chlA_mass_K$pvalue < 0.05) == TRUE),])

chlB_mass <- setNames(trait_means_keep$chlB_mass, trait_means_keep$species_names)
chlB_mass_K <- physig_iterations(chlB_mass, trees_object, "chlB_mass")
sd(chlB_mass_K$pvalue)
nrow(chlB_mass_K[which((chlB_mass_K$pvalue < 0.05) == TRUE),])

car_mass <- setNames(trait_means_keep$car_mass, trait_means_keep$species_names)
car_mass_K <- physig_iterations(car_mass, trees_object, "car_mass")
sd(car_mass_K$pvalue)
nrow(car_mass_K[which((car_mass_K$pvalue < 0.05) == TRUE),])

P_mass <- setNames(trait_means_keep$P_mass, trait_means_keep$species_names)
P_mass_K <- physig_iterations(P_mass, trees_object, "P_mass")
sd(P_mass_K$pvalue)
nrow(P_mass_K[which((P_mass_K$pvalue < 0.05) == TRUE),])

Cnorm <- setNames(trait_means_keep$Cnorm, trait_means_keep$species_names)
Cnorm_K <- physig_iterations(Cnorm, trees_object, "Cnorm")
sd(Cnorm_K$pvalue)
nrow(Cnorm_K[which((Cnorm_K$pvalue < 0.05) == TRUE),])

Nnorm <- setNames(trait_means_keep$Nnorm, trait_means_keep$species_names)
Nnorm_K <- physig_iterations(Nnorm, trees_object, "Nnorm")
sd(Nnorm_K$K)
mean(Nnorm_K$K)
sd(Nnorm_K$pvalue)
mean(Nnorm_K$pvalue)
nrow(Nnorm_K[which((Nnorm_K$pvalue < 0.05) == TRUE),])

lignin_norm <- setNames(trait_means_keep$lignin_norm, trait_means_keep$species_names)
lignin_norm_K <- physig_iterations(lignin_norm, trees_object, "lignin_norm")
sd(lignin_norm_K$K)
mean(lignin_norm_K$K)
sd(lignin_norm_K$pvalue)
mean(lignin_norm_K$pvalue)
nrow(lignin_norm_K[which((lignin_norm_K$pvalue < 0.05) == TRUE),])

P_norm <- setNames(trait_means_keep$P_norm, trait_means_keep$species_names)
P_norm_K <- physig_iterations(P_norm, trees_object, "P_norm")
sd(P_norm_K$K)
mean(P_norm_K$K)
sd(P_norm_K$pvalue)
mean(P_norm_K$pvalue)
nrow(P_norm_K[which((P_norm_K$pvalue < 0.05) == TRUE),])

chlA_norm <- setNames(trait_means_keep$chlA_norm, trait_means_keep$species_names)
chlA_norm_K <- physig_iterations(chlA_norm, trees_object, "chlA_norm")
mean(chlA_norm_K$K)
sd(chlA_norm_K$K)
mean(chlA_norm_K$pvalue)
sd(chlA_norm_K$pvalue)
nrow(chlA_norm_K[which((chlA_norm_K$pvalue < 0.05) == TRUE),])

car_norm <- setNames(trait_means_keep$car_norm, trait_means_keep$species_names)
car_norm_K <- physig_iterations(car_norm, trees_object, "car_norm")
mean(car_norm_K$K)
sd(car_norm_K$K)
mean(car_norm_K$pvalue)
sd(car_norm_K$pvalue)
nrow(car_norm_K[which((car_norm_K$pvalue < 0.05) == TRUE),])

chlB_norm <- setNames(trait_means_keep$chlB_norm, trait_means_keep$species_names)
chlB_norm_K <- physig_iterations(chlB_norm, trees_object, "chlB_norm")
mean(chlB_norm_K$K)
sd(chlB_norm_K$K)
mean(chlB_norm_K$pvalue)
sd(chlB_norm_K$pvalue)
nrow(chlB_norm_K[which((chlB_norm_K$pvalue < 0.05) == TRUE),])

hemicellulose_mass <- setNames(trait_means_keep$hemicellulose_mass, trait_means_keep$species_names)
hemicellulose_mass_K <- physig_iterations(hemicellulose_mass, trees_object, "hemicellulose_mass")
mean(hemicellulose_mass_K$K)
sd(hemicellulose_mass_K$K)
mean(hemicellulose_mass_K$pvalue)
sd(hemicellulose_mass_K$pvalue)
nrow(hemicellulose_mass_K[which((hemicellulose_mass_K$pvalue < 0.05) == TRUE),])

hemicellulose_norm <- setNames(trait_means_keep$hemicellulose_norm, trait_means_keep$species_names)
hemicellulose_norm_K <- physig_iterations(hemicellulose_norm, trees_object, "hemicellulose_norm")
mean(hemicellulose_norm_K$K)
sd(hemicellulose_norm_K$K)
mean(hemicellulose_norm_K$pvalue)
sd(hemicellulose_norm_K$pvalue)
nrow(hemicellulose_norm_K[which((hemicellulose_norm_K$pvalue < 0.05) == TRUE),])

LMA <- setNames(trait_means_keep$LMA, trait_means_keep$species_names)
LMA_K <- physig_iterations(LMA, trees_object, "LMA")
mean(LMA_K$K)
sd(LMA_K$K)
mean(LMA_K$pvalue)
sd(LMA_K$pvalue)
nrow(LMA_K[which((LMA_K$pvalue < 0.05) == TRUE),])


#plot over tree
jpeg("./output/trait_plots/Leaf_nutrients_phenogram_LMA.jpg", width = 8, height = 12, res = 400, units = "in")
obj<-contMap(pruned_consensus,LMA,plot=FALSE)
plot(obj,lwd=7)#,xlim=c(-0.2,3.6)
#phenogram(pruned_consensus,trait_means_keep$LMA,spread.labels=TRUE,spread.cost=c(10,0), fsize = 0.7)
dev.off()

jpeg("./output/trait_plots/Leaf_nutrients_phenogram_Pnorm.jpg", width = 8, height = 12, res = 400, units = "in")
obj<-contMap(pruned_consensus,P_norm,plot=FALSE)
plot(obj,lwd=7)#,xlim=c(-0.2,3.6)
#phenogram(pruned_consensus,trait_means_keep$LMA,spread.labels=TRUE,spread.cost=c(10,0), fsize = 0.7)
dev.off()

jpeg("./output/trait_plots/Leaf_nutrients_phenogram_Pmass.jpg", width = 8, height = 12, res = 400, units = "in")
obj<-contMap(pruned_consensus,P_mass,plot=FALSE)
plot(obj,lwd=7)#,xlim=c(-0.2,3.6)
#phenogram(pruned_consensus,trait_means_keep$LMA,spread.labels=TRUE,spread.cost=c(10,0), fsize = 0.7)
dev.off()

jpeg("./output/trait_plots/Leaf_nutrients_phenogram_hemicellulose_norm.jpg", width = 8, height = 12, res = 400, units = "in")
obj<-contMap(pruned_consensus,hemicellulose_norm,plot=FALSE)
plot(obj,lwd=7)#,xlim=c(-0.2,3.6)
#phenogram(pruned_consensus,trait_means_keep$LMA,spread.labels=TRUE,spread.cost=c(10,0), fsize = 0.7)
dev.off()

jpeg("./output/trait_plots/Leaf_nutrients_phenogram_cnorm.jpg", width = 8, height = 12, res = 400, units = "in")
obj<-contMap(pruned_consensus,Cnorm,plot=FALSE)
plot(obj,lwd=7)#,xlim=c(-0.2,3.6)
#phenogram(pruned_consensus,trait_means_keep$LMA,spread.labels=TRUE,spread.cost=c(10,0), fsize = 0.7)
dev.off()

jpeg("./output/trait_plots/Leaf_nutrients_phenogram_SLA.jpg", width = 8, height = 12, res = 400, units = "in")
obj<-contMap(pruned_consensus,SLA,plot=FALSE)
plot(obj,lwd=7)#,xlim=c(-0.2,3.6)
#phenogram(pruned_consensus,trait_means_keep$LMA,spread.labels=TRUE,spread.cost=c(10,0), fsize = 0.7)
dev.off()

jpeg("./output/trait_plots/Leaf_nutrients_phenogram_car_norm.jpg", width = 8, height = 12, res = 400, units = "in")
obj<-contMap(pruned_consensus,car_norm,plot=FALSE)
plot(obj,lwd=7)#,xlim=c(-0.2,3.6)
#phenogram(pruned_consensus,trait_means_keep$LMA,spread.labels=TRUE,spread.cost=c(10,0), fsize = 0.7)
dev.off()
car_norm


#do models 

trees_object 

trait_means_keep 

data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

myc_named <- setNames(data_spectra$myc, rownames(data_spectra$spectra))

simmap_tree <- readRDS("./analysis/pca_analysis/myc_simmap_trees.rds")

trait_means_modeling <- trait_means_keep[,c(1,2,4:6,13:16,29:30,32,40,42,44,46)]

trait_means_modeling <- as.data.frame(trait_means_modeling)
colnames(trait_means_modeling)
#exclude Narea?
trait_means_modeling <- trait_means_modeling[,-c(12,14)]

for (i in 1:length(simmap_tree)){ #test with 1 first
  
  output_aics_per <- data.frame(matrix(nrow=9, ncol = 9))
  colnames(output_aics_per)[1] <- "iteration"
  class(output_aics_per[,1]) <- "numeric"
  colnames(output_aics_per)[2] <- "trait"
  class(output_aics_per[,2]) <- "character"
  colnames(output_aics_per)[3] <- "BM1"
  class(output_aics_per[,3]) <- "numeric"
  colnames(output_aics_per)[4] <- "BMS"
  class(output_aics_per[,4]) <- "numeric"
  colnames(output_aics_per)[5] <- "OU1"
  class(output_aics_per[,5]) <- "numeric"
  colnames(output_aics_per)[6] <- "OUM"
  class(output_aics_per[,6]) <- "numeric"
  colnames(output_aics_per)[7] <- "OUMV"
  class(output_aics_per[,7]) <- "numeric"
  colnames(output_aics_per)[8] <- "OUMA"
  class(output_aics_per[,8]) <- "numeric"
  colnames(output_aics_per)[9] <- "OUMVA"
  class(output_aics_per[,9]) <- "numeric"
  
  #run ouwie models
  for (j in 2:10){
    #get data
    data <- data.frame(species=trait_means_modeling$species_names, myc=myc_named, X=as.numeric(trait_means_modeling[,j]))
    
    trimmed_tree <- drop.tip(simmap_tree[[i]], tip = trait_means_modeling$species_names[which(is.nan(trait_means_modeling[,j]))])
    i 
  
    #run models for each pc
    BM1_pc <- OUwie(trimmed_tree, data, model = "BM1", simmap.tree = TRUE) #single rate
    BMS_pc <- OUwie(trimmed_tree, data, model = "BMS", simmap.tree = TRUE) #multiple rates
    OU1_pc <- OUwie(trimmed_tree, data, model = "OU1", simmap.tree = TRUE) #single optimum
    OUM_pc <- OUwie(trimmed_tree, data, model = "OUM", simmap.tree = TRUE) #multiple optima
    OUMV_pc <- OUwie(trimmed_tree, data, model = "OUMV", simmap.tree = TRUE) #multiple optimum and multiple sigmas
    OUMA_pc <- OUwie(trimmed_tree, data, model = "OUMA", simmap.tree = TRUE) #multiple optimum and multiple alphas
    OUMVA_pc <- OUwie(trimmed_tree, data, model = "OUMVA", simmap.tree = TRUE) #multiple optimum and multiple sigmas and multiple alphas
    print(colnames(trait_means_modeling[j]))
    #save aic.w scores
    aics <- setNames(c(BM1_pc$AIC, BMS_pc$AIC,OU1_pc$AIC, OUM_pc$AIC, OUMV_pc$AIC, OUMA_pc$AIC, OUMVA_pc$AIC), c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"))
    weights <- aic.w(aics)
    
    output_aics_per$iteration[j-1] <- i
    output_aics_per$trait[j-1] <- colnames(trait_means_modeling[j])
    output_aics_per$BM1[j-1] <- weights[1]
    output_aics_per$BMS[j-1] <- weights[2]
    output_aics_per$OU1[j-1] <- weights[3]
    output_aics_per$OUM[j-1] <- weights[4]
    output_aics_per$OUMV[j-1] <- weights[5]
    output_aics_per$OUMA[j-1] <- weights[6]
    output_aics_per$OUMVA[j-1] <- weights[7]
    
  }
  print(i)
  if (i == 1){
    output_aics <- output_aics_per
    
  } else {
    output_aics <- rbind(output_aics, output_aics_per)
  }
}

saveRDS(output_aics, "./analysis/leaf_trait_models/leaf_trait_univariate_model_aics_92sp.rds")

#summarize output
best_ouwie_models <- output_aics %>%
  rowwise() %>%
  mutate(top_model = names(.)[which.max(c_across(BM1:OUMVA))+2], second_model = tail(head(names(cur_data())[order(c_across(BM1:OUMVA), decreasing = T)+2],2),1))#(.[3:9], 1, function(x) names(x)[maxn(2)(x)]))

best_ouwie_models

summarized_best_models <- best_ouwie_models %>% 
  dplyr::select(iteration, trait, top_model) %>% 
  dplyr::group_by(trait, top_model) %>% 
  dplyr::summarize(count = n())

summarized_best_models_df <- as.data.frame(summarized_best_models)

#plot the aic weights for each model for each pc axis (repeated iterations)
long_aics <- pivot_longer(output_aics, c(3:9))
long_aics

jpeg("./output/trait_plots/ouwie_aicw_boxplot_facet_by_trait_col.jpg", width = 12, height = 8, units = "in", res = 700)
ggplot(long_aics, aes(x = name, y = round(value, 2)))+
  #geom_point()+
  geom_boxplot(aes(fill = as.factor(name)))+
  #geom_jitter(color = "black", size = 0.5)+
  labs(y = "AICw", x = "Model")+
  facet_wrap(~as.factor(trait), ncol = 3)+
  scale_fill_discrete(name = "Model")+
  theme_bw()
dev.off()
