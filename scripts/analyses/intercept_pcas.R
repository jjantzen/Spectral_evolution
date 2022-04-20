#intercept model PCAs
library(mvMORPH)
library(phytools)
library(dplyr)
library(ggplot2)
library(OUwie)
source("./scripts/plot_phylo_pca_edited_function.R")
#load plot_pca_phylo function


#read in data
data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")
tree_myc <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

#Read in model results

#for file in folder
files <- list.files("./analysis/intercept_models/myc_dataset/models/", pattern=NULL, all.files=FALSE, full.names=TRUE)

#choose which models to use 
best_models <- readRDS("./analysis/intercept_models/myc_dataset/summarized_best_models_100trees_mycdataset.rds")

best_models <- best_models[,c(1,17)]

#subset files to match best models only
keep_files <- c()

name_formula <- gsub("_BM_iteration1.rds", "", files[1])

for (i in 1:nrow(best_models)){
  file_name <- paste0(name_formula, "_", best_models$best_model[i], "_iteration", best_models$iteration[i], ".rds")  
  keep_files[i] <- file_name
}

#read in models for keep files
model_list <- list()

for (i in 1:length(keep_files)){
  #read file
  model <- readRDS(keep_files[i])
  model_list[[i]] <- model
}

model_list[[i]]
i
#saveRDS(model_list, "./analysis/pca_analysis/best_intercept_models_for_pcas_92sp.rds")
model_list <- readRDS("./analysis/pca_analysis/best_intercept_models_for_pcas_92sp.rds")

#Map trait to trees - make simmap tree for predictor
myc_named <- setNames(data_spectra$myc, rownames(data_spectra$spectra))

#myc_simmaps <- lapply(tree_myc,findMRCA,tips=taxonomy$Species[which(taxonomy$Subclass == "Pinidae")], type="node") #specify nodes c(165)

simmap_tree <- lapply(tree_myc, make.simmap, x = myc_named, model="SYM", nsim=1)
plot(simmap_tree$STATE_241090000)

#assign states to nodes
for (i in 1:length(simmap_tree)){
  simmap_tree[[i]]$node.label<-getStates(simmap_tree[[i]],"nodes")  
}

#saveRDS(simmap_tree, "./analysis/pca_analysis/myc_simmap_trees.rds")
simmap_tree <- readRDS("./analysis/pca_analysis/myc_simmap_trees.rds")

#Do PCA for each of 100 models

for (i in 1:length(model_list)){ #test with 1 first
  #conduct PCA
  pca_best_model <- mvgls.pca(model_list[[i]], plot = FALSE)  
  
  #calculate variance explained for each pc axis
  var_explained <- 100 * pca_best_model$values / sum(pca_best_model$values)

  #create scree plot
  jpeg(paste0("./output/PCAs/iterations/myc_scree_plot_iteration", i, ".jpg"), res = 300, width = 10, height = 10, units = "in")
  q <- qplot(x = c(1:20), y = var_explained[c(1:20)]) +
    geom_line() +
    xlab("Principal Component") +
    ylab("Variance Explained") +
    ggtitle("Scree Plot") +
    geom_hline(yintercept = 0.1)+
    ylim(0, 100)+
    theme_bw()
  print(q)
  dev.off()

  #plot loadings
  vectors_df <- pca_best_model$vectors %>%
    as.data.frame()

  vectors_df$wavelength <- c(400:2400)
  vectors_for_plotting <- vectors_df %>% tidyr::gather(., pc_axis, value, -wavelength)
  vectors_for_plotting$pc_axis <- gsub("V", "", vectors_for_plotting$pc_axis)
  class(vectors_for_plotting$pc_axis) <- "numeric"

  pc_labs <- paste("PC", c(1:9))
  names(pc_labs) <- c(1:9)

  jpeg(paste0("./output/PCAs/iterations/myc_loadings_top10_iteration", i, ".jpg"), res = 400, width = 15, height = 15, units = "in")
  v <- ggplot(vectors_for_plotting[which(vectors_for_plotting$pc_axis %in% c(1:9)),], aes(x=wavelength,y=value)) +
    geom_bar(stat = "identity")+
    facet_wrap(vars(pc_axis), ncol = 3, labeller = labeller(pc_axis = pc_labs))+# scales = "free"
    xlab("Wavelength (nm)")+
    ylab("")+ #Vector unit
    theme_bw()+
    theme(strip.text.x = element_text(size = 15), axis.title=element_text(size=15))
  print(v)
  dev.off()

  #plot pcas
  col.group <- data_spectra$myc
  col.group <- gsub("AM", "orange", col.group)
  col.group <- gsub("EM", "blue", col.group)

  jpeg(paste0("./output/PCAs/iterations/phylo_pca_pc1_pc2_iteration", i, ".jpg"), res = 600, width = 10, height = 10, units = "in")
  pca12 <- plot_pca_phylo(pca_best_model, model_list[[i]], axes = c(1,2), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
  print(pca12)
  dev.off()

  jpeg(paste0("./output/PCAs/iterations/phylo_pca_pc3_pc4_iteration", i, ".jpg"), res = 600, width = 10, height = 10, units = "in")
  pca34 <- plot_pca_phylo(pca_best_model, model_list[[i]], axes = c(3,4), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
  print(pca34)
  dev.off()

  jpeg(paste0("./output/PCAs/iterations/phylo_pca_pc5_pc6_iteration", i, ".jpg"), res = 600, width = 10, height = 10, units = "in")
  pca56 <- plot_pca_phylo(pca_best_model, model_list[[i]], axes = c(5,6), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
  print(pca56)
  dev.off()

  # output_aics_per <- data.frame(matrix(nrow=10, ncol = 9))
  # colnames(output_aics_per)[1] <- "iteration"
  # class(output_aics_per[,1]) <- "numeric"
  # colnames(output_aics_per)[2] <- "PC_axis"
  # class(output_aics_per[,2]) <- "numeric"
  # colnames(output_aics_per)[3] <- "BM1"
  # class(output_aics_per[,3]) <- "numeric"
  # colnames(output_aics_per)[4] <- "BMS"
  # class(output_aics_per[,4]) <- "numeric"
  # colnames(output_aics_per)[5] <- "OU1"
  # class(output_aics_per[,5]) <- "numeric"
  # colnames(output_aics_per)[6] <- "OUM"
  # class(output_aics_per[,6]) <- "numeric"
  # colnames(output_aics_per)[7] <- "OUMV"
  # class(output_aics_per[,7]) <- "numeric"
  # colnames(output_aics_per)[8] <- "OUMA"
  # class(output_aics_per[,8]) <- "numeric"
  # colnames(output_aics_per)[9] <- "OUMVA"
  # class(output_aics_per[,9]) <- "numeric"
  # 
  # output_aov_per <- data.frame(matrix(nrow=10, ncol = 7))
  # colnames(output_aov_per)[1] <- "iteration"
  # class(output_aov_per[,1]) <- "numeric"
  # colnames(output_aov_per)[2] <- "PC_axis"
  # class(output_aov_per[,2]) <- "numeric"
  # colnames(output_aov_per)[3] <- "F_stat"
  # class(output_aov_per[,3]) <- "numeric"
  # colnames(output_aov_per)[4] <- "p_value"
  # class(output_aov_per[,4]) <- "numeric"
  # colnames(output_aov_per)[5] <- "df"
  # class(output_aov_per[,5]) <- "numeric"
  # colnames(output_aov_per)[6] <- "sum_sq"
  # class(output_aov_per[,6]) <- "numeric"
  # colnames(output_aov_per)[7] <- "mean_sq"
  # class(output_aov_per[,7]) <- "numeric"
  # 
  # 
  # #run ouwie models
  # for (j in 1: 10){
  #   #get data
  #   data <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named, X=as.numeric(pca_best_model$scores[,j]))
  #   
  #   #run models for each pc
  #   BM1_pc <- OUwie(simmap_tree[[i]], data, model = "BM1", simmap.tree = TRUE) #single rate
  #   BMS_pc <- OUwie(simmap_tree[[i]], data, model = "BMS", simmap.tree = TRUE) #multiple rates
  #   OU1_pc <- OUwie(simmap_tree[[i]], data, model = "OU1", simmap.tree = TRUE) #single optimum
  #   OUM_pc <- OUwie(simmap_tree[[i]], data, model = "OUM", simmap.tree = TRUE) #multiple optima
  #   OUMV_pc <- OUwie(simmap_tree[[i]], data, model = "OUMV", simmap.tree = TRUE) #multiple optimum and multiple sigmas
  #   OUMA_pc <- OUwie(simmap_tree[[i]], data, model = "OUMA", simmap.tree = TRUE) #multiple optimum and multiple alphas
  #   OUMVA_pc <- OUwie(simmap_tree[[i]], data, model = "OUMVA", simmap.tree = TRUE) #multiple optimum and multiple sigmas and multiple alphas
  #   print(j)
  #   #save aic.w scores
  #   aics <- setNames(c(BM1_pc$AIC, BMS_pc$AIC,OU1_pc$AIC, OUM_pc$AIC, OUMV_pc$AIC, OUMA_pc$AIC, OUMVA_pc$AIC), c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"))
  #   weights <- aic.w(aics)
  #   
  #   output_aics_per$iteration[j] <- i
  #   output_aics_per$PC_axis[j] <- j
  #   output_aics_per$BM1[j] <- weights[1]
  #   output_aics_per$BMS[j] <- weights[2]
  #   output_aics_per$OU1[j] <- weights[3]
  #   output_aics_per$OUM[j] <- weights[4]
  #   output_aics_per$OUMV[j] <- weights[5]
  #   output_aics_per$OUMA[j] <- weights[6]
  #   output_aics_per$OUMVA[j] <- weights[7]
  # 
  #   #figure out how to save aov results
  #   pc_aov <- aov(X ~ myc, data = data)
  #   summary_aov <- summary(pc_aov)
  #  
  #   output_aov_per$iteration[j] <- i
  #   output_aov_per$PC_axis[j] <- j
  #   output_aov_per$F_stat[j] <- summary_aov[[1]]$`F value`[1]
  #   output_aov_per$p_value[j] <- summary_aov[[1]]$`Pr(>F)`[1]
  #   output_aov_per$df[j] <- summary_aov[[1]]$Df
  #   output_aov_per$sum_sq[j] <- summary_aov[[1]]$`Sum Sq`
  #   output_aov_per$mean_sq[j] <- summary_aov[[1]]$`Mean Sq`
  #   
  # }
  # print(i)
  # if (i == 1){
  #   output_aics <- output_aics_per
  #   output_aovs <- output_aov_per
  #   
  # } else {
  #   output_aics <- rbind(output_aics, output_aics_per)
  #   output_aovs <- rbind(output_aovs, output_aov_per)
  # }
}

saveRDS(output_aics, "./analysis/pca_analysis/pc_univariate_model_aics_92sp.rds")
saveRDS(output_aovs, "./analysis/pca_analysis/pc_univariate_model_aovs_92sp.rds")
output_aovs[c(900:980), c(1:4)]
nrow(output_aovs)
tail(output_aics)
# output_aics <- readRDS("./analysis/pca_analysis/pc_univariate_model_aics_92sp.rds")
# output_aovs <- readRDS("./analysis/pca_analysis/pc_univariate_model_aovs_92sp.rds")




