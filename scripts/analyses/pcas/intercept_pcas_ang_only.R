#intercept model PCAs
library(mvMORPH)
library(phytools)
library(dplyr)
library(ggplot2)
library(OUwie)
library(phylolm)
source("./scripts/plot_phylo_pca_edited_function.R")
#load plot_pca_phylo function


#read in data
data_spectra <- readRDS("./data/for_analysis/ang_only_data_for_myc.rds")
tree_myc <- readRDS("./data/for_analysis/ang_only_trees_for_myc.rds")

#Read in model results

#for file in folder
files <- list.files("./analysis/intercept_models/myc_dataset/ang_only/models/", pattern=NULL, all.files=FALSE, full.names=TRUE)

#choose which models to use 
best_models <- readRDS("./analysis/intercept_models/myc_dataset/summarized_best_models_100trees_mycdataset.rds")

best_models <- best_models[,c(1,17)]

#subset files to match best models only
keep_files <- c()

#did this on cluster
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
model_list <- readRDS("./analysis/intercept_models/myc_dataset/ang_only/best_intercept_models_ang_only.rds")

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


#drop tips not in dataset
tree_drops <- simmap_tree[[1]]$tip.label[-which(simmap_tree[[1]]$tip.label %in% data_spectra$species)]

simmap_tree <- lapply(simmap_tree,drop.tip,tip=tree_drops)


#Do PCA for each of 100 models

for (i in 1:length(model_list)){ #test with 1 first
  #conduct PCA
  pca_best_model <- mvgls.pca(model_list[[i]], plot = FALSE)  
  
  #calculate variance explained for each pc axis
  var_explained <- 100 * pca_best_model$values / sum(pca_best_model$values)

  #create scree plot
  jpeg(paste0("./output/PCAs/iterations/myc_scree_plot_ang_only_iteration_", i, ".jpg"), res = 300, width = 10, height = 10, units = "in")
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

  jpeg(paste0("./output/PCAs/iterations/myc_loadings_top10_ang_only_iteration", i, ".jpg"), res = 400, width = 15, height = 15, units = "in")
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

  jpeg(paste0("./output/PCAs/iterations/phylo_pca_pc1_pc2_ang_only_iteration", i, ".jpg"), res = 600, width = 10, height = 10, units = "in")
  pca12 <- plot_pca_phylo(pca_best_model, model_list[[i]], axes = c(1,2), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
  print(pca12)
  dev.off()

  jpeg(paste0("./output/PCAs/iterations/phylo_pca_pc3_pc4_ang_only_iteration", i, ".jpg"), res = 600, width = 10, height = 10, units = "in")
  pca34 <- plot_pca_phylo(pca_best_model, model_list[[i]], axes = c(3,4), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
  print(pca34)
  dev.off()

  jpeg(paste0("./output/PCAs/iterations/phylo_pca_pc5_pc6_ang_only_iteration", i, ".jpg"), res = 600, width = 10, height = 10, units = "in")
  pca56 <- plot_pca_phylo(pca_best_model, model_list[[i]], axes = c(5,6), col = col.group, groups = data_spectra$myc, labels = rownames(data_spectra$spectra))
  print(pca56)
  dev.off()

  output_aics_per <- data.frame(matrix(nrow=10, ncol = 9))
  colnames(output_aics_per)[1] <- "iteration"
  class(output_aics_per[,1]) <- "numeric"
  colnames(output_aics_per)[2] <- "PC_axis"
  class(output_aics_per[,2]) <- "numeric"
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


  #run ouwie models
  for (j in 1: 10){
    #get data
    data <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named, X=as.numeric(pca_best_model$scores[,j]))

    #run models for each pc
    BM1_pc <- OUwie(simmap_tree[[i]], data, model = "BM1", simmap.tree = TRUE) #single rate
    BMS_pc <- OUwie(simmap_tree[[i]], data, model = "BMS", simmap.tree = TRUE) #multiple rates
    OU1_pc <- OUwie(simmap_tree[[i]], data, model = "OU1", simmap.tree = TRUE) #single optimum
    OUM_pc <- OUwie(simmap_tree[[i]], data, model = "OUM", simmap.tree = TRUE) #multiple optima
    OUMV_pc <- OUwie(simmap_tree[[i]], data, model = "OUMV", simmap.tree = TRUE) #multiple optimum and multiple sigmas
    OUMA_pc <- OUwie(simmap_tree[[i]], data, model = "OUMA", simmap.tree = TRUE) #multiple optimum and multiple alphas
    OUMVA_pc <- OUwie(simmap_tree[[i]], data, model = "OUMVA", simmap.tree = TRUE) #multiple optimum and multiple sigmas and multiple alphas
    print(j)
    #save aic.w scores
    aics <- setNames(c(BM1_pc$AIC, BMS_pc$AIC,OU1_pc$AIC, OUM_pc$AIC, OUMV_pc$AIC, OUMA_pc$AIC, OUMVA_pc$AIC), c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"))
    weights <- aic.w(aics)

    output_aics_per$iteration[j] <- i
    output_aics_per$PC_axis[j] <- j
    output_aics_per$BM1[j] <- weights[1]
    output_aics_per$BMS[j] <- weights[2]
    output_aics_per$OU1[j] <- weights[3]
    output_aics_per$OUM[j] <- weights[4]
    output_aics_per$OUMV[j] <- weights[5]
    output_aics_per$OUMA[j] <- weights[6]
    output_aics_per$OUMVA[j] <- weights[7]

    # #figure out how to save aov results
    # pc_aov <- aov(X ~ myc, data = data)
    # summary_aov <- summary(pc_aov)
    # 
    # output_aov_per$iteration[j] <- i
    # output_aov_per$PC_axis[j] <- j
    # output_aov_per$F_stat[j] <- summary_aov[[1]]$`F value`[1]
    # output_aov_per$p_value[j] <- summary_aov[[1]]$`Pr(>F)`[1]
    # output_aov_per$df[j] <- summary_aov[[1]]$Df
    # output_aov_per$sum_sq[j] <- summary_aov[[1]]$`Sum Sq`
    # output_aov_per$mean_sq[j] <- summary_aov[[1]]$`Mean Sq`
    # 
    
  }
  print(i)
  if (i == 1){
    output_aics <- output_aics_per
   # output_aovs <- output_aov_per

  } else {
    output_aics <- rbind(output_aics, output_aics_per)
   # output_aovs <- rbind(output_aovs, output_aov_per)
  }
}

saveRDS(output_aics, "./analysis/pca_analysis/pc_univariate_model_aics_ang_only.rds")
#saveRDS(output_aovs, "./analysis/pca_analysis/pc_univariate_model_aovs_92sp.rds")
output_aovs[c(900:980), c(1:4)]
nrow(output_aovs)
tail(output_aics)
# output_aics <- readRDS("./analysis/pca_analysis/pc_univariate_model_aics_92sp.rds")
# output_aovs <- readRDS("./analysis/pca_analysis/pc_univariate_model_aovs_92sp.rds")


for (i in 1:length(model_list)){ #test with 1 first
  #conduct PCA
  pca_best_model <- mvgls.pca(model_list[[i]], plot = FALSE)  
  
  output_aics_per <- data.frame(matrix(nrow=100, ncol = 5))
  colnames(output_aics_per)[1] <- "iteration"
  class(output_aics_per[,1]) <- "numeric"
  colnames(output_aics_per)[2] <- "BM1"
  class(output_aics_per[,2]) <- "numeric"
  colnames(output_aics_per)[3] <- "BMM"
  class(output_aics_per[,3]) <- "numeric"
  colnames(output_aics_per)[4] <- "OU1"
  class(output_aics_per[,4]) <- "numeric"
  colnames(output_aics_per)[5] <- "OUM"
  class(output_aics_per[,5]) <- "numeric"
  

  #run mvmorph models for pcs 1:6
  BM1 <- mvBM(simmap_tree[[i]], pca_best_model$scores[,1:6], model = "BM1") #single rate
  BMM <- mvBM(simmap_tree[[i]], pca_best_model$scores[,1:6], model = "BMM") #multiple rates
  OU1 <- mvOU(simmap_tree[[i]], pca_best_model$scores[,1:6], model = "OU1") #single optimum
  OUM <- mvOU(simmap_tree[[i]], pca_best_model$scores[,1:6], model = "OUM") #multiple optima
  
  #save aic.w scores
  aics <- setNames(c(BM1$AICc, BMM$AICc,OU1$AICc, OUM$AICc), c("BM1", "BMM", "OU1", "OUM"))
  weights <- aic.w(aics)

  output_aics_per$iteration[i] <- i
  output_aics_per$BM1[i] <- weights[1]
  output_aics_per$BMM[i] <- weights[2]
  output_aics_per$OU1[i] <- weights[3]
  output_aics_per$OUM[i] <- weights[4]
  
  print(i)
  if (i == 1){
    output_aics <- output_aics_per
    
  } else {
    output_aics <- rbind(output_aics, output_aics_per)
  }
}

df_sm <- pc_univariate_model_aics_92sp[which(pc_univariate_model_aics_92sp$PC_axis %in% c(1:6) & pc_univariate_model_aics_92sp$iteration == 1),c(1,5,6,7,8,9)]

colSums(df_sm[,-1])

#######do phylo analysis on pcs
for (i in 1:length(model_list)){ #test with 1 first
  #conduct PCA
  pca_best_model <- mvgls.pca(model_list[[i]], plot = FALSE)  
  
  output_gls_per <- data.frame(matrix(nrow=10, ncol = 6))
  colnames(output_gls_per)[1] <- "iteration"
  class(output_gls_per[,1]) <- "numeric"
  colnames(output_gls_per)[2] <- "PC_axis"
  class(output_gls_per[,2]) <- "numeric"
  colnames(output_gls_per)[3] <- "AIC"
  class(output_gls_per[,3]) <- "numeric"
  colnames(output_gls_per)[4] <- "p_value"
  class(output_gls_per[,4]) <- "numeric"
  # colnames(output_gls_per)[5] <- "lambda"
  # class(output_gls_per[,5]) <- "numeric"
  colnames(output_gls_per)[5] <- "parameter"
  class(output_gls_per[,5]) <- "numeric"
  colnames(output_gls_per)[6] <- "coefficient"
  class(output_gls_per[,6]) <- "numeric"
  
  output_gls_best <- data.frame(matrix(nrow=10, ncol = 6))
  colnames(output_gls_best)[1] <- "iteration"
  class(output_gls_best[,1]) <- "numeric"
  colnames(output_gls_best)[2] <- "PC_axis"
  class(output_gls_best[,2]) <- "numeric"
  colnames(output_gls_best)[3] <- "AIC"
  class(output_gls_best[,3]) <- "numeric"
  colnames(output_gls_best)[4] <- "p_value"
  class(output_gls_best[,4]) <- "numeric"
  # colnames(output_gls_per)[5] <- "lambda"
  # class(output_gls_per[,5]) <- "numeric"
  colnames(output_gls_best)[5] <- "parameter"
  class(output_gls_best[,5]) <- "numeric"
  colnames(output_gls_best)[6] <- "coefficient"
  class(output_gls_best[,6]) <- "numeric"
  
  #run models
  for (j in 1:6){ #10
    #get data
    data <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named, X=as.numeric(pca_best_model$scores[,j]))
  
    #do phylo model for each pc
    #phy_gls_lam <- gls(X~myc,data=data, correlation=corPagel(1,simmap_tree[[i]], form = ~species)) #model_list[[i]]$param
    initial_model <- model_list[[i]]$model
    if (initial_model == "OU"){
      initial_model <- "OUrandomRoot"
    }
    phy_gls_lam <- phylolm(X~myc,data=data, phy = simmap_tree[[i]], model = initial_model) #model_list[[i]]$param
    summary_gls <- summary(phy_gls_lam)
   
    # output_gls_per$iteration[j] <- i
    # output_gls_per$PC_axis[j] <- j
    # output_gls_per$AIC[j] <- summary_gls$AIC
    # output_gls_per$p_value[j] <- summary_gls$tTable[[8]]
    # output_gls_per$coefficient[j] <- summary_gls$tTable[[2]]
    # output_gls_per$lambda[j] <- summary_gls[[1]]$corStruct[[1]]
    
    output_gls_per$iteration[j] <- i
    output_gls_per$PC_axis[j] <- j
    output_gls_per$AIC[j] <- summary_gls$aic
    output_gls_per$p_value[j] <- summary_gls$coefficients[[8]]
    output_gls_per$coefficient[j] <- summary_gls$coefficients[[2]]
    if (is.null(summary_gls$optpar) == FALSE){
      output_gls_per$parameter[j] <- summary_gls$optpar  
    } else {
      output_gls_per$parameter[j] <- NA
    }
    output_gls_per$model[j] <- model_list[[i]]$model
    
    #testing model fit
    phy_gls_BM <- phylolm(X~myc,data=data, phy = simmap_tree[[i]], model = "BM") 
    phy_gls_lambda <- phylolm(X~myc,data=data, phy = simmap_tree[[i]], model = "lambda") 
    phy_gls_EB <- phylolm(X~myc,data=data, phy = simmap_tree[[i]], model = "EB") 
    phy_gls_OUrr <- phylolm(X~myc,data=data, phy = simmap_tree[[i]], model = "OUrandomRoot") 
    
    #compare aics
    aics_phylolm <- setNames(c(phy_gls_BM$aic, phy_gls_lambda$aic, phy_gls_EB$aic, phy_gls_OUrr$aic), c("BM", "lambda", "EB", "OUrr"))
    weights <- aic.w(aics_phylolm)
    
    #get output for "best" model via aic
    best_model <- names(weights)[which(weights == max(weights))]
    name <- paste0("phy_gls_", best_model)
    
    output_gls_best$iteration[j] <- i
    output_gls_best$PC_axis[j] <- j
    output_gls_best$AIC[j] <- summary(get(paste0("phy_gls_", best_model)))$aic
    output_gls_best$p_value[j] <- summary(get(paste0("phy_gls_", best_model)))$coefficients[[8]]
    output_gls_best$coefficient[j] <- summary(get(paste0("phy_gls_", best_model)))$coefficients[[2]]
    output_gls_best$parameter[j] <- summary(get(paste0("phy_gls_", best_model)))$optpar
    output_gls_best$model[j] <- best_model
    
    
    ###plot model assessment for each pc for each iteration
    #jpeg(paste0("./output/PCAs/iterations/model_assessment_pc", j, "_iteration", i, "_summary_phylolm", model_list[[i]]$model, ".jpg"), width = 8, height = 6, res = 500, units = "in")
    #par(mfrow=c(2,2))
    
    #plot residuals vs independent variable
    #print(plot(data$X, phy_gls_lam$residuals, ylab = "Residuals", xlab = "Variable"))
    
    #Distribution of residuals (normal)
    #qplot <- qqnorm(phy_gls_lam$residuals)
    #print(qplot)
    #print(qqline(phy_gls_lam$residuals, col = "red"))
    
    #density plot
    #print(plot(density(phy_gls_lam$residuals), main = "Density plot of residuals"))
    
    #Homoskedasticity (residuals vs fitted values)
    #print(plot(fitted(phy_gls_lam), phy_gls_lam$residuals, xlab = "Fitted", ylab = "Residuals"))
    #print(abline(0,0, col = "red"))
    #dev.off()
    
    ###and do for best model too
    #jpeg(paste0("./output/PCAs/iterations/model_assessment_pc", j, "_iteration", i, "_summary_phylolm", best_model, "_best_model.jpg"), width = 8, height = 6, res = 500, units = "in")
    #par(mfrow=c(2,2))
    
    #plot residuals vs independent variable
    #print(plot(data$X, get(paste0("phy_gls_", best_model))$residuals, ylab = "Residuals", xlab = "Variable"))
    
    #Distribution of residuals (normal)
    #qplot <- qqnorm(get(paste0("phy_gls_", best_model))$residuals)
    #print(qplot)
    #print(qqline(get(paste0("phy_gls_", best_model))$residuals, col = "red"))
    
    #density plot
    #print(plot(density(get(paste0("phy_gls_", best_model))$residuals, main = "Density plot of residuals")))
    
    #Homoskedasticity (residuals vs fitted values)
    #print(plot(fitted(get(paste0("phy_gls_", best_model))), get(paste0("phy_gls_", best_model))$residuals, xlab = "Fitted", ylab = "Residuals"))
    #print(abline(0,0, col = "red"))
    #dev.off()
    
    
  }
  
  print(i)
  if (i == 1){
    output_gls <- output_gls_per
    output_best <- output_gls_best
    
  } else {
    output_gls <- rbind(output_gls, output_gls_per)
    output_best <- rbind(output_best, output_gls_best)
  }
}

#saveRDS(output_gls, "./analysis/pca_analysis/pc_univariate_model_gls_parameters_ang_only.rds")
#saveRDS(output_best, "./analysis/pca_analysis/pc_univariate_model_gls_parameters_best_models_ang_only.rds")

output_gls <- readRDS("./analysis/pca_analysis/pc_univariate_model_gls_parameters_ang_only.rds")
output_best <- readRDS("./analysis/pca_analysis/pc_univariate_model_gls_parameters_best_models_ang_only.rds")


#plot distribution of pvalues for phylo models
jpeg("./output/PCAs/phylo_gls_pvalues_pcas_ang_only.jpg", width = 12, height = 8, units = "in", res = 600)
ggplot(output_gls, aes(x = as.factor(PC_axis), y = round(p_value, 3)))+
  #geom_point()+
  #geom_boxplot(aes(fill = as.factor(name)))+
  #geom_jitter(color = "black")+
  geom_jitter(aes(color = ifelse(p_value < 0.05, "red", "black"))) +
  scale_color_identity()+
  labs(y = "p-value", x = "PC Axis", title = "Phylogenetic GLS")+
  geom_hline(aes(yintercept = 0.05), colour = "red")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  theme_bw()
dev.off()

#plot distribution of pvalues for phylo models
jpeg("./output/PCAs/phylolm_int_model_pvalues_pcas_ang_only.jpg", width = 12, height = 8, units = "in", res = 600)
ggplot(output_gls, aes(x = as.factor(PC_axis), y = round(p_value, 3)))+
  #geom_point()+
  #geom_boxplot(aes(fill = as.factor(name)))+
  #geom_jitter(color = "black")+
  geom_jitter(aes(color = ifelse(p_value < 0.05, "red", "black"))) +
  scale_color_identity()+
  labs(y = "p-value", x = "PC Axis", title = "Phylogenetic LM - original model")+
  geom_hline(aes(yintercept = 0.05), colour = "red")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  theme_bw()
dev.off()

#plot distribution of pvalues for phylo models
#jpeg("./output/PCAs/phylolm_best_model_pvalues_pcas_ang_only.jpg", width = 6, height = 6, units = "in", res = 600)
pdf("./output/PCAs/phylolm_best_model_pvalues_pcas_ang_only.pdf", width = 6, height = 6)
ggplot(output_best[which(!is.na(output_best$p_value) & output_best$PC_axis %in% c(1:5)),], aes(x = as.factor(PC_axis), y = round(p_value, 3)))+
  #geom_point()+
  #geom_boxplot(aes(fill = as.factor(name)))+
  #geom_jitter(color = "black")+
  geom_jitter(aes(color = ifelse(p_value < 0.05, "red", "black")), size = 0.75) +
  scale_color_identity()+
  labs(y = "p-value", x = "PC Axis")+#, title = "Phylogenetic Linear Regressions"
  geom_hline(aes(yintercept = 0.05), colour = "red")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  theme_bw()
dev.off()
#thinking about PCAs
#keep axes with eigenvalue (variance) greater than 1

#less than 0.1 correlation for top axes
cor(pca_best_model$scores[,1], pca_best_model$scores[,2]) #0.072
cor(pca_best_model$scores[,1], pca_best_model$scores[,3]) #-0.076
cor(pca_best_model$scores[,1], pca_best_model$scores[,4]) #0.036
cor(pca_best_model$scores[,2], pca_best_model$scores[,3]) #0.0125
cor(pca_best_model$scores[,3], pca_best_model$scores[,4]) #0.010

#plot eigenvalues for pc axis

str(pca_best_model$values)
str(pca_best_model)
plot(colnames(data_spectra$spectra), pca_best_model$values)

pca_best_model$scores[,1]

plot(simmap_tree[[1]])

plotTree.wBars(simmap_tree[[1]],pca_best_model$values)

obj<-simmap_tree[[1]]

class(simmap_tree[[1]])


#Map trait to trees - make simmap tree for predictor
myc_named <- setNames(data_spectra$myc, rownames(data_spectra$spectra))

consensus_tree <- readRDS("./data/for_analysis/myc_consensus_tree.rds")

tips_to_drop <- consensus_tree$tip.label[-which(consensus_tree$tip.label %in% data_spectra$species)]

consensus_tree_pruned <- drop.tip(consensus_tree, tip = tips_to_drop)

simmap_consensus <- make.simmap(consensus_tree_pruned, myc_named, model="SYM", nsim=100)

summary_simmap <-describe.simmap(simmap_consensus,plot=TRUE,cex=0.7)

plot(summary_simmap$tree[1])
class(summary_simmap$tree)

# 
# #assign states to nodes
# for (i in 1:length(simmap_tree)){
#   simmap_tree[[i]]$node.label<-getStates(simmap_tree[[i]],"nodes")  
# }

class(summary_simmap)

colour_bars <- setNames(data_spectra$myc, data_spectra$species)
colour_bars <- gsub("AM", "orange", colour_bars)
colour_bars <- gsub("EM", "blue", colour_bars)

cols_tree<-setNames(c("blue","orange"),unique(myc_named))
#plot summary of 100 reps of stochastic mapping on consensus tree
plot(summary(simmap_consensus),colors=cols_tree,ftype="i")

#plot single map on consensus tree with scores for pc axis 1 and again for pc2
jpeg("./output/PCAs/consensus_plots/consensus_myc_pca1_iteration100_OUmodel.jpg", height = 10, width = 6, units = "in", res = 600)
plotTree.wBars(simmap_consensus[[1]],pca_best_model$scores[,1],method="plotSimmap",
               tip.labels=TRUE,fsize=0.7, col = colour_bars, colors = cols_tree)
dev.off()

jpeg("./output/PCAs/consensus_plots/consensus_myc_pca2_iteration100_OUmodel.jpg", height = 10, width = 6, units = "in", res = 600)

plotTree.wBars(simmap_consensus[[1]],pca_best_model$scores[,2],method="plotSimmap",
               tip.labels=TRUE,fsize=0.7, col = colour_bars, colors = cols_tree)
dev.off()


#plot just bars for subsequent axes and only one tree
sorted_scores <- pca_best_model$scores[,1][order(match(names(pca_best_model$scores[,1]),simmap_consensus[[1]]$tip.label))]

plotTree.barplot(simmap_consensus[[1]],sorted_scores, add=TRUE, args.barplot=list(xlab="PC 1",col=colour_bars, colors = cols_tree, mar=c(5.1,0,2.1,2.1)))
plotTree.barplot(simmap_consensus[[1]],pca_best_model$scores[,2],args.barplot=list(xlab="PC 2",mar=c(5.1,0,2.1,2.1),
                                          col=cols),args.plotTree=list(plot=FALSE),add=TRUE)
plotTree.barplot(simmap_consensus[[1]],pca_best_model$scores[,3],args.barplot=list(xlab="PC 3",mar=c(5.1,0,2.1,2.1),
                                          col=cols),args.plotTree=list(plot=FALSE),add=TRUE)





