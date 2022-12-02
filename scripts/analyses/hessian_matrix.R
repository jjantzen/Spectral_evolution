#intercept model PCAs
library(mvMORPH)
library(phytools)
library(dplyr)
library(ggplot2)
library(OUwie)
library(tidyr)
library(phylolm)
library(phytools)
source("./scripts/plot_phylo_pca_edited_function.R")
#load plot_pca_phylo function

#read in data
data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")
tree_myc <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

model_list <- readRDS("./analysis/pca_analysis/best_intercept_models_for_pcas_92sp.rds")

myc_named <- setNames(data_spectra$myc, rownames(data_spectra$spectra))

simmap_tree <- readRDS("./analysis/pca_analysis/myc_simmap_trees.rds")

#get best models
aics_ouwie <- readRDS("./analysis/pca_analysis/pc_univariate_model_aics_92sp.rds")

#get column names with highest and second highest value for columns 3:9
best_ouwie_models <- aics_ouwie %>%
  dplyr::rowwise() %>%
  dplyr::mutate(top_model = names(.)[which.max(c_across(BM1:OUMVA))+2], second_model = tail(head(names(cur_data())[order(c_across(BM1:OUMVA), decreasing = T)+2],2),1))#(.[3:9], 1, function(x) names(x)[maxn(2)(x)]))

best_ouwie_models

#run models and get output
for (i in 1:length(model_list)){ #test with 1 first
  #conduct PCA
  pca_best_model <- mvgls.pca(model_list[[i]], plot = FALSE)  
  
  output_aics_per <- data.frame(matrix(nrow=3, ncol = 27))
  colnames(output_aics_per)[1] <- "iteration"
  class(output_aics_per[,1]) <- "numeric"
  colnames(output_aics_per)[2] <- "PC_axis"
  class(output_aics_per[,2]) <- "numeric"
  colnames(output_aics_per)[3] <- "model"
  class(output_aics_per[,3]) <- "character"
  colnames(output_aics_per)[4] <- "AICc"
  class(output_aics_per[,4]) <- "numeric"
  colnames(output_aics_per)[5] <- "lnL"
  class(output_aics_per[,5]) <- "numeric"
  colnames(output_aics_per)[6] <- "AIC"
  class(output_aics_per[,6]) <- "numeric"
  colnames(output_aics_per)[7] <- "eigenval_1"
  class(output_aics_per[,7]) <- "numeric"
  colnames(output_aics_per)[8] <- "eigenval_2"
  class(output_aics_per[,8]) <- "numeric"
  colnames(output_aics_per)[9] <- "BIC"
  class(output_aics_per[,9]) <- "numeric"
  colnames(output_aics_per)[10] <- "SE_theta_single"
  class(output_aics_per[,10]) <- "numeric"
  colnames(output_aics_per)[11] <- "SE_alpha_single"
  class(output_aics_per[,11]) <- "numeric"
  colnames(output_aics_per)[12] <- "SE_sigsq_single"
  class(output_aics_per[,12]) <- "numeric"
  colnames(output_aics_per)[13] <- "alpha_single"
  class(output_aics_per[,13]) <- "numeric"
  colnames(output_aics_per)[14] <- "sigma.sq_single"
  class(output_aics_per[,14]) <- "numeric"
  colnames(output_aics_per)[15] <- "theta_single"
  class(output_aics_per[,15]) <- "numeric"
  colnames(output_aics_per)[16] <- "SE_theta_AM"
  class(output_aics_per[,16]) <- "numeric"
  colnames(output_aics_per)[17] <- "SE_alpha_AM"
  class(output_aics_per[,17]) <- "numeric"
  colnames(output_aics_per)[18] <- "SE_sigsq_AM"
  class(output_aics_per[,18]) <- "numeric"
  colnames(output_aics_per)[19] <- "SE_theta_EM"
  class(output_aics_per[,19]) <- "numeric"
  colnames(output_aics_per)[20] <- "SE_alpha_EM"
  class(output_aics_per[,20]) <- "numeric"
  colnames(output_aics_per)[21] <- "SE_sigsq_EM"
  class(output_aics_per[,21]) <- "numeric"
  colnames(output_aics_per)[22] <- "alpha_AM"
  class(output_aics_per[,22]) <- "numeric"
  colnames(output_aics_per)[23] <- "sigma.sq_AM"
  class(output_aics_per[,23]) <- "numeric"
  colnames(output_aics_per)[24] <- "theta_AM"
  class(output_aics_per[,24]) <- "numeric"
  colnames(output_aics_per)[25] <- "alpha_EM"
  class(output_aics_per[,25]) <- "numeric"
  colnames(output_aics_per)[26] <- "sigma.sq_EM"
  class(output_aics_per[,26]) <- "numeric"
  colnames(output_aics_per)[27] <- "theta_EM"
  class(output_aics_per[,27]) <- "numeric"
  
  #run ouwie models
  for (j in 1:3){
    #get data
    data <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named, X=as.numeric(pca_best_model$scores[,j]))
    
    #run models for each pc
    best_model <- best_ouwie_models$top_model[which(best_ouwie_models$iteration == i & best_ouwie_models$PC_axis == j)]
    best_model_output <- OUwie(simmap_tree[[i]], data, model = best_model, simmap.tree = TRUE, diagn = TRUE) 
    print(j)

    #save output of model
    output_aics_per$iteration[j] <- i
    output_aics_per$PC_axis[j] <- j
    output_aics_per$model[j] <- best_model
    #$solution.se #the standard error will be small, and the parameter estimate is considered stable
    
    if (best_model == "OU1"){
    output_aics_per$SE_theta_single[j] <- paste0(best_model_output$theta[2])
    output_aics_per$SE_alpha_single[j] <- paste0(best_model_output$solution.se[1])
    output_aics_per$SE_sigsq_single[j] <- paste0(best_model_output$solution.se[2])
    output_aics_per$eigenval_1[j] <- best_model_output$eigval[1]
    output_aics_per$eigenval_2[j] <- best_model_output$eigval[2]
    output_aics_per$AIC[j] <- best_model_output$AIC
    output_aics_per$AICc[j] <- best_model_output$AICc
    output_aics_per$lnL[j] <- best_model_output$loglik
    output_aics_per$BIC[j] <- best_model_output$BIC
    output_aics_per$alpha_single[j] <- paste0(best_model_output$solution[1]) 
    output_aics_per$sigma.sq_single[j] <- paste0(best_model_output$solution[2]) 
    output_aics_per$theta_single[j] <- paste0(best_model_output$theta[1]) 
    } else {
    output_aics_per$SE_theta_AM[j] <- best_model_output$theta[1,2]
    output_aics_per$SE_alpha_AM[j] <- best_model_output$solution.se[1,1]
    output_aics_per$SE_sigsq_AM[j] <- best_model_output$solution.se[2,1]
    output_aics_per$SE_theta_EM[j] <- best_model_output$theta[2,2]
    output_aics_per$SE_alpha_EM[j] <- best_model_output$solution.se[1,2]
    output_aics_per$SE_sigsq_EM[j] <- best_model_output$solution.se[2,2]
    output_aics_per$eigenval_1[j] <- best_model_output$eigval[1]
    output_aics_per$eigenval_2[j] <- best_model_output$eigval[2]
    output_aics_per$AIC[j] <- best_model_output$AIC
    output_aics_per$AICc[j] <- best_model_output$AICc
    output_aics_per$lnL[j] <- best_model_output$loglik
    output_aics_per$BIC[j] <- best_model_output$BIC
    output_aics_per$alpha_AM[j] <- best_model_output$solution[1,1]
    output_aics_per$sigma.sq_AM[j] <- best_model_output$solution[2,1]
    output_aics_per$theta_AM[j] <- best_model_output$theta[1,1]
    output_aics_per$alpha_EM[j] <- best_model_output$solution[1,2]
    output_aics_per$sigma.sq_EM[j] <- best_model_output$solution[2,2]
    output_aics_per$theta_EM[j] <- best_model_output$theta[2,1]
      
    }
    
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


saveRDS(output_aics, "./analysis/pca_analysis/pc_univariate_model_92sp_hessian_output.rds")

#$solution.se #the standard error will be small, and the parameter estimate is considered stable
#$eigval # If all the eigenvalues of the Hessian are positive, then the Hessian is positive definite, and all parameter estimates are considered reliable

output_aics <- readRDS("./analysis/pca_analysis/pc_univariate_model_92sp_hessian_output.rds")


#check the max of SE values

min(unique(!is.na(output_aics$SE_alpha_AM)))

colnames(output_aics)

#alpha SE
output_aics[which(!is.na(output_aics$SE_alpha_AM) & output_aics$SE_alpha_AM > 0.8),c(1,2,3, 17)] #big values of SE alpha AM #19-2, 43-3, 58-4

output_aics[which(!is.na(output_aics$SE_alpha_EM) & output_aics$SE_alpha_EM > 0.8),c(1,2,3, 20)] #big values of SE alpha EM #19-2, 43-3, 58-4, 69-2

output_aics[which(!is.na(output_aics$SE_alpha_single) & output_aics$SE_alpha_single > 0.8),c(1,2,3,11)] #big values of SE alpha single #82-6, 83-6, 92-2, 93-5 #also 66-6 and 96-6 but I think not

#theta SE (10, 16, 19)
nrow(output_aics[which(!is.na(output_aics$SE_theta_AM) & output_aics$SE_theta_AM > 0.8),c(1,2,3, 16)]) #big values of SE alpha AM #22 different rows
nrow(output_aics[which(!is.na(output_aics$SE_theta_EM) & output_aics$SE_theta_EM > 0.8),c(1,2,3, 19)]) #big values of SE alpha AM #105 different rows
output_aics[which(!is.na(output_aics$SE_theta_single) & output_aics$SE_theta_single > 0.8),c(1,2,3, 10)] #big values of SE alpha AM #none

#sigma.sq SE (12, 18, 21)
output_aics[which(!is.na(output_aics$SE_sigsq_AM) & output_aics$SE_sigsq_AM > 0.8),c(1,2,3, 18)] #big values of SE alpha AM #19-2, 43-3, 49-5, 58-4, 58-5
nrow(output_aics[which(!is.na(output_aics$SE_sigsq_EM) & output_aics$SE_sigsq_EM > 0.8),c(1,2,3, 21)]) #big values of SE alpha AM #40 different rows
output_aics[which(!is.na(output_aics$SE_sigsq_single) & output_aics$SE_sigsq_single > 0.8),c(1,2,3, 12)] #big values of SE alpha AM #82-6, 83-6, 92-2, 93-6

#add index for comparison later
output_aics_top3 <- output_aics[which(output_aics$PC_axis %in% c(1:3)),]

id <- rownames(output_aics_top3)
output_aics_top3 <- cbind(id=id, output_aics_top3)

#output_aics <- mutate(output_aics, ID = row_number())
colnames(output_aics_top3)

#check if any negative values in eigenval_1 or eigenval_2
eigenvals_neg <- output_aics_top3[which(output_aics_top3$eigenval_1 < 0 | output_aics_top3$eigenval_2 < 0),c(1,2,3,4, 8,9)] #only three for which a negative eigenval #27-4 (neg 2), 44-5 (neg 2) and 95-3 (neg 2)

#plot which theta points are unstable on origianl plot?
AM_unstable <- output_aics_top3[which(!is.na(output_aics_top3$SE_theta_AM) & output_aics_top3$SE_theta_AM > 1 & output_aics_top3$PC_axis %in% c(1,2,3)),c(1,2,3,4, 17)] #big values of SE alpha AM #22 different rows
EM_unstable <- output_aics_top3[which(!is.na(output_aics_top3$SE_theta_EM) & output_aics_top3$SE_theta_EM > 1 & output_aics_top3$PC_axis %in% c(1,2,3)),c(1,2,3,4, 20)] #big values of SE alpha AM #105 different rows

#should I rerun the unstable models with simpler model? eg sub OUMV for OUMVA
output_aics_top3[which(output_aics_top3$id %in% eigenvals_neg$id),]

output_aics_top3[which(output_aics_top3$id %in% AM_unstable$id),]

str(output_aics)

#colour by se
colnames(output_aics)

long_output_parameters <- output_aics %>% 
  dplyr::select(iteration, PC_axis, model, theta_AM, theta_EM, SE_theta_AM, SE_theta_EM, eigenval_2) %>% 
  pivot_longer(cols = c(theta_AM, theta_EM)) %>% 
  dplyr::filter(between(value, -10, 10))



jpeg("./analysis/pca_analysis/fig_SE_vs_thetas_OUM_reps.jpg", width = 8, height = 6, res = 600, units = "in")
ggplot() +
  geom_point(data = long_output_parameters[which(long_output_parameters$SE_theta_EM <1),], aes(x = iteration, y = value), colour = "blue")+
  geom_point(data = long_output_parameters[which(long_output_parameters$SE_theta_EM >=1),], aes(x = iteration, y = value), colour = "red")+
  facet_grid(name~PC_axis, scales = "free")
dev.off()

#exclude bad SE ones and plot thetas and optima








########plots########
output_parameters_original <- readRDS("./analysis/pca_analysis/summary_of_thetas_alphas_top3pcas.rds")

id_plot <- rownames(output_parameters_original)
output_parameters_all <- cbind(id=id_plot, output_parameters_original)



colnames(output_parameters_original) <- paste0(colnames(output_parameters_original), ".plot")

#combine dataframes

combo_output <- cbind(output_aics_top3, output_parameters_original)
colnames(combo_output)


long_output_parameters <- combo_output %>% 
  dplyr::select(id, iteration, PC_axis, model, theta_single, theta_AM, theta_EM, AM_theta.plot, EM_theta.plot) %>% 
  pivot_longer(cols = c(theta_AM, theta_EM, AM_theta.plot, EM_theta.plot))

#get rid of outliers
no_outliers <- long_output_parameters %>% 
  #dplyr::filter(!model %in% "OU1") %>% 
  dplyr::filter(between(value, -10, 10))

ggplot(no_outliers)+
  geom_boxplot(aes(x = name, y = value, fill = name))+
  facet_grid(PC_axis ~ model)

ggplot(combo_output)+
  geom_point(aes(x = id, y = theta_EM), colour = "blue")+
  geom_point(aes(x=id, y = EM_theta.plot), colour = "red")+
  facet_wrap(~PC_axis, nrow = 1, scales = "free")

colnames(no_outliers)

#difference between theta for two runs

differences_long_no_outliers <- combo_output %>% 
  dplyr::mutate(diff_AM = theta_AM-AM_theta.plot, diff_EM = theta_EM-EM_theta.plot) %>% 
  dplyr::select(id, iteration, PC_axis, model, diff_AM, diff_EM) %>% 
  pivot_longer(cols = c(diff_AM, diff_EM)) %>% 
  dplyr::filter(between(value, -10, 10))

differences_long_inc_outliers <- combo_output %>% 
  dplyr::mutate(diff_AM = theta_AM-AM_theta.plot, diff_EM = theta_EM-EM_theta.plot) %>% 
  dplyr::select(id, iteration, PC_axis, model, diff_AM, diff_EM) %>% 
  pivot_longer(cols = c(diff_AM, diff_EM))



ggplot(differences_long_no_outliers)+
  geom_point(aes(x = name, y = value, colour = name))+
  facet_wrap(~PC_axis, nrow = 1, scales = "free")


#get id for those where big difference

pos <- differences_long_inc_outliers$id[which(differences_long_inc_outliers$value > 0.5)]
negs <- differences_long_inc_outliers$id[which(differences_long_inc_outliers$value < -0.5)]

ids <- c(AM_unstable$id, EM_unstable$id)

ids_dif <- c(as.integer(pos), as.integer(negs))

ids_dif %in% ids #not just unstable parameters are highly different from other estimate

ids %in% ids_dif #most of unstable parameters are highly different from other estimates but vice versa isn't true

colnames(combo_output)

#rerun pca and models for ten trees and three pc axes to get distribution (but only single pca done per tree)
i <- 1
#run models and get output
for (i in 1:10){#length(model_list)){ #test with 1 first
  #conduct PCA
  pca_best_model <- mvgls.pca(model_list[[i]], plot = FALSE)  
  
  
  #run ouwie models
  for (j in 1:3){
    #repeat 10 times?
    
    output_param <- data.frame(matrix(nrow=10, ncol = 28))
    colnames(output_param)[1] <- "tree"
    class(output_param[,1]) <- "numeric"
    colnames(output_param)[2] <- "PC_axis"
    class(output_param[,2]) <- "numeric"
    colnames(output_param)[3] <- "model"
    class(output_param[,3]) <- "character"
    colnames(output_param)[4] <- "AICc"
    class(output_param[,4]) <- "numeric"
    colnames(output_param)[5] <- "lnL"
    class(output_param[,5]) <- "numeric"
    colnames(output_param)[6] <- "AIC"
    class(output_param[,6]) <- "numeric"
    colnames(output_param)[7] <- "eigenval_1"
    class(output_param[,7]) <- "numeric"
    colnames(output_param)[8] <- "eigenval_2"
    class(output_param[,8]) <- "numeric"
    colnames(output_param)[9] <- "BIC"
    class(output_param[,9]) <- "numeric"
    colnames(output_param)[10] <- "SE_theta_single"
    class(output_param[,10]) <- "numeric"
    colnames(output_param)[11] <- "SE_alpha_single"
    class(output_param[,11]) <- "numeric"
    colnames(output_param)[12] <- "SE_sigsq_single"
    class(output_param[,12]) <- "numeric"
    colnames(output_param)[13] <- "alpha_single"
    class(output_param[,13]) <- "numeric"
    colnames(output_param)[14] <- "sigma.sq_single"
    class(output_param[,14]) <- "numeric"
    colnames(output_param)[15] <- "theta_single"
    class(output_param[,15]) <- "numeric"
    colnames(output_param)[16] <- "SE_theta_AM"
    class(output_param[,16]) <- "numeric"
    colnames(output_param)[17] <- "SE_alpha_AM"
    class(output_param[,17]) <- "numeric"
    colnames(output_param)[18] <- "SE_sigsq_AM"
    class(output_param[,18]) <- "numeric"
    colnames(output_param)[19] <- "SE_theta_EM"
    class(output_param[,19]) <- "numeric"
    colnames(output_param)[20] <- "SE_alpha_EM"
    class(output_param[,20]) <- "numeric"
    colnames(output_param)[21] <- "SE_sigsq_EM"
    class(output_param[,21]) <- "numeric"
    colnames(output_param)[22] <- "alpha_AM"
    class(output_param[,22]) <- "numeric"
    colnames(output_param)[23] <- "sigma.sq_AM"
    class(output_param[,23]) <- "numeric"
    colnames(output_param)[24] <- "theta_AM"
    class(output_param[,24]) <- "numeric"
    colnames(output_param)[25] <- "alpha_EM"
    class(output_param[,25]) <- "numeric"
    colnames(output_param)[26] <- "sigma.sq_EM"
    class(output_param[,26]) <- "numeric"
    colnames(output_param)[27] <- "theta_EM"
    class(output_param[,27]) <- "numeric"
    colnames(output_param)[28] <- "rep"
    class(output_param)[28] <- "numeric"
    
    output_aics_per <- data.frame(matrix(nrow=10, ncol = 11))
    colnames(output_aics_per)[1] <- "tree"
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
    colnames(output_aics_per)[10] <- "top_model"
    class(output_aics_per)[10] <- "numeric"
    colnames(output_aics_per)[11] <- "rep"
    class(output_aics_per)[11] <- "numeric"
    
    for (k in 1:10){
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
      print(paste("PC", j))
      #save aic.w scores
      aics <- setNames(c(BM1_pc$AIC, BMS_pc$AIC,OU1_pc$AIC, OUM_pc$AIC, OUMV_pc$AIC, OUMA_pc$AIC, OUMVA_pc$AIC), c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"))
      weights <- aic.w(aics)
      
      output_aics_per$tree[k] <- i
      output_aics_per$PC_axis[k] <- j
      output_aics_per$BM1[k] <- weights[1]
      output_aics_per$BMS[k] <- weights[2]
      output_aics_per$OU1[k] <- weights[3]
      output_aics_per$OUM[k] <- weights[4]
      output_aics_per$OUMV[k] <- weights[5]
      output_aics_per$OUMA[k] <- weights[6]
      output_aics_per$OUMVA[k] <- weights[7]
      output_aics_per$rep[k] <- k
      
      
      
      # best_ouwie_models <- output_aics_per[k,] %>%
      #   dplyr::mutate(top_model = names(.)[which.max(c_across(BM1:OUMVA))+2], second_model = tail(head(names(cur_data())[order(c_across(BM1:OUMVA), decreasing = T)+2],2),1))#(.[3:9], 1, function(x) names(x)[maxn(2)(x)]))
      
      output_aics_per$top_model[k] <- colnames(output_aics_per)[apply(output_aics_per[k,3:9],1,which.max)+2]
      
      #based on all models, choose best model
      best_model <- output_aics_per$top_model[which(output_aics_per$tree == i & output_aics_per$PC_axis == j)]
      best_model_output <- OUwie(simmap_tree[[i]], data, model = best_model, simmap.tree = TRUE, diagn = TRUE) 
      print(j)
      
      #save output of model
      output_param$tree[k] <- i
      output_param$PC_axis[k] <- j
      output_param$model[k] <- best_model
      #$solution.se #the standard error will be small, and the parameter estimate is considered stable
      
      if (best_model == "OU1"){
        output_param$SE_theta_single[k] <- paste0(best_model_output$theta[2])
        output_param$SE_alpha_single[k] <- paste0(best_model_output$solution.se[1])
        output_param$SE_sigsq_single[k] <- paste0(best_model_output$solution.se[2])
        output_param$eigenval_1[k] <- best_model_output$eigval[1]
        output_param$eigenval_2[k] <- best_model_output$eigval[2]
        output_param$AIC[k] <- best_model_output$AIC
        output_param$AICc[k] <- best_model_output$AICc
        output_param$lnL[k] <- best_model_output$loglik
        output_param$BIC[k] <- best_model_output$BIC
        output_param$alpha_single[k] <- paste0(best_model_output$solution[1]) 
        output_param$sigma.sq_single[k] <- paste0(best_model_output$solution[2]) 
        output_param$theta_single[k] <- paste0(best_model_output$theta[1]) 
        output_param$rep[k] <- k
      } else {
        output_param$SE_theta_AM[k] <- best_model_output$theta[1,2]
        output_param$SE_alpha_AM[k] <- best_model_output$solution.se[1,1]
        output_param$SE_sigsq_AM[k] <- best_model_output$solution.se[2,1]
        output_param$SE_theta_EM[k] <- best_model_output$theta[2,2]
        output_param$SE_alpha_EM[k] <- best_model_output$solution.se[1,2]
        output_param$SE_sigsq_EM[k] <- best_model_output$solution.se[2,2]
        output_param$eigenval_1[k] <- best_model_output$eigval[1]
        output_param$eigenval_2[k] <- best_model_output$eigval[2]
        output_param$AIC[k] <- best_model_output$AIC
        output_param$AICc[k] <- best_model_output$AICc
        output_param$lnL[k] <- best_model_output$loglik
        output_param$BIC[k] <- best_model_output$BIC
        output_param$alpha_AM[k] <- best_model_output$solution[1,1]
        output_param$sigma.sq_AM[k] <- best_model_output$solution[2,1]
        output_param$theta_AM[k] <- best_model_output$solution[2,1]
        output_param$alpha_EM[k] <- best_model_output$solution[1,2]
        output_param$sigma.sq_EM[k] <- best_model_output$solution[2,2]
        output_param$theta_EM[k] <- best_model_output$solution[2,2]
        output_param$rep[k] <- k
        
      }
      
    }
    if (j == 1){
      output_aics <- output_aics_per
      output_params_pcs <- output_param
    } else {
      output_aics <-  rbind(output_aics, output_aics_per)
      output_params_pcs <- rbind(output_params_pcs, output_param)
    }
    
  }
  print(paste("tree", i))
  if (i == 1){
    output_aics_all <- output_aics
    output_params_all <- output_params_pcs
    
  } else {
    output_aics_all <- rbind(output_aics_all, output_aics)
    output_params_all <- rbind(output_params_all, output_params_pcs)
  }
}


saveRDS(output_params_all, "./analysis/pca_analysis/testing_replicability_ouwie_parameters.rds")
saveRDS(output_aics_all, "./analysis/pca_analysis/testing_replicability_ouwie_aics.rds")


output_params_all <- readRDS("./analysis/pca_analysis/testing_replicability_ouwie_parameters.rds")

output_params_all[which(output_params_all$tree == 3),]

#and with redoing pc for each

for (i in 1:10){#length(model_list)){ #test with 1 first
  #conduct PCA
  
  #run ouwie models
  for (j in 1:3){
    #repeat 10 times?
    
    output_param <- data.frame(matrix(nrow=10, ncol = 28))
    colnames(output_param)[1] <- "tree"
    class(output_param[,1]) <- "numeric"
    colnames(output_param)[2] <- "PC_axis"
    class(output_param[,2]) <- "numeric"
    colnames(output_param)[3] <- "model"
    class(output_param[,3]) <- "character"
    colnames(output_param)[4] <- "AICc"
    class(output_param[,4]) <- "numeric"
    colnames(output_param)[5] <- "lnL"
    class(output_param[,5]) <- "numeric"
    colnames(output_param)[6] <- "AIC"
    class(output_param[,6]) <- "numeric"
    colnames(output_param)[7] <- "eigenval_1"
    class(output_param[,7]) <- "numeric"
    colnames(output_param)[8] <- "eigenval_2"
    class(output_param[,8]) <- "numeric"
    colnames(output_param)[9] <- "BIC"
    class(output_param[,9]) <- "numeric"
    colnames(output_param)[10] <- "SE_theta_single"
    class(output_param[,10]) <- "numeric"
    colnames(output_param)[11] <- "SE_alpha_single"
    class(output_param[,11]) <- "numeric"
    colnames(output_param)[12] <- "SE_sigsq_single"
    class(output_param[,12]) <- "numeric"
    colnames(output_param)[13] <- "alpha_single"
    class(output_param[,13]) <- "numeric"
    colnames(output_param)[14] <- "sigma.sq_single"
    class(output_param[,14]) <- "numeric"
    colnames(output_param)[15] <- "theta_single"
    class(output_param[,15]) <- "numeric"
    colnames(output_param)[16] <- "SE_theta_AM"
    class(output_param[,16]) <- "numeric"
    colnames(output_param)[17] <- "SE_alpha_AM"
    class(output_param[,17]) <- "numeric"
    colnames(output_param)[18] <- "SE_sigsq_AM"
    class(output_param[,18]) <- "numeric"
    colnames(output_param)[19] <- "SE_theta_EM"
    class(output_param[,19]) <- "numeric"
    colnames(output_param)[20] <- "SE_alpha_EM"
    class(output_param[,20]) <- "numeric"
    colnames(output_param)[21] <- "SE_sigsq_EM"
    class(output_param[,21]) <- "numeric"
    colnames(output_param)[22] <- "alpha_AM"
    class(output_param[,22]) <- "numeric"
    colnames(output_param)[23] <- "sigma.sq_AM"
    class(output_param[,23]) <- "numeric"
    colnames(output_param)[24] <- "theta_AM"
    class(output_param[,24]) <- "numeric"
    colnames(output_param)[25] <- "alpha_EM"
    class(output_param[,25]) <- "numeric"
    colnames(output_param)[26] <- "sigma.sq_EM"
    class(output_param[,26]) <- "numeric"
    colnames(output_param)[27] <- "theta_EM"
    class(output_param[,27]) <- "numeric"
    colnames(output_param)[28] <- "rep"
    class(output_param)[28] <- "numeric"
    
    output_aics_per <- data.frame(matrix(nrow=10, ncol = 11))
    colnames(output_aics_per)[1] <- "tree"
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
    colnames(output_aics_per)[10] <- "top_model"
    class(output_aics_per)[10] <- "numeric"
    colnames(output_aics_per)[11] <- "rep"
    class(output_aics_per)[11] <- "numeric"
    
    for (k in 1:10){
      #get data
      pca_best_model <- mvgls.pca(model_list[[i]], plot = FALSE)  
      
      data <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named, X=as.numeric(pca_best_model$scores[,j]))
      
      #run models for each pc
      
      BM1_pc <- OUwie(simmap_tree[[i]], data, model = "BM1", simmap.tree = TRUE) #single rate
      BMS_pc <- OUwie(simmap_tree[[i]], data, model = "BMS", simmap.tree = TRUE) #multiple rates
      OU1_pc <- OUwie(simmap_tree[[i]], data, model = "OU1", simmap.tree = TRUE) #single optimum
      OUM_pc <- OUwie(simmap_tree[[i]], data, model = "OUM", simmap.tree = TRUE) #multiple optima
      OUMV_pc <- OUwie(simmap_tree[[i]], data, model = "OUMV", simmap.tree = TRUE) #multiple optimum and multiple sigmas
      OUMA_pc <- OUwie(simmap_tree[[i]], data, model = "OUMA", simmap.tree = TRUE) #multiple optimum and multiple alphas
      OUMVA_pc <- OUwie(simmap_tree[[i]], data, model = "OUMVA", simmap.tree = TRUE) #multiple optimum and multiple sigmas and multiple alphas
      print(paste("PC", j))
      #save aic.w scores
      aics <- setNames(c(BM1_pc$AIC, BMS_pc$AIC,OU1_pc$AIC, OUM_pc$AIC, OUMV_pc$AIC, OUMA_pc$AIC, OUMVA_pc$AIC), c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"))
      weights <- aic.w(aics)
      
      output_aics_per$tree[k] <- i
      output_aics_per$PC_axis[k] <- j
      output_aics_per$BM1[k] <- weights[1]
      output_aics_per$BMS[k] <- weights[2]
      output_aics_per$OU1[k] <- weights[3]
      output_aics_per$OUM[k] <- weights[4]
      output_aics_per$OUMV[k] <- weights[5]
      output_aics_per$OUMA[k] <- weights[6]
      output_aics_per$OUMVA[k] <- weights[7]
      output_aics_per$rep[k] <- k
      
      
      
      # best_ouwie_models <- output_aics_per[k,] %>%
      #   dplyr::mutate(top_model = names(.)[which.max(c_across(BM1:OUMVA))+2], second_model = tail(head(names(cur_data())[order(c_across(BM1:OUMVA), decreasing = T)+2],2),1))#(.[3:9], 1, function(x) names(x)[maxn(2)(x)]))
      
      output_aics_per$top_model[k] <- colnames(output_aics_per)[apply(output_aics_per[k,3:9],1,which.max)+2]
      
      #based on all models, choose best model
      best_model <- output_aics_per$top_model[which(output_aics_per$tree == i & output_aics_per$PC_axis == j)]
      best_model_output <- OUwie(simmap_tree[[i]], data, model = best_model, simmap.tree = TRUE, diagn = TRUE) 
      print(paste0("rep", k))
      
      #save output of model
      output_param$tree[k] <- i
      output_param$PC_axis[k] <- j
      output_param$model[k] <- best_model
      #$solution.se #the standard error will be small, and the parameter estimate is considered stable
      
      if (best_model == "OU1"){
        output_param$SE_theta_single[k] <- paste0(best_model_output$theta[2])
        output_param$SE_alpha_single[k] <- paste0(best_model_output$solution.se[1])
        output_param$SE_sigsq_single[k] <- paste0(best_model_output$solution.se[2])
        output_param$eigenval_1[k] <- best_model_output$eigval[1]
        output_param$eigenval_2[k] <- best_model_output$eigval[2]
        output_param$AIC[k] <- best_model_output$AIC
        output_param$AICc[k] <- best_model_output$AICc
        output_param$lnL[k] <- best_model_output$loglik
        output_param$BIC[k] <- best_model_output$BIC
        output_param$alpha_single[k] <- paste0(best_model_output$solution[1]) 
        output_param$sigma.sq_single[k] <- paste0(best_model_output$solution[2]) 
        output_param$theta_single[k] <- paste0(best_model_output$theta[1]) 
        output_param$rep[k] <- k
      } else {
        output_param$SE_theta_AM[k] <- best_model_output$theta[1,2]
        output_param$SE_alpha_AM[k] <- best_model_output$solution.se[1,1]
        output_param$SE_sigsq_AM[k] <- best_model_output$solution.se[2,1]
        output_param$SE_theta_EM[k] <- best_model_output$theta[2,2]
        output_param$SE_alpha_EM[k] <- best_model_output$solution.se[1,2]
        output_param$SE_sigsq_EM[k] <- best_model_output$solution.se[2,2]
        output_param$eigenval_1[k] <- best_model_output$eigval[1]
        output_param$eigenval_2[k] <- best_model_output$eigval[2]
        output_param$AIC[k] <- best_model_output$AIC
        output_param$AICc[k] <- best_model_output$AICc
        output_param$lnL[k] <- best_model_output$loglik
        output_param$BIC[k] <- best_model_output$BIC
        output_param$alpha_AM[k] <- best_model_output$solution[1,1]
        output_param$sigma.sq_AM[k] <- best_model_output$solution[2,1]
        output_param$theta_AM[k] <- best_model_output$solution[2,1]
        output_param$alpha_EM[k] <- best_model_output$solution[1,2]
        output_param$sigma.sq_EM[k] <- best_model_output$solution[2,2]
        output_param$theta_EM[k] <- best_model_output$solution[2,2]
        output_param$rep[k] <- k
        
      }
      
    }
    if (j == 1){
      output_aics <- output_aics_per
      output_params_pcs <- output_param
    } else {
      output_aics <-  rbind(output_aics, output_aics_per)
      output_params_pcs <- rbind(output_params_pcs, output_param)
    }
    
  }
  print(paste("tree", i))
  if (i == 1){
    output_aics_all <- output_aics
    output_params_all <- output_params_pcs
    
  } else {
    output_aics_all <- rbind(output_aics_all, output_aics)
    output_params_all <- rbind(output_params_all, output_params_pcs)
  }
}



saveRDS(output_params_all, "./analysis/pca_analysis/testing_replicability_ouwie_parameters_pca_redo.rds")
saveRDS(output_aics_all, "./analysis/pca_analysis/testing_replicability_ouwie_aics_pca_redo.rds")

output_params_all <- readRDS("./analysis/pca_analysis/testing_replicability_ouwie_parameters_pca_redo.rds")


#compare with original parameters

output_parameters_original$AM_theta[1] == output_params_all$theta_AM[1]


####################

long_output_parameters <- output_aics %>% 
  dplyr::select(iteration, PC_axis, model, theta_AM, theta_EM, SE_theta_AM, SE_theta_EM, eigenval_2) %>% 
  pivot_longer(cols = c(theta_AM, theta_EM)) %>% #, theta_single
  dplyr::filter(between(value, -10, 10))

long_output_parameters_ou1 <- output_aics %>% 
  dplyr::select(iteration, PC_axis, model, theta_single, SE_theta_single, eigenval_2) %>% 
  #pivot_longer(cols = c(theta_AM, theta_EM)) %>% #, theta_single
  dplyr::filter(between(theta_single, -10, 10))  

colnames(output_aics)

mean_theta_pc1_stable <- output_aics %>% 
  dplyr::select(iteration, PC_axis, model, theta_AM, theta_EM, SE_theta_AM, SE_theta_EM) %>% 
  dplyr::filter(PC_axis == 1) %>% 
  dplyr::filter(SE_theta_EM < 1) %>% 
  dplyr::summarise(mean_AM = mean(theta_AM), mean_EM = mean(theta_EM), sd_AM = sd(theta_AM), sd_EM = sd(theta_EM))

mean_theta_pc1

#seems like some outliers
output_aics %>% 
  dplyr::select(PC_axis, theta_AM, theta_EM) %>% 
  dplyr::filter(PC_axis == 3)

long_output_parameters <- pivot_longer(output_aics, cols = c(theta_AM, theta_EM, alpha_AM, alpha_EM))

#get data in format for plotting
no_outliers_ou1 <- long_output_parameters_ou1 %>% 
  #dplyr::filter(!model %in% "OU1") %>% 
  dplyr::filter(between(theta_single, -10, 10)) %>% 
  dplyr::filter(SE_theta_single < 1)

no_outliers_ou1$theta_single <- as.numeric(no_outliers_ou1$theta_single)

no_outliers <- long_output_parameters %>% 
  #dplyr::filter(!model %in% "OU1") %>% 
  dplyr::filter(between(value, -10, 10))%>% 
  dplyr::filter(SE_theta_EM < 1)

# only_outliers <- long_output_parameters %>% 
#   dplyr::filter(!between(value, -10, 10))

long_output_parameters %>% 
  dplyr::filter(model %in% "OUMVA") %>% 
  dplyr::select(PC_axis) %>% 
  unique()



#compare each theta (how different from different runs?)

jpeg("./output/PCAs/consensus_plots/distribution_theta_alpha_pc_ouwie_models_pc_123_myc_no_outliers_inc_ou1.jpg", height = 10, width = 14, units = "in", res = 600)
some_models <- ggplot(no_outliers, aes(model, value))+ #
  geom_boxplot()+
  xlab("Model")+
  ylab("Parameter value")+
  labs(title = "No outliers and no OU1 models")+
  facet_grid(PC_axis ~name, scales = "free")+#, labeller = as_labeller(model_names))+
  #geom_hline(yintercept = 100, colour = "red")+
  #geom_hline(yintercept = -100, colour = "red")+
  theme_bw()+
  theme(axis.ticks.x=element_blank())

some_models
dev.off()

mean_thetas_ou1 <- no_outliers_ou1 %>% 
  dplyr::select(-c(iteration, model, SE_theta_single, eigenval_2)) %>% 
  #dplyr::filter(!is.na(theta_single)) %>% 
  dplyr::group_by(PC_axis) %>% 
  dplyr::summarise(mean = mean(theta_single), sd = sd(theta_single))

mean_thetas <- no_outliers %>% 
  dplyr::select(-c(iteration, model, SE_theta_AM, SE_theta_EM, eigenval_2)) %>% 
  dplyr::filter(name %in% c("theta_AM", "theta_EM")) %>% 
  dplyr::group_by(PC_axis, name) %>% 
  dplyr::summarise(mean = mean(value), sd = sd(value))

mean_thetas

#get simplified model names
no_outliers$model_simp <- "OUM"
no_outliers$model_simp[which(no_outliers$model %in% "OU1")] <- "OU1"


# mean_thetas_inc_ou1 <- no_outliers %>% 
#   dplyr::select(-c(iteration, model, SE_theta_AM, SE_theta_EM, eigenval_2)) %>% 
#   dplyr::filter(name %in% c("theta_AM", "theta_EM")) %>% 
#   dplyr::group_by(PC_axis, name, model_simp) %>% 
#   dplyr::summarise(mean = mean(value), sd = sd(value))

consensus_tree <- readRDS("./data/for_analysis/myc_consensus_tree.rds")

#read in data
data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")
tree_myc <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

tips_to_drop <- consensus_tree$tip.label[-which(consensus_tree$tip.label %in% data_spectra$species)]

consensus_tree_pruned <- drop.tip(consensus_tree, tip = tips_to_drop)


consensus_tree_pruned$tip.label <- gsub("_", " ", consensus_tree_pruned$tip.label)
myc_named <- setNames(data_spectra$myc, rownames(data_spectra$spectra))

simmap_consensus <- make.simmap(consensus_tree_pruned, myc_named, model="SYM", nsim=100)
summary_simmap <-describe.simmap(simmap_consensus,plot=TRUE,cex=0.7)

simmap_consensus[[1]]$tip.label <- gsub(" ", "_", simmap_consensus[[1]]$tip.label)

combo <- readRDS("./analysis/pca_analysis/scores_for_all_reps_and_pcs.rds")

#get data as named list
pc1_all<-setNames(combo$V1,combo$species_factor)
names(pc1_all) <- gsub(" ", "_", names(pc1_all))

pc2_all<-setNames(combo$V2,combo$species_factor)
names(pc2_all) <- gsub(" ", "_", names(pc2_all))

pc3_all<-setNames(combo$V3,combo$species_factor)
names(pc3_all) <- gsub(" ", "_", names(pc3_all))

cols_tree<-setNames(c("blue","orange"),unique(myc_named))


ss<-getStates(simmap_consensus[[1]],"tips")
colors<-setNames(c("orange","blue"),c("AM","EM"))

boxcols<-setNames(sapply(ss,function(pc1_all,y) y[which(names(y)==pc1_all)],
                         y=colors),names(ss))

#sort colours to match tree and data
sorted_boxcols <- boxcols[order(match(names(boxcols),simmap_consensus[[1]]$tip.label))]

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

trans_orange <- t_col("orange", 90, name = "trans_orange")
trans_blue <- t_col("blue", 95, name = "trans_blue")
trans_green <- t_col("brown", 90, name = "trans_green")

#plot including ou1 for pc2
#make plot
jpeg("./output/PCAs/consensus_plots/consensus_myc_pca123_across_reps_boxplot_scores_with_optima_inc_ou1_excl_unstable.jpg", height = 9, width = 8, units = "in", res = 600)
#pdf("./output/PCAs/consensus_plots/consensus_myc_pca123_across_reps_boxplot_scores_with_optima_inc_ou1_excl_unstable.pdf", height = 9, width = 8)

#par(mfrow=c(1,4))
#layout(matrix(c(1,1,1,2,2,3,3,4,4), nrow = 1, ncol = 9, byrow = TRUE))
layout(matrix(c(1,1,1,1,1,2,2,3,3,4,4), nrow = 1, ncol = 11, byrow = TRUE))

plot(simmap_consensus[[1]],cols_tree,ftype='i', fsize = 0.6, xlim = c(0,500), mar=c(5.1,1.1,2.1,0.1))
par(mar=c(5.1,0.1,2.1,1.1))

legend(10, 15,   # Coordinates (x also accepts keywords)
       c("AM", "EM", "Shared"),#legend, # Vector with the name of each group
       #c("orange", "blue", "brown"),#fill,   # Creates boxes in the legend with the specified colors
       col = c("orange", "blue", "brown"),#par("col"), # Color of lines or symbols
       border = "black", # Fill box border color
       lty = c(1,2,4),
       #border = c("orange", "blue", "brown"),
       #lwd = 1,         # Line type and width
       #pch,              # Add pch symbols to legend lines or boxes
       bty = "n",        # Box type (bty = "n" removes the box)
       bg = par("bg"),    # Background color of the legend
       box.lwd = par("lwd"), # Legend box line width
       box.lty = par("lty"), # Legend box line type
       box.col = par("fg"),  # Legend box line color
       cex = 1.5,          # Legend size
       horiz = FALSE,     # Horizontal (TRUE) or vertical (FALSE) legend
       title = "Optima",      # Legend title
       title.adj = 0.1, 
       x.intersp = 0.5, 
       y.intersp = 1,
       xjust = 0
)

legend(10, 25,   # Coordinates (x also accepts keywords)
       c("AM", "EM"),#legend, # Vector with the name of each group
       c("orange", "blue"),#fill,   # Creates boxes in the legend with the specified colors
       #col = c("orange", "blue", "brown"),#par("col"), # Color of lines or symbols
       border = "black", # Fill box border color
       #lty = c(1,2,4),
       #border = c("orange", "blue", "brown"),
       #lwd = 1,         # Line type and width
       #pch,              # Add pch symbols to legend lines or boxes
       bty = "n",        # Box type (bty = "n" removes the box)
       bg = par("bg"),    # Background color of the legend
       box.lwd = par("lwd"), # Legend box line width
       box.lty = par("lty"), # Legend box line type
       box.col = par("fg"),  # Legend box line color
       cex = 1.5,          # Legend size
       horiz = FALSE,     # Horizontal (TRUE) or vertical (FALSE) legend
       title = "Scores",      # Legend title
       title.adj = 0.25, 
       x.intersp = 0.5, 
       y.intersp = 1,
       xjust = 0
)

boxplot(pc1_all~factor(names(pc1_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "", boxlwd = 0.5, medlwd = 1, whisklwd = 1, staplelwd = 1, outlwd = 0.5, outcex  = 0.5)
axis(1, cex.axis = 1.5)
abline(v = 0)
title(xlab="PC 1 scores", cex.lab = 1.5)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 1 & no_outliers$name == "theta_AM" & no_outliers$model_simp %in% "OUM")]), col = "orange", lty = 1)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 1 & no_outliers$name == "theta_EM" & no_outliers$model_simp %in% "OUM")]), col = "blue", lty = 2)
rect(xleft=(mean_thetas$mean[1]+mean_thetas$sd[1]), xright=(mean_thetas$mean[1]-mean_thetas$sd[1]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_orange, border = NA)
rect(xleft=(mean_thetas$mean[2]+mean_thetas$sd[2]), xright=(mean_thetas$mean[2]-mean_thetas$sd[2]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_blue, border = NA)


boxplot(pc2_all~factor(names(pc2_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "", boxlwd = 0.5, medlwd = 1, whisklwd = 1, staplelwd = 1, outlwd = 0.5, outcex  = 0.5)
axis(1, cex.axis = 1.5)
abline(v = 0)
title(xlab="PC 2 scores", cex.lab = 1.5)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 2 & no_outliers$name == "theta_AM" & no_outliers$model_simp %in% "OUM")]), col = "orange", lty = 1)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 2 & no_outliers$name == "theta_EM" & no_outliers$model_simp %in% "OUM")]), col = "blue", lty = 2)
abline(v = mean(no_outliers_ou1$theta_single[which(no_outliers_ou1$PC_axis == 2)]), col = "brown", lty = 4)
rect(xleft=(mean_thetas$mean[3]+mean_thetas_inc_ou1$sd[3]), xright=(mean_thetas_inc_ou1$mean[3]-mean_thetas_inc_ou1$sd[3]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_orange, border = NA)
rect(xleft=(mean_thetas$mean[4]+mean_thetas_inc_ou1$sd[4]), xright=(mean_thetas_inc_ou1$mean[4]-mean_thetas_inc_ou1$sd[4]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_blue, border = NA)
rect(xleft=(mean_thetas_ou1$mean[1]+mean_thetas_ou1$sd[1]), xright=(mean_thetas_ou1$mean[1]-mean_thetas_ou1$sd[1]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_green, border = NA)

boxplot(pc3_all~factor(names(pc3_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "", boxlwd = 0.5, medlwd = 1, whisklwd = 1, staplelwd = 1, outlwd = 0.5, outcex  = 0.5)
axis(1, cex.axis = 1.5)
abline(v = 0)
title(xlab="PC 3 scores", cex.lab = 1.5)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 3 & no_outliers$name == "theta_AM" & no_outliers$model_simp %in% "OUM")]), col = "orange", lty = 1)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 3 & no_outliers$name == "theta_EM" & no_outliers$model_simp %in% "OUM")]), col = "blue", lty = 2)
rect(xleft=(mean_thetas$mean[5]+mean_thetas$sd[5]), xright=(mean_thetas$mean[5]-mean_thetas$sd[5]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_orange, border = NA)
rect(xleft=(mean_thetas$mean[6]+mean_thetas$sd[6]), xright=(mean_thetas$mean[6]-mean_thetas$sd[6]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_blue, border = NA)


dev.off()

