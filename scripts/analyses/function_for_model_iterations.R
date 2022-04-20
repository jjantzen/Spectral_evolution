#function for iterations of models

library(phytools)
library(mvMORPH)


#example expression
expression <- "myc * woody * persistence"

form <- as.formula(trait ~ woody * order)

#run model over number of trees and all evolutionary models

model_iterations <- function(expression, data_obj, trees, models){
  #create output dataframe
  df_output <- data.frame(matrix(nrow=100, ncol = 16))
  colnames(df_output) <- c("iteration", "BM.mserr", "BM.tuning", "OU.alpha", "OU.mserr", "OU.tuning", 
                           "EB.rate", "EB.mserr", "EB.tuning", "lambda.lambda", "lambda.mserr", "lambda.tuning",
                           "GIC_BM", "GIC_OU", "GIC_EB", "GIC_lambda")
  for(p in c(1:ncol(df_output))) {
    df_output[,p] <- as.numeric(df_output[,p])
  }
  #run iterations
  for (i in 1:length(trees)){
    #run for each of 4 models of evolution
    for (j in 1:length(models)){
      #specify model
      model <- mvgls(form, data=data_obj, tree=trees[i][[1]], model=models[j], error = TRUE)
      #assign output
      if (j == 1){
        df_output$iteration[i] <- i
        df_output$BM.mserr[i] <- model$mserr
        df_output$BM.tuning[i] <- model$tuning
        df_output$GIC_BM[i] <- GIC(model)$GIC
      } else if (j == 2){
        df_output$OU.alpha[i] <- model$param
        df_output$OU.mserr[i] <- model$mserr
        df_output$OU.tuning[i] <- model$tuning
        df_output$GIC_OU[i] <- GIC(model)$GIC
      } else if (j == 3){
        df_output$EB.rate[i] <- model$param
        df_output$EB.mserr[i] <- model$mserr
        df_output$EB.tuning[i] <- model$tuning
        df_output$GIC_EB[i] <- GIC(model)$GIC
      } else if (j == 4){
        df_output$lambda.lambda[i] <- model$param
        df_output$lambda.mserr[i] <- model$mserr
        df_output$lambda.tuning[i] <- model$tuning
        df_output$GIC_lambda[i] <- GIC(model)$GIC
      } 
      #save model itself
      saveRDS(model, paste0("./analysis/testing_models/", expression, "_model_", models[j], "_iteration", i, ".rds"))
    }
  }
  saveRDS(df_output, paste0("./analysis/testing_models/model_parameters_" , expression, ".rds"))
}

####
#test function for two trees and myc models

#read spectra
#data_spectra <- readRDS("./data/tidy/new_spectra_matched_trees.rds") #is this the right final object?

trait_data <- readRDS("./data/tidy/new_combo_data_matched_spectra.rds")

#read spectra
spectra_matrix <- readRDS("./data/tidy/spectra_not_reordered_to_tree.rds")

data_spectra <- list(trait=spectra_matrix, woody = trait_data$Woody, order = trait_data$Order)

str(data_spectra)

#prepare trees
new_trees <- readRDS("./data/tidy/new_trees_matched_spectra.rds")

trees <- new_trees[1:2]

#prepare set of models to run
models <- c("BM", "OU", "EB", "lambda")

#example expression
form <- as.formula(trait ~ woody * order)

#run model over number of trees and all evolutionary models

testing_version1 <- model_iterations(form, data_spectra, trees, models)

#looking at output from cluster
model_parameters_woody_x_order$lambda.lambda
model_parameters_woody_x_order$OU.alpha
woody_x_order_model_lambda_iteration2$param

#get total tree lengths
max(nodeHeights(trees[[1]])) #276.2251
max(nodeHeights(trees[[2]])) #152.6423


#do manova for woody and order
manova_I <- manova.gls(woody_x_order_model_lambda_iteration2, test = "Wilks", type = "I")
