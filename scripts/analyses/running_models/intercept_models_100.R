####code for running intercept models

library(phytools)
library(mvMORPH)

#prepare datasets

#read spectra
data_spectra <- readRDS("./data/tidy/new_spectra_matched_trees.rds") #is this the right final object?

#read trees
trees <- new_trees #(read them in - decide on rescaled trees or original height or both)

#prepare set of models to run
models <- c("BM", "OU", "EB", "lambda")

#iterate over 100 trees

#create output dataframe
df_output <- data.frame(matrix(nrow=100, ncol = 16))
colnames(df_output) <- c("iteration", "BM.mserr", "BM.tuning", "OU.alpha", "OU.mserr", "OU.tuning", 
                         "EB.rate", "EB.mserr", "EB.tuning", "lambda.lambda", "lambda.mserr", "lambda.tuning",
                         "GIC_BM", "GIC_OU", "GIC_EB", "GIC_lambda")

class(df_output[,c(1:16)]) <- "numeric"

#run loop over 100 trees
for (i in 1:length(trees)){
  #run for each of 4 models of evolution
  for (j in 1:length(models)){
    model <- mvgls(trait ~ 1, data=data_spectra, tree=trees[i], model=models[j], error = TRUE)
    #assign output
    if (j == 1){
      df_output$iteration[i] <- i
      df_output$BM.mserr[i] <- model$mserr
      df_output$BM.tuning[i] <- model$tuning
      df_output$GIC_BM[i] <- GIC(model)$GIC
    } if (j == 2){
      df_output$OU.alpha[i] <- model$param
      df_output$OU.mserr[i] <- model$mserr
      df_output$OU.tuning[i] <- model$tuning
      df_output$GIC_OU[i] <- GIC(model)$GIC
    } if (j == 3){
      df_output$EB.rate[i] <- model$param
      df_output$EB.mserr[i] <- model$mserr
      df_output$EB.tuning[i] <- model$tuning
      df_output$GIC_EB[i] <- GIC(model)$GIC
    } if (j == 4){
      df_output$lambda.lambda[i] <- model$param
      df_output$lambda.mserr[i] <- model$mserr
      df_output$lambda.tuning[i] <- model$tuning
      df_output$GIC_lambda[i] <- GIC(model)$GIC
    } 
    #save model itself
    saveRDS(model, paste0("./models/intercept_model_", models[j], "_iteration", i, ".rds"))
  }
}

saveRDS(df_output, "./model_parameters_intercept.rds")

###