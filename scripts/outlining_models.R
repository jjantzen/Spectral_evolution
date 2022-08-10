#Outlining models to run using trait data

#prepare datasets

#make lists of parameters
traits <- c("woody", "shade", "drought", "AM", "EM", "ErM", "NM", xxx) #add others as they are chosen
#include leaf persistence (evergreen vs deciduous), more detailed growth form (raunkauer?), pH, soil texture, habitat description

#prepare trees
trees <- new_trees #(read them in - decide on rescaled trees or original height or both)

#prepare set of models to run
models <- c("BM", "OU", "EB", "lambda")

#outline of models to run
#First, run the intercept models on the spectral data alone

#Second, run the next simplest versions of the models (single predictor)

#make for loop to loop over different models for each tree and each variable

model_fit <- mvgls(trait ~ traits[k], data=data_spectra, tree=trees[i], model=models[j], error = TRUE) 

#Third, run the more complicated models including interactions

#make for loop to loop over different models for each tree

#probably need to specify the predictors individually rather than indexed
model_fit <- mvgls(trait ~ traits[k] + traits[k+1], data=data_spectra, tree=trees[i], model=models[j], error = TRUE) 


####Myc hypotheses

#simplest models first?
model_fit <- mvgls(trait ~ myc, data=data_spectra, tree=trees[i], model=models[j], error = TRUE) 

model_fit <- mvgls(trait ~ woody, data=data_spectra, tree=trees[i], model=models[j], error = TRUE) 

model_fit <- mvgls(trait ~ lp, data=data_spectra, tree=trees[i], model=models[j], error = TRUE) 

model_fit <- mvgls(trait ~ pH, data=data_spectra, tree=trees[i], model=models[j], error = TRUE)

model_fit <- mvgls(trait ~ hd, data=data_spectra, tree=trees[i], model=models[j], error = TRUE)

#additive models
model_fit <- mvgls(trait ~ myc + woody + lp, data=data_spectra, tree=trees[i], model=models[j], error = TRUE) 
model_fit <- mvgls(trait ~ myc*wooody + myc*lp, data=data_spectra, tree=trees[i], model=models[j], error = TRUE) 

#interaction models
model_fit <- mvgls(trait ~ myc*woody*lp, data=data_spectra, tree=trees[i], model=models[j], error = TRUE) 

#non myc explanation
model_fit <- mvgls(trait ~ pH*woody + pH*lp, data=data_spectra, tree=trees[i], model=models[j], error = TRUE) 

model_fit <- mvgls(trait ~ myc, data=data_spectra, tree=trees[i], model=models[j], error = TRUE) 


####code for running intercept models
#iterate over 100 trees

#create output dataframe
df_output <- data.frame(matrix(nrow=100, ncol = 12))
colnames(df_output) <- c("iteration", "BM.mserr", "BM.tuning", "OU.alpha", "OU.mserr", "OU.tuning", 
                         "EB.rate", "EB.mserr", "EB.tuning", "lambda.lambda", "lambda.mserr", "lambda.tuning")

class(df_output[,c(1:12)]) <- "numeric"

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
    } if (j == 2){
      df_output$OU.alpha[i] <- model$param
      df_output$OU.mserr[i] <- model$mserr
      df_output$OU.tuning[i] <- model$tuning
    } if (j == 3){
      df_output$EB.rate[i] <- model$param
      df_output$EB.mserr[i] <- model$mserr
      df_output$EB.tuning[i] <- model$tuning
    } if (j == 4){
      df_output$lambda.lambda[i] <- model$param
      df_output$lambda.mserr[i] <- model$mserr
      df_output$lambda.tuning[i] <- model$tuning
    } 
    #save model itself
    saveRDS(model, paste0("./models/intercept_model_", models[j], "_iteration", i, ".rds"))
  }
}

saveRDS(df_output, "./model_parameters_intercept.rds")

###

models <- c("BM", "OU", "EB", "lambda")

fit_bm_intercept$mserr
fit_eb_intercept$param
fit_ou_intercept$param
fit_lambda_intercept$param
