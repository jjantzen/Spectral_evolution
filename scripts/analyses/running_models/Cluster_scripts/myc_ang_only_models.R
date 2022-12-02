#myc script

#read libraries
library(phytools)
library(mvMORPH)

taskid <- as.numeric(commandArgs(trailingOnly=TRUE))

#save object
#data_spectra <- readRDS("/home/jjantzen/scratch/R/input/iterations_myc/myc_data_list_92sp_binary_for_analysis.rds")
data_spectra <- readRDS("/home/jjantzen/scratch/R/input/iterations_myc/ang_only_data_for_myc.rds")

new_trees <- readRDS("/home/jjantzen/scratch/R/input/iterations_myc/ang_only_trees_for_myc.rds")

trees <- new_trees[[taskid]]
str(trees)
class(trees)

str(data_spectra)

#create output dataframe
df_output <- data.frame(matrix(nrow=4, ncol = 4))
colnames(df_output) <- c("iteration", "model", "GIC", "parameter")

#run for each of 4 models of evolution
models <- c("BM", "OU", "EB", "lambda")

for (j in 1:length(models)){
  j
  model <- mvgls(spectra ~ myc, data = data_spectra, tree = trees, model=models[j], error = TRUE)
  #assign output
  df_output$model[j] <- models[j]
  df_output$iteration[j] <- taskid
  df_output$parameter[j] <- model$param
  df_output$GIC[j] <- GIC(model)$GIC
  #save model itself
  saveRDS(model, paste0("/home/jjantzen/scratch/R/output/myc_models/ang_only/myc_model_ang_only_", models[j], "_iteration", taskid, ".rds"))
}

saveRDS(df_output, paste0("/home/jjantzen/scratch/R/output/myc_models/ang_only/model_parameters_myc_ang_only_iteration", taskid, ".rds"))


