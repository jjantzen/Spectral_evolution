#myc script
#read libraries
library(phytools)
library(mvMORPH)

#get task id
taskid <- as.numeric(commandArgs(trailingOnly=TRUE))
taskid
#read list of models
model_list <- readRDS("/home/jjantzen/scratch/R/output/myc_models/92sp_binary/best_myc_models_for_manovas_92sp.rds")

#choose which models to use 
best_models <- readRDS("/home/jjantzen/scratch/R/input/iterations_myc/best_models_myc_models_92sp.rds")
best_models <- best_models[,c(1,2)]
best_models
###manova
#create output dataframe
df_output <- data.frame(matrix(nrow=1, ncol = 4))
colnames(df_output) <- c("iteration", "model", "pvalue", "test_stat")

#manova
manova_output <- manova.gls(model_list[[taskid]], test = "Wilks", type = "II", nperm = 1000, nbcores = 4L) #need to request 4 cores from cluster
df_output$iteration[1] <- taskid
df_output$model[1] <- best_models$best_model[taskid]
df_output$pvalue[1] <- manova_output$pvalue
df_output$test_stat[1] <- manova_output$stat
saveRDS(manova_output, paste0("/home/jjantzen/scratch/R/output/myc_models/92sp_binary/manova_output_iteration", taskid, "_", best_models$best_model[taskid], "_model_.rds"))
  
saveRDS(df_output, paste0("/home/jjantzen/scratch/R/output/myc_models/92sp_binary/manova_summary_dataframe_", taskid, ".rds"))

