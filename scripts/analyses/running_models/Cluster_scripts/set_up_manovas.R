#read libraries
library(phytools)
library(mvMORPH)

#for file in folder
files <- list.files("/home/jjantzen/scratch/R/output/myc_lp_models/92sp_binary", pattern="^myc_lp_model_.*rds", all.files=FALSE, full.names=TRUE)

#choose which models to use 
best_models <- readRDS("/home/jjantzen/scratch/R/input/iterations_myc/best_models_myc_lp_models_92sp_binary.rds")

best_models <- best_models[,c(1,2)]

#subset files to match best models only
keep_files <- c()

name_formula <- gsub("_BM_iteration1.rds", "", files[1])
name_formula

for (i in 1:100){
  file_name <- paste0(name_formula, "_", best_models$best_model[i], "_iteration", best_models$iteration[i], ".rds")  
  keep_files[i] <- file_name
  
}
#keep_files <- keep_files[-100]
keep_files

#read in models for keep files
model_list <- list()

for (i in 1:length(keep_files)){
  #read file
  model <- readRDS(keep_files[i])
  model_list[[i]] <- model
}

saveRDS(model_list, "/home/jjantzen/scratch/R/output/myc_lp_models/92sp_binary/best_myc_lp_models_for_manovas_92sp.rds")
