#read output from myc models and summarize parameters
library(dplyr)
library(mvMORPH)
library(stringr)


#for file in folder
#files <- list.files("./analysis/myc_gf_models/manovas/", pattern="*manova_summary*", all.files=FALSE, full.names=TRUE, include.dirs = FALSE)

files <- list.files("./analysis/myc_lp_interaction_models/manovas/", pattern="*manova_output*", all.files=FALSE, full.names=TRUE, include.dirs = FALSE)

files <- list.files("./analysis/myc_gf_models/manovas/", pattern="*manova_output*", all.files=FALSE, full.names=TRUE, include.dirs = FALSE)

files <- files[-4]

df_output <- data.frame(matrix(nrow=100, ncol = 8))
colnames(df_output) <- c("iteration", "model", "pvalue_myc", "pvalue_lp", "pvalue_myc:lp", "test_stat_myc", "test_stat_lp", "test_stat_myc:lp")

for (i in 1:length(files)){
  #read file
  parameters <- readRDS(files[i])
  #take output values and put into single dataframe
  
  item <- files[i]
  item <- gsub(paste0("./analysis/myc_lp_interaction_models/manovas/manova_output_iteration"), "", item)
  item <- gsub("_model_.rds", "", item)
  model <- word(item, 2, sep = "_")
  iteration <- word(item, 1, sep = "_")
  
  df_output$iteration[i] <- iteration
  
  df_output$model[i] <- model
  df_output$pvalue_myc[i] <- parameters$pvalue[1]
  df_output$pvalue_lp[i] <- parameters$pvalue[2]
  df_output$`pvalue_myc:lp`[i] <- parameters$pvalue[3]
  
  df_output$test_stat_myc[i] <- parameters$stat[1]
  df_output$test_stat_lp[i] <- parameters$stat[2]
  df_output$`test_stat_myc:lp`[i] <- parameters$stat[3]
}


missing_reps <- c(1:100)[-which(c(1:100) %in% df_output$iteration)]




example <- readRDS(files[15])
example

#sort output dataframe
output_sorted <- df_output[order(as.numeric(df_output$iteration)),]

#save output
saveRDS(df_output, "./analysis/myc_lp_interaction_models/summarized_manovas_myc_lp_models_92sp_binary_interaction_terms.rds")

df_output <- readRDS("./analysis/myc_gf_models/summarized_manovas_myc_lp_models_92sp_binary_interaction_terms.rds")


#get mean of pvalue 
min(output_sorted$pvalue_myc[which(!is.na(output_sorted$pvalue_myc))])

min(output_sorted$pvalue_lp[which(!is.na(output_sorted$pvalue_lp))])

min(output_sorted$`pvalue_myc:lp`[which(!is.na(output_sorted$`pvalue_myc:lp`))])

mean(output_sorted$pvalue_myc[which(!is.na(output_sorted$pvalue_myc))])
sd(output_sorted$pvalue_myc[which(!is.na(output_sorted$pvalue_myc))])

mean(output_sorted$pvalue_lp[which(!is.na(output_sorted$pvalue_lp))])
sd(output_sorted$pvalue_lp[which(!is.na(output_sorted$pvalue_lp))])

mean(output_sorted$`pvalue_myc:lp`[which(!is.na(output_sorted$`pvalue_myc:lp`))])
sd(output_sorted$`pvalue_myc:lp`[which(!is.na(output_sorted$`pvalue_myc:lp`))])


mean(output_sorted$test_stat_myc[which(!is.na(output_sorted$test_stat_myc))])
sd(output_sorted$test_stat_myc[which(!is.na(output_sorted$test_stat_myc))])

mean(output_sorted$test_stat_lp[which(!is.na(output_sorted$test_stat_lp))])
sd(output_sorted$test_stat_lp[which(!is.na(output_sorted$test_stat_lp))])

mean(output_sorted$`test_stat_myc:lp`[which(!is.na(output_sorted$`test_stat_myc:lp`))])
sd(output_sorted$`test_stat_myc:lp`[which(!is.na(output_sorted$`test_stat_myc:lp`))])






output_sorted_by_pvalue <- output[order(output$pvalue),]

mean(output_sorted$pvalue)
sd(output_sorted$pvalue)
mean(output_sorted$test_stat)
sd(output_sorted$test_stat)

print(manova_output_iteration79_OU_model_)
