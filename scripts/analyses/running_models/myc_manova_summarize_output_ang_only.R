#read output from myc models and summarize parameters
library(dplyr)
library(mvMORPH)

#for file in folder
files <- list.files("./analysis/myc_models/92sp_myc/ang_only/manovas/", pattern="*manova_summary*", all.files=FALSE, full.names=TRUE, include.dirs = FALSE)

for (i in 1:length(files)){
  #read file
  parameters <- readRDS(files[i])
  #take output values and put into single dataframe
  if (i == 1){
    output <- parameters
  } else {
    output <- rbind(output,parameters)
  }
}

example <- readRDS(files[1])
example

#sort output dataframe
output_sorted <- output[order(output$iteration),]

#save output
saveRDS(output_sorted, "./analysis/myc_models/92sp_myc/ang_only/summarized_manovas_myc_models_ang_only.rds")

#output_sorted <- readRDS("./analysis/myc_models/92sp_myc/summarized_manovas_myc_models_92sp_binary.rds")


#get mean of pvalue 
mean(output_sorted$pvalue)
sd(output_sorted$pvalue)
mean(output_sorted$test_stat)
sd(output_sorted$test_stat)

print(manova_output_iteration15_OU_model_)

nrow(output_sorted[which((output_sorted$pvalue < 0.05) == TRUE),])


##################if needed#############
#for each iteration, get best model (model for lowest GIC)

best_models <- data.frame(matrix(nrow=100, ncol = 3))
colnames(best_models)[1] <- "iteration"
class(best_models[,1]) <- "numeric"
colnames(best_models)[2] <- "best_model"
class(best_models[,2]) <- "character"
colnames(best_models)[3] <- "parameter"
class(best_models[,3]) <- "numeric"

for (i in 1:length(unique(output_sorted$iteration))){
  data <- output_sorted[which(output_sorted$iteration == i),]
  best_model <- data$model[which(data$GIC == min(data$GIC))]
  parameter <- data$parameter[which(data$GIC == min(data$GIC))]
  best_models$iteration[i] <- i
  best_models$best_model[i] <- best_model
  best_models$parameter[i] <- parameter
}

saveRDS(best_models, "./analysis/myc_models/best_models_myc_models.rds")

#summarize by number of model
summary <- best_models %>% 
  group_by(best_model) %>% 
  summarize(n = n(), parameter = mean(parameter))

saveRDS(output_sorted, "./analysis/myc_models/summary_myc_models.rds")


#get stats for results table
stats <- output_sorted %>% 
  group_by(model) %>% 
  summarize(mean_param = mean(parameter), sd2_param = sd(parameter), mean_GIC =  mean(GIC), sd2_GIC = sd(GIC))#, deltaGIC = (mean(GIC)-minGIC_EB))


#calculate phylogenetic halflife for each tree for the ou models
best_models

trees_ou <- best_models$iteration[which(best_models$best_model == "OU")]

trees_eb <- best_models$iteration[which(best_models$best_model == "EB")]

#get trees
trees_all <- readRDS("./data/for_analysis/final_trees_matched_spectra.rds")

ou_trees <- trees_all[trees_ou]
eb_trees <- trees_all[trees_eb]

length(eb_trees)

#get mean max height for these trees
best_models_ou <- best_models[which(best_models$best_model == "OU"),]

for (i in 1:length(ou_trees)){
  best_models_ou$tree_height[i] <- max(nodeHeights(ou_trees[[i]]))
  best_models_ou$halflife[i] <- log(2)/best_models_ou$parameter[i]
}

mean(best_models_ou$tree_height)
sd(best_models_ou$tree_height)

mean(best_models_ou$halflife)
sd(best_models_ou$halflife)

#get mean max height for these trees
best_models_eb <- best_models[which(best_models$best_model == "EB"),]

for (i in 1:length(eb_trees)){
  best_models_eb$tree_height[i] <- max(nodeHeights(eb_trees[[i]]))
  #$halflife[i] <- log(2)/best_models_ou$parameter[i]
}

mean(best_models_eb$tree_height)
sd(best_models_eb$tree_height)

# max(nodeHeights(ou_trees[[1]]))
# 
# halflife <- log(2)/alpha #in units of branch lengths 

##############################
deltaGIC <- GIC - min(GIC)

minGIC_OU <- -2434426-117713
minGIC_EB <- -2560278-16510


-2293811 - 10779
-2293811 + 10779

-2560278 - 16510
-2560278 + 16510

-2434426 - 117713
-2434426 + 117713

-2392054 - 27590
-2392054 + 27590


#compare with Artuso paper
(-40009-70) - (-40056.9-72) #48.9 #3.37–131.60

#min
-40009 - 70 #40079
#max 
-40009 + 70 #-39939

#lower lower
-40079 - (-40056.9 - 72) #49.9

#lower higher
-40079 - (-40056.9 + 72) #-94.1

#higher higher
-39939 - (-40056.9 - 72) #189.9

#higher lower
-39939 - (-40056.9 + 72) # 45.9

3.37 + (-40056.9 -72) #40132.27
3.37 - (-40056.9 +72) #39988.27

(-40009+70) - (-40056.9-72) #118.9 #3.37–131.63


