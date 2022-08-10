#read output from myc models and summarize parameters
library(dplyr)
library(phytools)
library(ggplot2)
library(ggpubr)

#for file in folder
files <- list.files("./analysis/myc_models/92sp_myc/ang_only/parameters/", pattern="model_parameters", all.files=FALSE, full.names=TRUE)

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

#sort output dataframe
output_sorted <- output[order(output$iteration),]

#save output
saveRDS(output_sorted, "./analysis/myc_models/92sp_myc/ang_only/summarized_parameters_myc_models_ang_only.rds")

output_sorted <- readRDS("./analysis/myc_models/92sp_myc/ang_only/summarized_parameters_myc_models_ang_only.rds")

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


saveRDS(best_models, "./analysis/myc_models/92sp_myc/ang_only/best_models_myc_models_ang_only.rds")
best_models <- readRDS("./analysis/myc_models/92sp_myc/ang_only/best_models_myc_models_ang_only.rds")

#summarize by number of model
summary <- best_models %>% 
  group_by(best_model) %>% 
  summarize(n = n(), parameter = mean(parameter))

saveRDS(summary, "./analysis/myc_models/92sp_myc/ang_only/summary_myc_models_ang_only.rds")

summary <- readRDS("./analysis/myc_models/92sp_myc/ang_only/summary_myc_models_ang_only.rds")


#get stats for results table
stats <- output_sorted %>% 
  dplyr::group_by(model) %>% 
  dplyr::summarize(mean_param = mean(parameter), sd2_param = sd(parameter), mean_GIC =  mean(GIC), sd2_GIC = sd(GIC))#, deltaGIC = (mean(GIC)-minGIC_EB))

saveRDS(stats, "./analysis/myc_models/92sp_myc/ang_only/stats_myc_models_ang_only.rds")

stats <- readRDS("./analysis/myc_models/92sp_myc/ang_only/stats_myc_models_ang_only.rds")

#calculate phylogenetic halflife for each tree for the ou models
best_models

trees_ou <- best_models$iteration[which(best_models$best_model == "OU")]

trees_eb <- best_models$iteration[which(best_models$best_model == "EB")]

trees_bm <- best_models$iteration[which(best_models$best_model == "BM")]

trees_lambda <- best_models$iteration[which(best_models$best_model == "lambda")]

#get trees
trees_all <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

ou_trees <- trees_all[trees_ou]
eb_trees <- trees_all[trees_eb]
bm_trees <- trees_all[trees_bm]
lambda_trees <- trees_all[trees_lambda]

length(eb_trees)

#get mean max height for these trees
best_models_ou <- best_models[which(best_models$best_model == "OU"),]

for (i in 1:length(ou_trees)){
  best_models_ou$tree_height[i] <- max(nodeHeights(ou_trees[[i]]))
  best_models_ou$halflife[i] <- log(2)/best_models_ou$parameter[i]
}

mean(best_models_ou$tree_height[-c(2)])
sd(best_models_ou$tree_height[-c(2)])
nrow(best_models_ou)
mean(best_models_ou$halflife[-c(2)])
sd(best_models_ou$halflife[-c(2)])

mean(best_models_ou$parameter[-c(2)])
sd(best_models_ou$parameter[-c(2)])


#get mean max height for these trees
best_models_eb <- best_models[which(best_models$best_model == "EB"),]

for (i in 1:length(eb_trees)){
  best_models_eb$tree_height[i] <- max(nodeHeights(eb_trees[[i]]))
  #$halflife[i] <- log(2)/best_models_ou$parameter[i]
}

nrow(best_models_eb)

mean(best_models_eb$tree_height)
sd(best_models_eb$tree_height)

mean(best_models_eb$parameter[-which(best_models_eb$iteration %in% c(98))])
sd(best_models_eb$parameter[-which(best_models_eb$iteration %in% c(98))])




best_models_bm <- best_models[which(best_models$best_model == "BM"),]

for (i in 1:length(bm_trees)){
  best_models_bm$tree_height[i] <- max(nodeHeights(bm_trees[[i]]))
  #$halflife[i] <- log(2)/best_models_ou$parameter[i]
}

mean(best_models_bm$tree_height)
sd(best_models_bm$tree_height)

best_models_lambda <- best_models[which(best_models$best_model == "lambda"),]

for (i in 1:length(lambda_trees)){
  best_models_lambda$tree_height[i] <- max(nodeHeights(lambda_trees[[i]]))
  #$halflife[i] <- log(2)/best_models_ou$parameter[i]
}

outliers <- c(75, 82, 8, 3)
mean(best_models_lambda$tree_height)#[-which(best_models_lambda$iteration %in% outliers)])
sd(best_models_lambda$tree_height)#[-which(best_models_lambda$iteration %in% outliers)])

mean(best_models_lambda$parameter)

mean(best_models_lambda$parameter[-which(best_models_lambda$iteration %in% outliers)])
sd(best_models_lambda$parameter[-which(best_models_lambda$iteration %in% outliers)])



#plot distribution of parameters

model_names <- c(
  "rate" = "Early Burst",
  "lambda" = "Lambda",
  "alpha" = "Ornstein-Uhlenbeck"
)

best_models$best_model <- gsub("EB", "rate", best_models$best_model)
best_models$best_model <- gsub("OU", "alpha", best_models$best_model)

output_sorted$model <- gsub("OU", "alpha", output_sorted$model)
output_sorted$model <- gsub("EB", "rate", output_sorted$model)


all_models <- ggplot(output_sorted[-which(output_sorted$model == "BM"),], aes(model, parameter))+
  geom_boxplot()+
  xlab("Model")+
  ylab("Parameter value")+
  labs(title = "All models")+
  facet_wrap(~model, scales = "free", labeller = as_labeller(model_names))+
  theme_bw()+
  theme(axis.ticks.x=element_blank())

best_only

best_only <- ggplot(best_models[which(best_models$parameter != "BM" & !is.na(best_models$parameter)),], aes(best_model, parameter))+ #[-which(best_models$best_model %in% c("BM", "lambda")),]
  geom_boxplot()+
  xlab("Model")+
  ylab("Parameter value")+
  labs(title = "Best models")+
  facet_wrap(~best_model, scales = "free", labeller = as_labeller(model_names))+
  theme_bw()+
  theme( axis.ticks.x=element_blank())


pdf("./output/model_assessment/myc_ang_only_parameter_distributions_binary.pdf", width = 10, height = 8)
par(mfrow=c(2,1))
best_only
all_models
dev.off()

jpeg("./output/model_assessment/myc_ang_only_parameter_distributions_binary.jpg", width = 6, height = 6, units = "in", res = 500)
ggarrange(best_only, all_models, nrow = 2)
# best_only
# all_models
dev.off()

########
keep_files <- c()

name_formula <- gsub("_BM_iteration1.rds", "", files[1])
name_formula

for (i in 1:100){#unique(best_models$iteration)){
  file_name <- paste0(name_formula, "_", best_models$best_model[i], "_iteration", best_models$iteration[i], ".rds")  
  keep_files[i] <- file_name
  
}

keep_files <- keep_files[-100]
