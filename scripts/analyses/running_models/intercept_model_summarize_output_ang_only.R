#read output from myc models and summarize parameters
library(dplyr)
library(tidyr)
library(mvMORPH)
#library(cowplot)
library(ggpubr)

#for file in folder
files <- list.files("./analysis/intercept_models/myc_dataset/ang_only/parameters/", pattern=NULL, all.files=FALSE, full.names=TRUE)

for (i in 1:length(files)){
  #read file
  parameters <- readRDS(files[i])
  #take output values and put into single dataframe
  if (i == 1){
    output <- parameters[1,]
  } else {
    output <- rbind(output,parameters[1,])
  }
}


strsplit(files[[1]], "_")

#sort output dataframe
output_sorted <- output[order(output$iteration),]
nrow(output_sorted)

#save output
saveRDS(output_sorted, "./analysis/intercept_models/myc_dataset/ang_only/summarized_parameters_intercept_models_ang_only.rds")

nrow(output_sorted)

# #for each iteration, get best model (model for lowest GIC)
# 
# best_models <- data.frame(matrix(nrow=100, ncol = 3))
# colnames(best_models)[1] <- "iteration"
# class(best_models[,1]) <- "numeric"
# colnames(best_models)[2] <- "best_model"
# class(best_models[,2]) <- "character"
# colnames(best_models)[3] <- "parameter"
# class(best_models[,3]) <- "numeric"
# 
# for (i in 1:length(unique(output_sorted$iteration))){
#   data <- output_sorted[which(output_sorted$iteration == i),]
#   best_model <- data$model[which(data$GIC == min(data$GIC))]
#   parameter <- data$parameter[which(data$GIC == min(data$GIC))]
#   best_models$iteration[i] <- i
#   best_models$best_model[i] <- best_model
#   best_models$parameter[i] <- parameter
# }
# 
# saveRDS(best_models, "./analysis/intercept_models/best_models_intercept_models.rds")
# 
# #summarize by number of model
# summary <- best_models %>% 
#   group_by(best_model) %>% 
#   summarize(n = n(), parameter = mean(parameter))
# 
# saveRDS(summary, "./analysis/intercept_models/summary_intercept_models.rds")

#read summary with gics
model_parameters_gics <- output_sorted

#model_parameters_gics <- readRDS("./analysis/intercept_models/model_parameters_gics.rds")

nrow(model_parameters_gics)
tail(model_parameters_gics)
model_parameters_gics$BM.mserr

#trim to only filled rows
#model_parameters_gics_filled <- model_parameters_gics[which(complete.cases(model_parameters_gics)),]
model_parameters_gics_filled <- model_parameters_gics


#get best models for each row
for (i in 1:nrow(model_parameters_gics_filled)){
  min_val <- min(model_parameters_gics_filled[i,c(13:16)])
  model_parameters_gics_filled$best_model[i] <- colnames(model_parameters_gics_filled)[which(model_parameters_gics_filled[i,] == min_val)]
  model_parameters_gics_filled$best_model[i] <- gsub("GIC_", "", model_parameters_gics_filled$best_model[i])
}
     
saveRDS(model_parameters_gics_filled, "./analysis/intercept_models/myc_dataset/ang_only/summarized_best_models_100trees_mycdataset.rds")

#get summary stats overall
summary <- model_parameters_gics_filled %>% 
  group_by(best_model) %>% 
  summarize(n = n(), EB.rate = mean(EB.rate), OU.alpha = mean(OU.alpha), lambda = mean(lambda.lambda))

saveRDS(summary, "./analysis/intercept_models/myc_dataset/ang_only/summary_intercept_models_mycdataset.rds")

#calculate phylogenetic halflife for each tree for the ou models
lambda_models <- model_parameters_gics_filled$iteration[which(model_parameters_gics_filled$best_model == "lambda")]

ou_models <- model_parameters_gics_filled$iteration[which(model_parameters_gics_filled$best_model == "OU")]

eb_models <- model_parameters_gics_filled$iteration[which(model_parameters_gics_filled$best_model == "EB")]

bm_models <- model_parameters_gics_filled$iteration[which(model_parameters_gics_filled$best_model == "BM")]

#get trees
#trees_all <- readRDS("./data/for_analysis/final_trees_matched_spectra.rds")
trees_all <- readRDS("./data/for_analysis/ang_only_trees_for_myc.rds")

lambda_trees <- trees_all[lambda_models]
ou_trees <- trees_all[ou_models]
eb_trees <- trees_all[eb_models]
bm_trees <- trees_all[bm_models]

length(eb_trees)

#get mean max height for these trees
model_parameters_gics_filled$halflife <- NA
model_parameters_gics_filled$tree_height <- NA

for (i in 1:length(ou_trees)){
  iteration <- ou_models[i]
  model_parameters_gics_filled$tree_height[which(model_parameters_gics_filled$iteration == iteration)] <- max(nodeHeights(ou_trees[[i]]))
  model_parameters_gics_filled$halflife[which(model_parameters_gics_filled$iteration == iteration)] <- log(2)/model_parameters_gics_filled$OU.alpha[i]
}

for (i in 1:length(eb_trees)){
  iteration <- eb_models[i]
  model_parameters_gics_filled$tree_height[which(model_parameters_gics_filled$iteration == iteration)] <- max(nodeHeights(eb_trees[[i]]))
}

for (i in 1:length(bm_trees)){
  iteration <- bm_models[i]
  model_parameters_gics_filled$tree_height[which(model_parameters_gics_filled$iteration == iteration)] <- max(nodeHeights(bm_trees[[i]]))
}

for (i in 1:length(lambda_trees)){
  iteration <- lambda_models[i]
  model_parameters_gics_filled$tree_height[which(model_parameters_gics_filled$iteration == iteration)] <- max(nodeHeights(lambda_trees[[i]]))
}

#saveRDS(model_parameters_gics_filled, "./analysis/intercept_models/myc_dataset/ang_only/summarized_best_models_100trees_mycdataset.rds")

model_parameters_gics_filled <- readRDS("./analysis/intercept_models/myc_dataset/ang_only/summarized_best_models_100trees_angonly.rds")

#get stats by model
stats_BM <- model_parameters_gics_filled[which(model_parameters_gics_filled$best_model == "BM"),] %>% 
  group_by(best_model) %>% 
  summarize(mean_GIC =  mean(GIC_BM), sd2_GIC = sd(GIC_BM), mean_tree_height = mean(tree_height), sd_tree_height = sd(tree_height))

stats_EB <- model_parameters_gics_filled[which(model_parameters_gics_filled$best_model == "EB"),] %>% 
  group_by(best_model) %>% 
  summarize(mean_r = mean(EB.rate), sd2_r = sd(EB.rate), mean_GIC =  mean(GIC_EB), sd2_GIC = sd(GIC_EB), mean_tree_height = mean(tree_height), sd_tree_height = sd(tree_height))

stats_OU <- model_parameters_gics_filled[which(model_parameters_gics_filled$best_model == "OU"),] %>% 
  group_by(best_model) %>% 
  summarize(mean_alpha = mean(OU.alpha), sd2_param = sd(OU.alpha), mean_GIC =  mean(GIC_OU), sd2_GIC = sd(GIC_OU), mean_tree_height = mean(tree_height), sd_tree_height = sd(tree_height), mean_halflife = mean(halflife), sd_halflife = sd(halflife))

stats_lambda <- model_parameters_gics_filled[which(model_parameters_gics_filled$best_model == "lambda"),] %>% 
  group_by(best_model) %>% 
  summarize(mean_lambda = mean(lambda.lambda), sd2_lambda = sd(lambda.lambda), mean_GIC =  mean(GIC_lambda), sd2_GIC = sd(GIC_lambda), mean_tree_height = mean(tree_height), sd_tree_height = sd(tree_height))

#put this into results presentation

#need to recalculate mean GICs over all trees/model type (instead of just for trees where model is best)

summarized_parameters_intercept_models_mycdataset
summarized_best_models_100trees_mycdataset
summary_intercept_models_mycdataset

#make long dataset

long_all_trees <- pivot_longer(model_parameters_gics_filled[c(1,4,7,10)], cols = c(2:4))

#get "best models"

BM_list <- pivot_longer(model_parameters_gics_filled[which(model_parameters_gics_filled$best_model == "BM"),c(1,17:19)], cols = c(2))

OU_list <- pivot_longer(model_parameters_gics_filled[which(model_parameters_gics_filled$best_model == "OU"),c(1,4,17:19)], cols = c(2))

EB_list <- pivot_longer(model_parameters_gics_filled[which(model_parameters_gics_filled$best_model == "EB"),c(1,7,17:19)], cols = c(2))

lambda_list <- pivot_longer(model_parameters_gics_filled[which(model_parameters_gics_filled$best_model == "lambda"),c(1,10,17:19)], cols = c(2))

best_model_list <- rbind(OU_list, EB_list, lambda_list)

nrow(best_model_list)

colnames(model_parameters_gics_filled)

#figure out how to plot distribution

long_all_trees$model <- long_all_trees$name

model_names <- c(
  "EB" = "Early Burst",
  "lambda" = "Lambda",
  "OU" = "Ornstein-Uhlenbeck"
)

model_names_all <- c(
  "EB.rate" = "Early Burst",
  "lambda.lambda" = "Lambda",
  "OU.alpha" = "Ornstein-Uhlenbeck"
)

names(model_names) <- c("EB", "lambda", "OU")

unique(long_all_trees$name)

long_all_trees$name <- gsub("OU.alpha", "alpha", long_all_trees$name)
long_all_trees$name <- gsub("EB.rate", "rate", long_all_trees$name)
long_all_trees$name <- gsub("lambda.lambda", "lambda", long_all_trees$name)

all_models <- ggplot(long_all_trees, aes(name, value))+
  geom_boxplot()+
  xlab("Model")+
  ylab("Parameter value")+
  labs(title = "All models")+
  facet_wrap(~model, scales = "free", labeller = as_labeller(model_names_all))+
  theme_bw()+
  theme(axis.ticks.x=element_blank())

best_model_list$name <- gsub("OU.alpha", "alpha", best_model_list$name)
best_model_list$name <- gsub("EB.rate", "rate", best_model_list$name)
best_model_list$name <- gsub("lambda.lambda", "lambda", best_model_list$name)


best_only <- ggplot(best_model_list, aes(name, value))+
  geom_boxplot()+
  xlab("Model")+
  ylab("Parameter value")+
  labs(title = "Best models")+
  facet_wrap(~best_model, scales = "free", labeller = as_labeller(model_names))+
  theme_bw()+
  theme( axis.ticks.x=element_blank())
  

pdf("./output/model_assessment/intercept_92sp_parameter_distributions_ang_only.pdf", width = 10, height = 8)
par(mfrow=c(2,1))
best_only
all_models
dev.off()

jpeg("./output/model_assessment/intercept_92sp_parameter_distributions_ang_only.jpg", width = 6, height = 6, units = "in", res = 500)
ggarrange(best_only, all_models, nrow = 2)
# best_only
# all_models
dev.off()

