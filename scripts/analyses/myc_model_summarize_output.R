#read output from myc models and summarize parameters
library(dplyr)
library(ggplot2)

#for file in folder
#files <- list.files("./analysis/myc_models/models/", pattern=NULL, all.files=FALSE, full.names=TRUE)
files <- list.files("./analysis/myc_models/92sp_myc/parameters/", pattern=NULL, all.files=FALSE, full.names=TRUE)


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
saveRDS(output_sorted, "./analysis/myc_models/92sp_myc/summarized_parameters_myc_models_92sp.rds")

output_sorted <- readRDS("./analysis/myc_models/92sp_myc/summarized_parameters_myc_models_92sp.rds")

#output_sorted <- readRDS("./analysis/myc_models/summarized_parameters_myc_models.rds")

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

saveRDS(best_models, "./analysis/myc_models/92sp_myc/best_models_myc_models_92sp.rds")
best_models <- readRDS("./analysis/myc_models/92sp_myc/best_models_myc_models_92sp.rds")


#summarize by number of model
summary <- best_models %>% 
  group_by(best_model) %>% 
  summarize(n = n(), parameter = mean(parameter))

saveRDS(summary, "./analysis/myc_models/92sp_myc/summary_myc_models_92sp.rds")


#get stats for results table
stats <- output_sorted %>% 
  group_by(model) %>% 
  summarize(mean_param = mean(parameter), sd2_param = sd(parameter), mean_GIC =  mean(GIC), sd2_GIC = sd(GIC))#, deltaGIC = (mean(GIC)-minGIC_EB))


#calculate phylogenetic halflife for each tree for the ou models
best_models

trees_ou <- best_models$iteration[which(best_models$best_model == "OU")]

trees_eb <- best_models$iteration[which(best_models$best_model == "EB")]

trees_bm <- best_models$iteration[which(best_models$best_model == "BM")]

trees_lambda <- best_models$iteration[which(best_models$best_model == "lambda")]

#get trees
# trees_all <- readRDS("./data/for_analysis/final_trees_matched_spectra.rds")
trees_all <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

ou_trees <- trees_all[trees_ou]
eb_trees <- trees_all[trees_eb]
bm_trees <- trees_all[trees_bm]
lambda_trees <- trees_all[trees_lambda]

length(lambda_trees)

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

#exclude 26 and 43 trees
mean(best_models_ou$tree_height[-which(best_models_ou$iteration %in% c(26, 43))])
sd(best_models_ou$tree_height[-which(best_models_ou$iteration %in% c(26, 43))])

mean(best_models_ou$halflife[-which(best_models_ou$iteration %in% c(26, 43))])
sd(best_models_ou$halflife[-which(best_models_ou$iteration %in% c(26, 43))])

mean(best_models_ou$parameter[-which(best_models_ou$iteration %in% c(26, 43))])
sd(best_models_ou$parameter[-which(best_models_ou$iteration %in% c(26, 43))])



#get mean max height for these trees
best_models_eb <- best_models[which(best_models$best_model == "EB"),]

for (i in 1:length(eb_trees)){
  best_models_eb$tree_height[i] <- max(nodeHeights(eb_trees[[i]]))
  #$halflife[i] <- log(2)/best_models_ou$parameter[i]
}

mean(best_models_eb$tree_height[-which(best_models_eb$iteration %in% c(4))])
sd(best_models_eb$tree_height[-which(best_models_eb$iteration %in% c(4))])

mean(best_models_eb$parameter[-which(best_models_eb$iteration %in% c(4))])
sd(best_models_eb$parameter[-which(best_models_eb$iteration %in% c(4))])

#BM
best_models_bm <- best_models[which(best_models$best_model == "BM"),]

for (i in 1:length(bm_trees)){
  best_models_bm$tree_height[i] <- max(nodeHeights(bm_trees[[i]]))
  #$halflife[i] <- log(2)/best_models_ou$parameter[i]
}

mean(best_models_bm$tree_height)
sd(best_models_bm$tree_height)

#lambda
best_models_lambda <- best_models[which(best_models$best_model == "lambda"),]

for (i in 1:length(lambda_trees)){
  best_models_lambda$tree_height[i] <- max(nodeHeights(lambda_trees[[i]]))
  #$halflife[i] <- log(2)/best_models_ou$parameter[i]
}

mean(best_models_lambda$tree_height[-which(best_models_lambda$iteration %in% c(48))])
sd(best_models_lambda$tree_height[-which(best_models_lambda$iteration %in% c(48))])

mean(best_models_lambda$parameter[-which(best_models_lambda$iteration %in% c(48))])
sd(best_models_lambda$parameter[-which(best_models_lambda$iteration %in% c(48))])

# max(nodeHeights(ou_trees[[1]]))
# 
# halflife <- log(2)/alpha #in units of branch lengths 

#make plot of distribution of parameters (one plot) and tree lengths(second tree)

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

best_only <- ggplot(best_models[-which(best_models$best_model == "BM"),], aes(best_model, parameter))+
  geom_boxplot()+
  xlab("Model")+
  ylab("Parameter value")+
  labs(title = "Best models")+
  facet_wrap(~best_model, scales = "free", labeller = as_labeller(model_names))+
  theme_bw()+
  theme( axis.ticks.x=element_blank())


pdf("./output/model_assessment/myc_only_92sp_parameter_distributions.pdf", width = 10, height = 8)
par(mfrow=c(2,1))
best_only
all_models
dev.off()

jpeg("./output/model_assessment/myc_only_92sp_parameter_distributions.jpg", width = 6, height = 6, units = "in", res = 500)
ggarrange(best_only, all_models, nrow = 2)
# best_only
# all_models
dev.off()

best_models_myc_models_92sp

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


