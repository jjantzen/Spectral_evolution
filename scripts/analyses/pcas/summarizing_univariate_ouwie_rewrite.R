#rewriting summarizing ouwie output scripts

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

#load data
aics_ouwie <- readRDS("./analysis/pca_analysis/pc_univariate_model_aics_92sp.rds")
aovs_univariate <- readRDS("./analysis/pca_analysis/pc_univariate_model_aovs_92sp.rds")

#Summarize output from univariate models and ouwie

#get column names with highest and second highest value for columns 3:9

best_ouwie_models <- aics_ouwie %>%
  rowwise() %>%
  mutate(top_model = names(.)[which.max(c_across(BM1:OUMVA))+2], second_model = tail(head(names(cur_data())[order(c_across(BM1:OUMVA), decreasing = T)+2],2),1))#(.[3:9], 1, function(x) names(x)[maxn(2)(x)]))

best_ouwie_models

best_ouwie_models %>% 
  dplyr::select(iteration, PC_axis, top_model) %>% 
  dplyr::filter(PC_axis %in% c(1:3)) %>% 
  dplyr::group_by(PC_axis, top_model) %>% 
  dplyr::summarize(count = n())


#plot the aic weights for each model for each pc axis (repeated iterations)
long_aics <- pivot_longer(aics_ouwie, c(3:9))
long_aics

#write function for labeling facets
pc_names <- as_labeller(c(`1` = "PC 1", `2` = "PC 2", `3` = "PC 3", `4` = "PC 4", `5` = "PC 5",
                        `6` = "PC 6", `7` = "PC 7", `8` = "PC 8", `9` = "PC 9", `10` = "PC 10"))

#plot boxplots for ouwie models

jpeg("./output/PCAs/ouwie/boxplot_facet_by_pc_model_col.jpg", width = 12, height = 8, units = "in", res = 600)
ggplot(long_aics[which(long_aics$PC_axis %in% c(1:6)),], aes(x = name, y = round(value, 2)))+
  #geom_point()+
  geom_boxplot(aes(fill = as.factor(name)))+
  #geom_jitter(color = "black", size = 0.5)+
  labs(y = "AICw", x = "Model")+
  facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  scale_fill_discrete(name = "Model")+
  theme_bw()
dev.off()

jpeg("./output/PCAs/ouwie/boxplot_facet_by_pc_model_col_plus_jitter.jpg", width = 12, height = 8, units = "in", res = 600)
ggplot(long_aics[which(long_aics$PC_axis %in% c(1:6)),], aes(x = name, y = round(value, 2)))+
  #geom_point()+
  geom_boxplot(aes(fill = as.factor(name)))+
  geom_jitter(color = "black", size = 0.5)+
  labs(y = "AICw", x = "Model")+
  facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  scale_fill_discrete(name = "Model")+
  theme_bw()
dev.off()

#plot anova results 

jpeg("./output/PCAs/ouwie/aov_pvalues_pcas.jpg", width = 4, height = 4, units = "in", res = 600)
pdf("./output/PCAs/ouwie/aov_pvalues_pcas.pdf", width = 4, height = 4)
ggplot(aovs_univariate[which(aovs_univariate$PC_axis %in% c(1:6)),], aes(x = as.factor(PC_axis), y = round(p_value, 3)))+
  #geom_point()+
  #geom_boxplot(aes(fill = as.factor(name)))+
  # geom_jitter(color = "black")+
  # labs(y = "p-value", x = "PC Axis")+
  geom_jitter(aes(color = ifelse(p_value < 0.05, "red", "black")), size = 0.75) +
  scale_color_identity()+
  labs(y = "p-value", x = "PC Axis", title = "Non-phylogenetic ANOVA")+
  geom_hline(aes(yintercept = 0.05), colour = "red")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  #scale_fill_discrete(name = "Model")+
  theme_bw()
dev.off()

##do this
#An eigendecomposition of the Hessian can provide an indication of whether the search returned the maximum likelihood estimates. If all the eigenvalues of the Hessian are positive, then the Hessian is positive definite, and all parameter estimates are considered reliable.



model_list <- readRDS("./analysis/pca_analysis/best_intercept_models_for_pcas_92sp.rds")





##############Alternative plotting #####################
#faceted by pc axis, violin plot, coloured by model
jpeg("./output/PCAs/ouwie/violin_facet_by_pc_model_col.jpg", width = 10, height = 8, units = "in", res = 600)
ggplot(long_aics, aes(x = name, y = round(value, 2)))+
  #geom_point()+
  geom_flat_violin(aes(fill = as.factor(name)))+
  labs(y = "AICw", x = "Model")+
  facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  scale_fill_discrete(name = "Model")+
  theme_bw()
dev.off()

jpeg("./output/PCAs/ouwie/boxplot_model_col.jpg", width = 10, height = 8, units = "in", res = 500)
ggplot(long_aics, aes(x = as.factor(PC_axis), y = round(value, 2)))+
  #geom_point()+
  geom_boxplot(aes(fill = as.factor(name)))+
  labs(y = "AICw", x = "PC Axis")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  scale_fill_discrete(name = "Model")+
  theme_bw()
dev.off()

jpeg("./output/PCAs/ouwie/violin_model_col.jpg", width = 10, height = 8, units = "in", res = 600)
ggplot(long_aics, aes(x = as.factor(PC_axis), y = round(value, 2)))+
  #geom_point()+
  geom_flat_violin(aes(fill = as.factor(name)))+
  labs(y = "AICw", x = "PC Axis")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  scale_fill_discrete(name = "Model")+
  theme_bw()
dev.off()

jpeg("./output/PCAs/ouwie/testing_log.jpg", width = 10, height = 8, units = "in", res = 600)
ggplot(long_aics[which(long_aics$PC_axis == 1),], aes(x = as.factor(name), y = round(value, 2)))+
  #geom_point()+
  geom_boxplot(fill = "red")+
  geom_flat_violin(fill = "blue")+
  geom_jitter(colour = "green")+
  labs(y = "AICw", x = "Model")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  #scale_fill_discrete(name = "Model")+
  theme_bw()+
  scale_y_sqrt()
dev.off()

ggplot(long_aics[which(long_aics$PC_axis == 1),], aes(x = as.factor(name), y = round(value, 2)))+
  #geom_point()+
  geom_boxplot(fill = "red")+
  geom_flat_violin(fill = "blue")+
  geom_jitter(colour = "green")+
  labs(y = "AICw", x = "Model")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  #scale_fill_discrete(name = "Model")+
  theme_bw()


