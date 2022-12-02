#combining phylom and regression figs into one

library(phytools)
library(dplyr)
library(ggplot2)
library(OUwie)
library(phylolm)
library(ggforce)
library(ggpubr)
source("./scripts/plot_half_violin.R")


#read data

#aovs
aovs_univariate <- readRDS("./analysis/pca_analysis/pc_univariate_model_aovs_92sp.rds")

adjusted_pvalues <- readRDS("./analysis/pca_analysis/mult_corr_anova_ppc_axes.rds")


#phylolm
output_best <- readRDS("./analysis/pca_analysis/pc_univariate_model_gls_parameters_92sp_best_models.rds")

output_keep <- output_best %>% 
  dplyr::filter(PC_axis %in% c(1:5))

#do multiple correction on reg aovs
adjusted_pvalues <- aovs_univariate %>% 
  dplyr::filter(PC_axis %in% c(1:5)) %>% 
  #dplyr::group_by(iteration) %>% 
  mutate(pval.adj = p.adjust (p_value, method=p.adjust.methods[1]))
#p.adjust(p_value, method = p.adjust.methods[1])


#do multiple correction on phylolm values
adjusted_pvalues_phylolm <- output_keep %>% 
  dplyr::filter(PC_axis %in% c(1:5)) %>% 
  #dplyr::group_by(iteration) %>% 
  mutate(pval.adj = p.adjust (p_value, method=p.adjust.methods[1]))
#p.adjust(p_value, method = p.adjust.methods[1])



# ggplot(aovs_univariate[which(aovs_univariate$PC_axis %in% c(1:5)),], aes(as.factor(PC_axis),round(p_value, 3))) +
#   #geom_violin(position = position_nudge(x = -0.5, y = 0)) +
#   #geom_boxplot(position = position_nudge(x = -0.5, y = 0)) +
#   geom_violin()+
#   geom_sina(alpha=0.5)+
#   # geom_jitter(aes(as.factor(PC_axis),round(p_value, 3), color = ifelse(p_value < 0.05, "red", "black")), size = 0.75) +
#   # scale_color_identity()+
#   # labs(y = "p-value", x = "PC Axis", title = "Non-phylogenetic ANOVA")+
#   # geom_hline(aes(yintercept = 0.05), colour = "red")+
#   theme(legend.position="none")

non_phylo_plot <- ggplot(aovs_univariate[which(aovs_univariate$PC_axis %in% c(1:5)),], aes(x = as.factor(PC_axis), y = round(p_value, 3)))+
  #geom_point()+
  #geom_boxplot(aes(fill = as.factor(name)))+
  # geom_jitter(color = "black")+
  # labs(y = "p-value", x = "PC Axis")+
  geom_jitter(aes(color = ifelse(p_value < 0.05, "red", "black")), size = 0.75) +
  scale_color_identity()+
  labs(y = "p-value", x = "PC Axis")+ #, title = "Non-Phylogenetic ANOVA"
  ylim(0,1)+
  geom_hline(aes(yintercept = 0.05), colour = "red")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  #scale_fill_discrete(name = "Model")+
  theme_bw()


#plot distribution of pvalues for phylo models
phylo_plot <- ggplot(output_keep[which(!is.na(output_keep$p_value)),], aes(x = as.factor(PC_axis), y = round(p_value, 3)))+
  #geom_point()+
  #geom_boxplot(aes(fill = as.factor(name)))+
  #geom_jitter(color = "black")+
  geom_jitter(aes(color = ifelse(p_value < 0.05, "red", "black")), size = 0.75) +
  scale_color_identity()+
  labs(y = "p-value", x = "PC Axis")+ #, title = "Phylogenetic Linear Regression"
  geom_hline(aes(yintercept = 0.05), colour = "red")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  theme_bw()

#plot adjusted pvalues
non_phylo_plot_corr <- ggplot(adjusted_pvalues[which(adjusted_pvalues$PC_axis %in% c(1:5)),], aes(x = as.factor(PC_axis), y = round(pval.adj, 3)))+
  #geom_point()+
  #geom_boxplot(aes(fill = as.factor(name)))+
  # geom_jitter(color = "black")+
  # labs(y = "p-value", x = "PC Axis")+
  geom_jitter(aes(color = ifelse(pval.adj < 0.05, "red", "black")), size = 0.75) +
  scale_color_identity()+
  labs(y = "p-value", x = "PC Axis", title = "Non-Phylogenetic ANOVA")+
  ylim(0,1)+
  geom_hline(aes(yintercept = 0.05), colour = "red")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  #scale_fill_discrete(name = "Model")+
  theme_bw()


#plot distribution of pvalues for phylo models
phylo_plot_corr <- ggplot(adjusted_pvalues_phylolm[which(!is.na(adjusted_pvalues_phylolm$p_value)),], aes(x = as.factor(PC_axis), y = round(pval.adj, 3)))+
  #geom_point()+
  #geom_boxplot(aes(fill = as.factor(name)))+
  #geom_jitter(color = "black")+
  geom_jitter(aes(color = ifelse(pval.adj < 0.05, "red", "black")), size = 0.75) +
  scale_color_identity()+
  labs(y = "p-value", x = "PC Axis", title = "Phylogenetic Linear Regression")+
  geom_hline(aes(yintercept = 0.05), colour = "red")+
  #facet_wrap(~as.factor(PC_axis), ncol = 3, labeller = pc_names)+
  theme_bw()

jpeg("./output/PCAs/ouwie/combined_anova_pvalues_pcas_iterations.jpg", width = 8, height = 4, units = "in", res = 600)
#jpeg("./output/PCAs/ouwie/combined_anova_pvalues_pcas_corrected_pval_iterations.jpg", width = 8, height = 4, units = "in", res = 600)
#pdf("./output/PCAs/ouwie/combined_anova_pvalues.pdf", width = 8, height = 4)
#pdf("./output/PCAs/ouwie/combined_anova_pvalues_pcas_corrected_pval.pdf", width = 8, height = 4)

ggarrange(non_phylo_plot,phylo_plot, labels = c("a", "b"), nrow = 1)

dev.off()
