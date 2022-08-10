#plotting optima from models with phylogeny (like barplot/boxplots)
#for pc 1:3

library(OUwie)
library(dplyr)
library(mvMORPH)

#read in models
model_list <- readRDS("./analysis/pca_analysis/best_intercept_models_for_pcas_92sp.rds")

#read in data
data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")
tree_myc <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

consensus_tree <- readRDS("./data/for_analysis/myc_consensus_tree.rds")
tips_to_drop <- consensus_tree$tip.label[-which(consensus_tree$tip.label %in% data_spectra$species)]
consensus_tree_pruned <- drop.tip(consensus_tree, tip = tips_to_drop)
consensus_tree_pruned$tip.label <- gsub("_", " ", consensus_tree_pruned$tip.label)

simmap_tree <- readRDS("./analysis/pca_analysis/myc_simmap_trees.rds")

myc_named <- setNames(data_spectra$myc, rownames(data_spectra$spectra))

aics_ouwie <- readRDS("./analysis/pca_analysis/pc_univariate_model_aics_92sp.rds")

best_ouwie_models <- aics_ouwie %>%
  rowwise() %>%
  dplyr::mutate(top_model = names(.)[which.max(c_across(BM1:OUMVA))+2], second_model = tail(head(names(cur_data())[order(c_across(BM1:OUMVA), decreasing = T)+2],2),1))#(.[3:9], 1, function(x) names(x)[maxn(2)(x)]))

best_ouwie_models

#select top 3 pcs for each iteration
best_ouwie_models_top3 <- best_ouwie_models %>% 
  dplyr::select(iteration, PC_axis, top_model, second_model) %>% 
  dplyr::filter(PC_axis %in% c(1:3))


# #confirming that OU1 was only for pc2
best_ouwie_models_top3 %>%
  dplyr::filter(PC_axis == 2) %>%
  #dplyr::select(top_model) %>%
  dplyr::group_by(top_model) %>% 
  dplyr::summarize(count = n())

#list of models to choose from 
top_model_list <- unique(best_ouwie_models_top3$top_model)

#get matching models
for (i in 1:length(model_list)){ #test with 1 first
  #conduct PCA
  pca_best_model <- mvgls.pca(model_list[[i]], plot = FALSE)  
  
  output_parameters <- data.frame(matrix(nrow=3, ncol = 7))
  colnames(output_parameters)[1] <- "iteration"
  class(output_parameters[,1]) <- "numeric"
  colnames(output_parameters)[2] <- "PC_axis"
  class(output_parameters[,2]) <- "numeric"
  colnames(output_parameters)[3] <- "model"
  class(output_parameters[,3]) <- "character"
  colnames(output_parameters)[4] <- "AM_theta"
  class(output_parameters[,4]) <- "numeric"
  colnames(output_parameters)[5] <- "EM_theta"
  class(output_parameters[,5]) <- "numeric"
  colnames(output_parameters)[6] <- "AM_alpha"
  class(output_parameters[,6]) <- "numeric"
  colnames(output_parameters)[7] <- "EM_alpha"
  class(output_parameters[,7]) <- "numeric"

  
  #run ouwie models
  for (j in 1:3){
    #get data
    data <- data.frame(species=rownames(pca_best_model$scores),myc=myc_named, X=as.numeric(pca_best_model$scores[,j]))
    
    #run models for each pc
    model_for_pc <- best_ouwie_models_top3$top_model[which(best_ouwie_models_top3$iteration == i & best_ouwie_models_top3$PC_axis == j)]
    best_model <- OUwie(simmap_tree[[i]], data, model = model_for_pc, simmap.tree = TRUE) #single rate
    
    output_parameters$iteration[j] <- i
    output_parameters$PC_axis[j] <- j
    output_parameters$model[j] <- best_model$model
    output_parameters$AM_theta[j] <- best_model$theta[1]
    output_parameters$EM_theta[j] <- best_model$theta[2]
    output_parameters$AM_alpha[j] <- best_model$solution[1]
    output_parameters$EM_alpha[j] <- best_model$solution[3]
 
    print(j)
    
  }
  print(i)
  if (i == 1){
    output_parameters_all <- output_parameters
    
  } else {
    output_parameters_all <- rbind(output_parameters_all, output_parameters)
    
  }
}


#saveRDS(output_parameters_all, "./analysis/pca_analysis/summary_of_thetas_alphas_top3pcas.rds")

output_parameters_all <- readRDS("./analysis/pca_analysis/summary_of_thetas_alphas_top3pcas.rds")

#get mean optima for each species (first getting rid of ou1 models)

mean_theta_pc1 <- output_parameters_all %>% 
  dplyr::select(PC_axis, AM_theta, EM_theta) %>% 
  dplyr::filter(PC_axis == 1) %>% 
  dplyr::summarise(mean_AM = mean(AM_theta), mean_EM = mean(EM_theta), sd_AM = sd(AM_theta), sd_EM = sd(EM_theta))

mean_theta_pc1

#seems like some outliers
output_parameters_all %>% 
  dplyr::select(PC_axis, AM_theta, EM_theta) %>% 
  dplyr::filter(PC_axis == 3)

str(output_parameters_all)
long_output_parameters <- pivot_longer(output_parameters_all, cols = c(AM_theta, EM_theta, AM_alpha, EM_alpha))

#plot distribution of theta and alpha

jpeg("./output/PCAs/consensus_plots/distribution_theta_alpha_pc_ouwie_models_pc_123_myc.jpg", height = 10, width = 14, units = "in", res = 600)
all_models <- ggplot(long_output_parameters[-which(long_output_parameters$model == "OU1"),], aes(model, value))+ #
  geom_boxplot()+
  xlab("Model")+
  ylab("Parameter value")+
  labs(title = "All models")+
  facet_wrap(~name, scales = "free")+#, labeller = as_labeller(model_names))+
  #geom_hline(yintercept = 100, colour = "red")+
  #geom_hline(yintercept = -100, colour = "red")+
  theme_bw()+
  theme(axis.ticks.x=element_blank())

all_models
dev.off()

#get rid of outliers of theta

no_outliers_ou1 <- long_output_parameters %>% 
  dplyr::filter(!model %in% "OU1") %>% 
  dplyr::filter(between(value, -10, 10))

no_outliers <- long_output_parameters %>% 
  #dplyr::filter(!model %in% "OU1") %>% 
  dplyr::filter(between(value, -10, 10))

only_outliers <- long_output_parameters %>% 
  dplyr::filter(!between(value, -10, 10))

long_output_parameters %>% 
  dplyr::filter(model %in% "OUMVA") %>% 
  dplyr::select(PC_axis) %>% 
  unique()


jpeg("./output/PCAs/consensus_plots/distribution_theta_alpha_pc_ouwie_models_pc_123_myc_no_outliers.jpg", height = 10, width = 14, units = "in", res = 600)
some_models <- ggplot(no_outliers_ou1, aes(model, value))+ #
  geom_boxplot()+
  xlab("Model")+
  ylab("Parameter value")+
  labs(title = "No outliers and no OU1 models")+
  facet_grid(PC_axis ~name, scales = "free")+#, labeller = as_labeller(model_names))+
  #geom_hline(yintercept = 100, colour = "red")+
  #geom_hline(yintercept = -100, colour = "red")+
  theme_bw()+
  theme(axis.ticks.x=element_blank())

some_models
dev.off()

jpeg("./output/PCAs/consensus_plots/distribution_theta_alpha_pc_ouwie_models_pc_123_myc_no_outliers_inc_ou1.jpg", height = 10, width = 14, units = "in", res = 600)
some_models <- ggplot(no_outliers, aes(model, value))+ #
  geom_boxplot()+
  xlab("Model")+
  ylab("Parameter value")+
  labs(title = "No outliers and no OU1 models")+
  facet_grid(PC_axis ~name, scales = "free")+#, labeller = as_labeller(model_names))+
  #geom_hline(yintercept = 100, colour = "red")+
  #geom_hline(yintercept = -100, colour = "red")+
  theme_bw()+
  theme(axis.ticks.x=element_blank())

some_models
dev.off()

#get means of no outliers

mean_thetas <- no_outliers_ou1 %>% 
  dplyr::select(-c(iteration, model)) %>% 
  dplyr::filter(name %in% c("AM_theta", "EM_theta")) %>% 
  dplyr::group_by(PC_axis, name) %>% 
  dplyr::summarise(mean = mean(value), sd = sd(value))

mean_thetas

#get simplified model names
no_outliers$model_simp <- "OUM"
no_outliers$model_simp[which(no_outliers$model %in% "OU1")] <- "OU1"


mean_thetas_inc_ou1 <- no_outliers %>% 
  dplyr::select(-c(iteration, model)) %>% 
  dplyr::filter(name %in% c("AM_theta", "EM_theta")) %>% 
  dplyr::group_by(PC_axis, name, model_simp) %>% 
  dplyr::summarise(mean = mean(value), sd = sd(value))



#read in pc data and trees
combo <- readRDS("./analysis/pca_analysis/scores_for_all_reps_and_pcs.rds")

consensus_tree <- readRDS("./data/for_analysis/myc_consensus_tree.rds")

#read in data
data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")
tree_myc <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

tips_to_drop <- consensus_tree$tip.label[-which(consensus_tree$tip.label %in% data_spectra$species)]

consensus_tree_pruned <- drop.tip(consensus_tree, tip = tips_to_drop)


consensus_tree_pruned$tip.label <- gsub("_", " ", consensus_tree_pruned$tip.label)
myc_named <- setNames(data_spectra$myc, rownames(data_spectra$spectra))

simmap_consensus <- make.simmap(consensus_tree_pruned, myc_named, model="SYM", nsim=100)
summary_simmap <-describe.simmap(simmap_consensus,plot=TRUE,cex=0.7)

simmap_consensus[[1]]$tip.label <- gsub(" ", "_", simmap_consensus[[1]]$tip.label)

#get data as named list
pc1_all<-setNames(combo$V1,combo$species_factor)
names(pc1_all) <- gsub(" ", "_", names(pc1_all))

pc2_all<-setNames(combo$V2,combo$species_factor)
names(pc2_all) <- gsub(" ", "_", names(pc2_all))

pc3_all<-setNames(combo$V3,combo$species_factor)
names(pc3_all) <- gsub(" ", "_", names(pc3_all))

#get colour states

cols_tree<-setNames(c("blue","orange"),unique(myc_named))


ss<-getStates(simmap_consensus[[1]],"tips")
colors<-setNames(c("orange","blue"),c("AM","EM"))

boxcols<-setNames(sapply(ss,function(pc1_all,y) y[which(names(y)==pc1_all)],
                         y=colors),names(ss))

#sort colours to match tree and data
sorted_boxcols <- boxcols[order(match(names(boxcols),simmap_consensus[[1]]$tip.label))]

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

trans_orange <- t_col("orange", 90, name = "trans_orange")
trans_blue <- t_col("blue", 95, name = "trans_blue")
trans_green <- t_col("brown", 90, name = "trans_green")

  
#make plot
jpeg("./output/PCAs/consensus_plots/consensus_myc_pca123_across_reps_boxplot_scores_with_optima.jpg", height = 10, width = 14, units = "in", res = 600)

#par(mfrow=c(1,4))
layout(matrix(c(1,1,1,2,2,3,3,4,4), nrow = 1, ncol = 9, byrow = TRUE))

plot(simmap_consensus[[1]],cols_tree,ftype='i', fsize = 0.75, xlim = c(0,500), mar=c(5.1,1.1,2.1,0.1))
par(mar=c(5.1,0.1,2.1,1.1))


boxplot(pc1_all~factor(names(pc1_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "")
axis(1)
abline(v = 0)
title(xlab="PC 1 scores")
abline(v = mean(no_outliers_ou1$value[which(no_outliers_ou1$PC_axis == 1 & no_outliers_ou1$name == "AM_theta")]), col = "orange")
abline(v = mean(no_outliers_ou1$value[which(no_outliers_ou1$PC_axis == 1 & no_outliers_ou1$name == "EM_theta")]), col = "blue")
rect(xleft=(mean_thetas$mean[1]+mean_thetas$sd[1]), xright=(mean_thetas$mean[1]-mean_thetas$sd[1]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_orange, border = NA)
rect(xleft=(mean_thetas$mean[2]+mean_thetas$sd[2]), xright=(mean_thetas$mean[2]-mean_thetas$sd[2]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_blue, border = NA)


boxplot(pc2_all~factor(names(pc2_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "")
axis(1)
abline(v = 0)
title(xlab="PC 2 scores")
abline(v = mean(no_outliers_ou1$value[which(no_outliers_ou1$PC_axis == 2 & no_outliers_ou1$name == "AM_theta")]), col = "orange")
abline(v = mean(no_outliers_ou1$value[which(no_outliers_ou1$PC_axis == 2 & no_outliers_ou1$name == "EM_theta")]), col = "blue")
rect(xleft=(mean_thetas$mean[3]+mean_thetas$sd[3]), xright=(mean_thetas$mean[3]-mean_thetas$sd[3]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_orange, border = NA)
rect(xleft=(mean_thetas$mean[4]+mean_thetas$sd[4]), xright=(mean_thetas$mean[4]-mean_thetas$sd[4]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_blue, border = NA)


boxplot(pc3_all~factor(names(pc3_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "")
axis(1)
abline(v = 0)
title(xlab="PC 3 scores")
abline(v = mean(no_outliers_ou1$value[which(no_outliers_ou1$PC_axis == 3 & no_outliers_ou1$name == "AM_theta")]), col = "orange")
abline(v = mean(no_outliers_ou1$value[which(no_outliers_ou1$PC_axis == 3 & no_outliers_ou1$name == "EM_theta")]), col = "blue")
rect(xleft=(mean_thetas$mean[5]+mean_thetas$sd[5]), xright=(mean_thetas$mean[5]-mean_thetas$sd[5]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_orange, border = NA)
rect(xleft=(mean_thetas$mean[6]+mean_thetas$sd[6]), xright=(mean_thetas$mean[6]-mean_thetas$sd[6]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_blue, border = NA)


dev.off()



#plot including ou1 for pc2
#make plot
jpeg("./output/PCAs/consensus_plots/consensus_myc_pca123_across_reps_boxplot_scores_with_optima_inc_ou1.jpg", height = 9, width = 8, units = "in", res = 600)
#pdf("./output/PCAs/consensus_plots/consensus_myc_pca123_across_reps_boxplot_scores_with_optima_inc_ou1.pdf", height = 9, width = 8)

#par(mfrow=c(1,4))
#layout(matrix(c(1,1,1,2,2,3,3,4,4), nrow = 1, ncol = 9, byrow = TRUE))
layout(matrix(c(1,1,1,1,1,2,2,3,3,4,4), nrow = 1, ncol = 11, byrow = TRUE))

plot(simmap_consensus[[1]],cols_tree,ftype='i', fsize = 0.6, xlim = c(0,500), mar=c(5.1,1.1,2.1,0.1))
par(mar=c(5.1,0.1,2.1,1.1))

legend(10, 15,   # Coordinates (x also accepts keywords)
       c("AM", "EM", "Shared"),#legend, # Vector with the name of each group
       #c("orange", "blue", "brown"),#fill,   # Creates boxes in the legend with the specified colors
       col = c("orange", "blue", "brown"),#par("col"), # Color of lines or symbols
       border = "black", # Fill box border color
       lty = c(1,2,4),
       #border = c("orange", "blue", "brown"),
       #lwd = 1,         # Line type and width
       #pch,              # Add pch symbols to legend lines or boxes
       bty = "n",        # Box type (bty = "n" removes the box)
       bg = par("bg"),    # Background color of the legend
       box.lwd = par("lwd"), # Legend box line width
       box.lty = par("lty"), # Legend box line type
       box.col = par("fg"),  # Legend box line color
       cex = 1.5,          # Legend size
       horiz = FALSE,     # Horizontal (TRUE) or vertical (FALSE) legend
       title = "Optima",      # Legend title
       title.adj = 0.1, 
       x.intersp = 0.5, 
       y.intersp = 1,
       xjust = 0
)

legend(10, 25,   # Coordinates (x also accepts keywords)
       c("AM", "EM"),#legend, # Vector with the name of each group
       c("orange", "blue"),#fill,   # Creates boxes in the legend with the specified colors
       #col = c("orange", "blue", "brown"),#par("col"), # Color of lines or symbols
       border = "black", # Fill box border color
       #lty = c(1,2,4),
       #border = c("orange", "blue", "brown"),
       #lwd = 1,         # Line type and width
       #pch,              # Add pch symbols to legend lines or boxes
       bty = "n",        # Box type (bty = "n" removes the box)
       bg = par("bg"),    # Background color of the legend
       box.lwd = par("lwd"), # Legend box line width
       box.lty = par("lty"), # Legend box line type
       box.col = par("fg"),  # Legend box line color
       cex = 1.5,          # Legend size
       horiz = FALSE,     # Horizontal (TRUE) or vertical (FALSE) legend
       title = "Scores",      # Legend title
       title.adj = 0.25, 
       x.intersp = 0.5, 
       y.intersp = 1,
       xjust = 0
)

boxplot(pc1_all~factor(names(pc1_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "", boxlwd = 0.5, medlwd = 1, whisklwd = 1, staplelwd = 1, outlwd = 0.5, outcex  = 0.5)
axis(1, cex.axis = 1.5)
abline(v = 0)
title(xlab="PC 1 scores", cex.lab = 1.5)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 1 & no_outliers$name == "AM_theta" & no_outliers$model_simp %in% "OUM")]), col = "orange", lty = 1)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 1 & no_outliers$name == "EM_theta" & no_outliers$model_simp %in% "OUM")]), col = "blue", lty = 2)
rect(xleft=(mean_thetas_inc_ou1$mean[1]+mean_thetas_inc_ou1$sd[1]), xright=(mean_thetas_inc_ou1$mean[1]-mean_thetas_inc_ou1$sd[1]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_orange, border = NA)
rect(xleft=(mean_thetas_inc_ou1$mean[2]+mean_thetas_inc_ou1$sd[2]), xright=(mean_thetas_inc_ou1$mean[2]-mean_thetas_inc_ou1$sd[2]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_blue, border = NA)


boxplot(pc2_all~factor(names(pc2_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "", boxlwd = 0.5, medlwd = 1, whisklwd = 1, staplelwd = 1, outlwd = 0.5, outcex  = 0.5)
axis(1, cex.axis = 1.5)
abline(v = 0)
title(xlab="PC 2 scores", cex.lab = 1.5)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 2 & no_outliers$name == "AM_theta" & no_outliers$model_simp %in% "OUM")]), col = "orange", lty = 1)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 2 & no_outliers$name == "EM_theta" & no_outliers$model_simp %in% "OUM")]), col = "blue", lty = 2)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 2 & no_outliers$name == "AM_theta" & no_outliers$model_simp %in% "OU1")]), col = "brown", lty = 4)
rect(xleft=(mean_thetas_inc_ou1$mean[4]+mean_thetas_inc_ou1$sd[4]), xright=(mean_thetas_inc_ou1$mean[4]-mean_thetas_inc_ou1$sd[4]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_orange, border = NA)
rect(xleft=(mean_thetas_inc_ou1$mean[6]+mean_thetas_inc_ou1$sd[6]), xright=(mean_thetas_inc_ou1$mean[6]-mean_thetas_inc_ou1$sd[6]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_blue, border = NA)
rect(xleft=(mean_thetas_inc_ou1$mean[3]+mean_thetas_inc_ou1$sd[3]), xright=(mean_thetas_inc_ou1$mean[3]-mean_thetas_inc_ou1$sd[3]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_green, border = NA)

boxplot(pc3_all~factor(names(pc3_all),levels=simmap_consensus[[1]]$tip.label), col = sorted_boxcols, horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(simmap_consensus[[1]])), xlab = "", boxlwd = 0.5, medlwd = 1, whisklwd = 1, staplelwd = 1, outlwd = 0.5, outcex  = 0.5)
axis(1, cex.axis = 1.5)
abline(v = 0)
title(xlab="PC 3 scores", cex.lab = 1.5)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 3 & no_outliers$name == "AM_theta" & no_outliers$model_simp %in% "OUM")]), col = "orange", lty = 1)
abline(v = mean(no_outliers$value[which(no_outliers$PC_axis == 3 & no_outliers$name == "EM_theta" & no_outliers$model_simp %in% "OUM")]), col = "blue", lty = 2)
rect(xleft=(mean_thetas_inc_ou1$mean[7]+mean_thetas_inc_ou1$sd[7]), xright=(mean_thetas_inc_ou1$mean[7]-mean_thetas_inc_ou1$sd[7]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_orange, border = NA)
rect(xleft=(mean_thetas_inc_ou1$mean[8]+mean_thetas_inc_ou1$sd[8]), xright=(mean_thetas_inc_ou1$mean[8]-mean_thetas_inc_ou1$sd[8]), ybottom=1, ytop=Ntip(simmap_consensus[[1]]), col=trans_blue, border = NA)


dev.off()

#get reconstructed optimum spectra

#score*eigenvectors

combo_loadings <- readRDS("./analysis/pca_analysis/vectors_all_phylo_pc1to6.rds")

combo_pc1_loadings <- combo_loadings %>% 
  dplyr::filter(pc_axis == 1) 

#get distribution of scores
output_parameters_all

pc1_thetas <- output_parameters_all %>% 
  dplyr::select(PC_axis, AM_theta, EM_theta, iteration) %>% 
  dplyr::filter(PC_axis == 1)


#example code 
X = iris[,1:4]
mu = colMeans(X) #trait mean (wavelength mean) - but does this differ for phyloPCs?

Xpca = prcomp(X)

nComp = 2
Xhat = Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp]) #score something loading
Xhat = scale(Xhat, center = -mu, scale = FALSE)

Xhat[1,]


#get reconstructed spectra
for (i in 1:100){
  reconstructed_AM <- as.data.frame(matrix(ncol = 3, nrow = 2001), stringsAsFactors = FALSE)
  reconstructed_EM <- as.data.frame(matrix(ncol = 3, nrow = 2001), stringsAsFactors = FALSE)
  colnames(reconstructed_AM) <- c("wavelength", "value", "iteration")
  colnames(reconstructed_EM) <- c("wavelength", "value", "iteration")
  loadings_df <- combo_pc1_loadings[which(combo_pc1_loadings$iteration == i),]
  loadings <- loadings_df[,c(1,3)]
  scores_AM <- pc1_thetas$AM_theta[which(pc1_thetas$iteration == i)]
  scores_EM <- pc1_thetas$EM_theta[which(pc1_thetas$iteration == i)]
  reconstructed_AM$wavelength <- loadings$wavelength
  reconstructed_EM$wavelength <- loadings$wavelength
  reconstructed_AM$value <- loadings$value*scores_AM
  reconstructed_EM$value <- loadings$value*scores_EM
  reconstructed_AM$iteration <- i
  reconstructed_EM$iteration <- i
  
  #get into right format for plotting
  reconstructed_AM_t <- as.data.frame(t(reconstructed_AM))
  colnames(reconstructed_AM_t) <- reconstructed_AM_t[1,]
  reconstructed_AM_t <- reconstructed_AM_t[-1,]
  reconstructed_AM_t <- reconstructed_AM_t[-2,]
  
  reconstructed_EM_t <- as.data.frame(t(reconstructed_EM))
  colnames(reconstructed_EM_t) <- reconstructed_EM_t[1,]
  reconstructed_EM_t <- reconstructed_EM_t[-1,]
  reconstructed_EM_t <- reconstructed_EM_t[-2,]
  
  if (i == 1) {
    reconstructed_AM_all <- reconstructed_AM_t
    reconstructed_EM_all <- reconstructed_EM_t
  } else {
    reconstructed_AM_all <- rbind(reconstructed_AM_all, reconstructed_AM_t)
    reconstructed_EM_all <- rbind(reconstructed_EM_all, reconstructed_EM_t)
  }
}

#get spectra into spectral object for plotting
#convert to spectra object
library(spectrolab)

#transpose

spec_for_plot_AM <- as_spectra(reconstructed_AM_all)
spec_for_plot_EM <- as_spectra(reconstructed_EM_all)


plot(abs(spec_for_plot_EM), col = "red")
plot(spec_for_plot_AM, add = TRUE)

#get means
AM_mean <- mean(spec_for_plot_AM)

EM_mean <- mean(spec_for_plot_EM)

colours_touse <- c("orange", "blue")

colours_trans <- c(make.transparent("orange", 0.15), make.transparent("blue", 0.15))

#plot means
plot(AM_mean, lwd = 0.75, lty = 1, col = "blue", main = "Spectra by mycorrhizal association")
plot(EM_mean, lwd = 0.75, lty = 1, col = "red", add = TRUE)

jpeg("./output/spectra_by_myc.jpg", res = 600, width = 10, height = 10, units = "in")
plot(AM_mean, lwd = 0.75, lty = 1, col = rgb(1,0,0), ylim = c(0,0.6), main = "Spectra by mycorrhizal association (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance")
plot(EM_mean, lwd = 0.75, lty = 1, col = rgb(0,0,1), add = TRUE)
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], total_prob = 0.95, col = rgb(1,0,0,0.1), border = FALSE, add = TRUE) #main = "Spectra by mycorrhizal association (80% quantile)"
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], total_prob = 0.95, col = rgb(0,0,1,0.1), border = FALSE, add = TRUE)
legend(2000, 0.6, legend=c("AM", "EM"), fill=c(rgb(1,0,0,0.25), rgb(0,0,1,0.25)), cex=2)#lty=1,col
dev.off()
