#phylogenetic signal of traits

library(ape)
library(caper)
library(phytools)
source("./external_code/delta_statistic-master/code.R")

#read data  - need to read in the actual final data
#combo_data <- readRDS("./data/for_analysis/final_data.rds") # not final data

combo_data <- readRDS("./data/for_analysis/trait_data_100_taxa.rds")

colnames(combo_data)

combo_data$Species

myc_data <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

#read trees
new_trees <- readRDS("./data/for_analysis/final_trees_matched_spectra.rds")

#tips to drop (for myc trees)
exclude <- new_trees[[1]]$tip.label[-which(new_trees[[1]]$tip.label %in% dimnames(myc_data$spectra)[[1]])]

#drop tips
#myc_trees <- drop.tip(new_trees[[1]], exclude)

myc_trees <- lapply(new_trees,drop.tip,tip=exclude)

#saveRDS(myc_trees, "./data/for_analysis/myc_trees_for_analysis_plotting.rds")
colnames(combo_data)

#calculate phylo signal for binary traits
woody_data <- combo_data[which(!is.na(combo_data$Woody)),c(1,2)]

df_output_woody <- data.frame(matrix(nrow=100, ncol = 4))
colnames(df_output_woody)[1] <- "iteration"
class(df_output_woody[,1]) <- "numeric"
colnames(df_output_woody)[2] <- "Estimated D"
class(df_output_woody[,2]) <- "numeric"
colnames(df_output_woody)[3] <- "Prob of E(D) from no (random) phylo structure"
class(df_output_woody[,3]) <- "numeric"
colnames(df_output_woody)[4] <- "Prob of E(D) from BM phylo structure"
class(df_output_woody[,4]) <- "numeric"

for (i in 1:length(new_trees)){
  out <- phylo.d(data = woody_data, phy = new_trees[[i]], names.col = Species, binvar = Woody, permut = 999) #D = -0.5893118, P of E(D) BM = 0.978979, P of E(D) non BM = 0
  df_output_woody$iteration[i] <- i
  df_output_woody[i,2] <- out$DEstimate
  df_output_woody[i,3] <- out$Pval1
  df_output_woody[i,4] <- out$Pval0
}

df_output_woody
saveRDS(df_output_woody, "./analysis/trait_phylosig/woody_phylosig.rds")

# AM_data <- combo_data[which(!is.na(combo_data$Myc_AM)),c(1,16)]
# phylo.d(data = AM_data, phy = tree1, names.col = Species, binvar = Myc_AM, permut = 999) #D = -0.161692, P of E(D) BM = 0.6816817
# 
# EM_data <- combo_data[which(!is.na(combo_data$Myc_EM)),c(1,17)]
# phylo.d(data = EM_data, phy = tree1, names.col = Species, binvar = Myc_EM, permut = 999) #D = -0.1952848, P of E(D) BM = 0.7547548
# 
# ErM_data <- combo_data[which(!is.na(combo_data$Myc_ErM)),c(1,19)]
# phylo.d(data = ErM_data, phy = tree1, names.col = Species, binvar = Myc_ErM, permut = 999) #D = -2.915203, P of E(D) BM = 1
# 
# NM_data <- combo_data[which(!is.na(combo_data$Myc_NM)),c(1,18)]
# phylo.d(data = NM_data, phy = tree1, names.col = Species, binvar = Myc_NM, permut = 999) #D = 0.4114759, P of E(D) BM = 0.07807808, P of E(D) non BM = 0.002002002

coarse_data <- combo_data[which(!is.na(combo_data$Coarse_soil)),c(1,4)]

colnames(combo_data)

df_output_coarse <- data.frame(matrix(nrow=100, ncol = 4))
colnames(df_output_coarse)[1] <- "iteration"
class(df_output_coarse[,1]) <- "numeric"
colnames(df_output_coarse)[2] <- "Estimated D"
class(df_output_coarse[,2]) <- "numeric"
colnames(df_output_coarse)[3] <- "Prob of E(D) from no (random) phylo structure"
class(df_output_coarse[,3]) <- "numeric"
colnames(df_output_coarse)[4] <- "Prob of E(D) from BM phylo structure"
class(df_output_coarse[,4]) <- "numeric"


for (i in 1:length(new_trees)){
  out <- phylo.d(data = coarse_data, phy = new_trees[[i]], names.col = Species, binvar = Coarse_soil, permut = 999) #D = -0.5893118, P of E(D) BM = 0.978979, P of E(D) non BM = 0
  df_output_coarse$iteration[i] <- i
  df_output_coarse[i,2] <- out$DEstimate
  df_output_coarse[i,3] <- out$Pval1
  df_output_coarse[i,4] <- out$Pval0
}
saveRDS(df_output_coarse, "./analysis/trait_phylosig/coarse_phylosig.rds")
#phylo.d(data = coarse_data, phy = tree, names.col = Species, binvar = Coarse_soil, permut = 999) #D = 0.8142746, P of E(D) BM = 0, P of E(D) non BM = 0.1161161

fine_data <- combo_data[which(!is.na(combo_data$Fine_soil)),c(1,5)]

colnames(combo_data)

df_output_fine <- data.frame(matrix(nrow=100, ncol = 4))
colnames(df_output_fine)[1] <- "iteration"
class(df_output_fine[,1]) <- "numeric"
colnames(df_output_fine)[2] <- "Estimated D"
class(df_output_fine[,2]) <- "numeric"
colnames(df_output_fine)[3] <- "Prob of E(D) from no (random) phylo structure"
class(df_output_fine[,3]) <- "numeric"
colnames(df_output_fine)[4] <- "Prob of E(D) from BM phylo structure"
class(df_output_fine[,4]) <- "numeric"

for (i in 1:length(new_trees)){
  out <- phylo.d(data = fine_data, phy = new_trees[[i]], names.col = Species, binvar = Fine_soil, permut = 999) #D = -0.5893118, P of E(D) BM = 0.978979, P of E(D) non BM = 0
  df_output_fine$iteration[i] <- i
  df_output_fine[i,2] <- out$DEstimate
  df_output_fine[i,3] <- out$Pval1
  df_output_fine[i,4] <- out$Pval0
}
saveRDS(df_output_fine, "./analysis/trait_phylosig/fine_phylosig.rds")
#phylo.d(data = fine_data, phy = tree, names.col = Species, binvar = Fine_soil, permut = 999) #D = 0.9781816, P of E(D) BM = 0, P of E(D) non BM = 0.4364364

drought_data <- combo_data[which(!is.na(combo_data$Drought_bin)),c(1,7)]

colnames(combo_data)

df_output_drought <- data.frame(matrix(nrow=100, ncol = 4))
colnames(df_output_drought)[1] <- "iteration"
class(df_output_drought[,1]) <- "numeric"
colnames(df_output_drought)[2] <- "Estimated D"
class(df_output_drought[,2]) <- "numeric"
colnames(df_output_drought)[3] <- "Prob of E(D) from no (random) phylo structure"
class(df_output_drought[,3]) <- "numeric"
colnames(df_output_drought)[4] <- "Prob of E(D) from BM phylo structure"
class(df_output_drought[,4]) <- "numeric"

for (i in 1:length(new_trees)){
  out <- phylo.d(data = drought_data, phy = new_trees[[i]], names.col = Species, binvar = Drought_bin, permut = 999) #D = -0.5893118, P of E(D) BM = 0.978979, P of E(D) non BM = 0
  df_output_drought$iteration[i] <- i
  df_output_drought[i,2] <- out$DEstimate
  df_output_drought[i,3] <- out$Pval1
  df_output_drought[i,4] <- out$Pval0
}
saveRDS(df_output_drought, "./analysis/trait_phylosig/drought_phylosig.rds")
#phylo.d(data = drought_data, phy = tree, names.col = Species, binvar = Drought_bin, permut = 999) #D = 0.8100299, P of E(D) BM = 0, P of E(D) non BM = 0.08408408

lp_data <- combo_data[which(!is.na(combo_data$leaf_persistence)),c(1,6)]

colnames(combo_data)

df_output_lp <- data.frame(matrix(nrow=100, ncol = 4))
colnames(df_output_lp)[1] <- "iteration"
class(df_output_lp[,1]) <- "numeric"
colnames(df_output_lp)[2] <- "Estimated D"
class(df_output_lp[,2]) <- "numeric"
colnames(df_output_lp)[3] <- "Prob of E(D) from no (random) phylo structure"
class(df_output_lp[,3]) <- "numeric"
colnames(df_output_lp)[4] <- "Prob of E(D) from BM phylo structure"
class(df_output_lp[,4]) <- "numeric"

for (i in 1:length(myc_trees)){
  out <- phylo.d(data = lp_data, phy = new_trees[[i]], names.col = Species, binvar = leaf_persistence, permut = 999) #D = -0.5893118, P of E(D) BM = 0.978979, P of E(D) non BM = 0
  df_output_lp$iteration[i] <- i
  df_output_lp[i,2] <- out$DEstimate
  df_output_lp[i,3] <- out$Pval1
  df_output_lp[i,4] <- out$Pval0
}
saveRDS(df_output_lp, "./analysis/trait_phylosig/leaf_persistence_phylosig.rds")
#phylo.d(data = lp_data, phy = tree, names.col = Species, binvar = leaf_persistence, permut = 999) #D = -0.6889767, P of E(D) BM = 0.96997, P of E(D) non BM = 0

str(myc_data)

#myc_data2 <- as.data.frame(myc_data[c(5,2)])

myc_data2 <- as.data.frame(myc_data$myc)
myc_data2 <- cbind(myc_data2, myc_data$species)
colnames(myc_data2) <- c("myc", "species")

df_output_myc <- data.frame(matrix(nrow=100, ncol = 4))
colnames(df_output_myc)[1] <- "iteration"
class(df_output_myc[,1]) <- "numeric"
colnames(df_output_myc)[2] <- "Estimated D"
class(df_output_myc[,2]) <- "numeric"
colnames(df_output_myc)[3] <- "Prob of E(D) from no (random) phylo structure"
class(df_output_myc[,3]) <- "numeric"
colnames(df_output_myc)[4] <- "Prob of E(D) from BM phylo structure"
class(df_output_myc[,4]) <- "numeric"

for (i in 1:length(myc_trees)){
  out <- phylo.d(data = myc_data2, phy = myc_trees[[i]], names.col = species, binvar = myc, permut = 999) #D = -0.5893118, P of E(D) BM = 0.978979, P of E(D) non BM = 0
  df_output_myc$iteration[i] <- i
  df_output_myc[i,2] <- out$DEstimate
  df_output_myc[i,3] <- out$Pval1
  df_output_myc[i,4] <- out$Pval0
}
saveRDS(df_output_myc, "./analysis/trait_phylosig/myc_AMEM_phylosig_92sp.rds")
#phylo.d(data = myc_data2, phy = tree, names.col = species, binvar = myc, permut = 999) #D = -0.5454606, P of E(D) BM = 0.992993, P of E(D) non BM = 0



#too many states
# gf_data2 <- as.data.frame(myc_data[c(5,3)])
# phylo.d(data = gf_data2, phy = tree, names.col = species, binvar = gf, permut = 999) #D = -0.5454606, P of E(D) BM = 0.992993, P of E(D) non BM = 0

#lp_data2 <- as.data.frame(myc_data[c(5,4)])

lp_data2 <- as.data.frame(myc_data$lp)
lp_data2 <- cbind(lp_data2, myc_data$species)
colnames(lp_data2) <- c("lp", "species")

df_output_lp2 <- data.frame(matrix(nrow=100, ncol = 4))
colnames(df_output_lp2)[1] <- "iteration"
class(df_output_lp2[,1]) <- "numeric"
colnames(df_output_lp2)[2] <- "Estimated D"
class(df_output_lp2[,2]) <- "numeric"
colnames(df_output_lp2)[3] <- "Prob of E(D) from no (random) phylo structure"
class(df_output_lp2[,3]) <- "numeric"
colnames(df_output_lp2)[4] <- "Prob of E(D) from BM phylo structure"
class(df_output_lp2[,4]) <- "numeric"

for (i in 1:length(myc_trees)){
  out <- phylo.d(data = lp_data2, phy = myc_trees[[i]], names.col = species, binvar = lp, permut = 999) #D = -0.5893118, P of E(D) BM = 0.978979, P of E(D) non BM = 0
  df_output_lp2$iteration[i] <- i
  df_output_lp2[i,2] <- out$DEstimate
  df_output_lp2[i,3] <- out$Pval1
  df_output_lp2[i,4] <- out$Pval0
}
saveRDS(df_output_lp2, "./analysis/trait_phylosig/leaf_persistence_myc_phylosig_92sp.rds")

#phylo.d(data = lp_data2, phy = tree, names.col = species, binvar = lp, permut = 999) #D = -0.2812493, P of E(D) BM = 0.7807808, P of E(D) non BM = 0

# med_data <- combo_data[which(!is.na(combo_data$`Adapted to Medium Textured Soils`)),c(1,40)]
# phylo.d(data = med_data, phy = tree1, names.col = Species, binvar = `Adapted to Medium Textured Soils`, permut = 999) #D = 0.756396, P of E(D) BM = 0.1201201, P of E(D) non BM = 0.1941942


# #only works for continuous traits?
# min_pH_data <- setNames(as.numeric(combo_data$`pH, Minimum`), combo_data$Species)
# 
# phylosig(tree1, min_pH_data, method="K", test=TRUE, nsim=999) #K = 0.00414716, pvalue = 0.984985
# phylosig(tree1, min_pH_data, method="lambda", test=TRUE, nsim=999) #lambda = 0.108554, pvalue = 1
# 
# max_pH_data <- setNames(as.numeric(combo_data$`pH, Maximum`), combo_data$Species)
# 
# phylosig(tree1, max_pH_data, method="K", test=TRUE, nsim=999) #K = 0.00265125, pvalue = 0.988989
# phylosig(tree1, max_pH_data, method="lambda", test=TRUE, nsim=999) #lambda = 0.227, pvalue = 0.0754164
# 

#which taxa are missing?

combo_data$Species

#using Borges code

#trait <- c(PASTE_YOUR_TRAIT_VECTOR_HERE) #trait is list in order of that in tree
#trait, tree,lambda0,se,sim,thin,burn

#shade
df_output_shade <- data.frame(matrix(nrow=100, ncol = 3))
colnames(df_output_shade)[1] <- "iteration"
class(df_output_shade[,1]) <- "numeric"
colnames(df_output_shade)[2] <- "deltaA"
class(df_output_shade[,2]) <- "numeric"
colnames(df_output_shade)[3] <- "Pvalue"
class(df_output_shade[,3]) <- "numeric"

for (j in 1:length(new_trees)){
  deltaA_shade <- delta(as.character(combo_data$Shade),new_trees[[j]],0.1,0.0589,10000,10,100)
  
  random_delta <- rep(NA,100)
  for (i in 1:100){
    rtrait <- sample(combo_data$Shade)
    random_delta[i] <- delta(rtrait,new_trees[[j]],0.1,0.0589,10000,10,100)
  }
  p_value <- sum(random_delta>deltaA_shade)/length(random_delta)
  
  df_output_shade$iteration[j] <- j
  df_output_shade[j,2] <- deltaA_shade
  df_output_shade[j,3] <- p_value
}

saveRDS(df_output_shade, "./analysis/trait_phylosig/shade_phylosig.rds")


#growth form
myc_data$gf
as.data.frame(myc_data[c(5,2)])


df_output_gf <- data.frame(matrix(nrow=100, ncol = 3))
colnames(df_output_gf)[1] <- "iteration"
class(df_output_gf[,1]) <- "numeric"
colnames(df_output_gf)[2] <- "deltaA"
class(df_output_gf[,2]) <- "numeric"
colnames(df_output_gf)[3] <- "Pvalue"
class(df_output_gf[,3]) <- "numeric"

for (j in 1:length(myc_trees)){
  deltaA_gf <- delta(as.character(myc_data$gf),myc_trees[[j]],0.1,0.0589,10000,10,100)
  
  random_delta <- rep(NA,100)
  for (i in 1:100){
    rtrait <- sample(myc_data$gf)
    random_delta[i] <- delta(rtrait,myc_trees[[j]],0.1,0.0589,10000,10,100)
  }
  p_value <- sum(random_delta>deltaA_gf)/length(random_delta)
  
  df_output_gf$iteration[j] <- j
  df_output_gf[j,2] <- deltaA_gf
  df_output_gf[j,3] <- p_value
}

plot(myc_trees[[42]]) #run 42 failed

saveRDS(df_output_gf, "./analysis/trait_phylosig/gf_myctrees_phylosig_92sp.rds")

sd(df_output_gf$Pvalue, na.rm = TRUE)
mean(df_output_gf$Pvalue, na.rm = TRUE)

nrow(df_output_gf[which((df_output_gf$Pvalue < 0.05) == TRUE),])
nrow(df_output_lp2[which((df_output_lp2$`Prob of E(D) from BM phylo structure` > 0.95) == TRUE),])

df_output_lp2


#single rep example with plot
deltaA_shade <- delta(as.character(combo_data$Shade),new_trees[[1]],0.1,0.0589,10000,10,100)

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(combo_data$Shade)
  random_delta[i] <- delta(rtrait,new_trees[[1]],0.1,0.0589,10000,10,100)
}
p_value <- sum(random_delta>deltaA_shade)/length(random_delta)
boxplot(random_delta, ylim = c(0,3))
abline(h=deltaA_gf,col="red")
p_value #significant phylo signal in shade 

