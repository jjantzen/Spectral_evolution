#Code for doing a PLS-DA (non-phylogenetic comparison) for myc types

##########################
#Florence script#

#_#_#_#_#_#_#_#_#_#_#_#_#
#### load libraries ####
#_#_#_#_#_#_#_#_#_#_#_#_#

library(spectrolab)
library(caret)
library(pls)
library(tidyverse)
library(klaR)
library(hsdar)
library(ape)
library(vegan)
library(agricolae)
library(plotly)
library(reshape2)
library(corrplot)
library(RColorBrewer)
source("./scripts/pick_pls_ncomp_caret.R")   

#_#_#_#_#_#_#_#_#_#
#### load data #####
#_#_#_#_#_#_#_#_#_#

myc_data <- readRDS("./data/for_analysis/myc_data_list_for_analysis.rds")

new_spectra <- readRDS("./data/tidy/new_spectra_matched_trees.rds")

### format df ###
# refls<-data %>% 
#   dplyr::filter(property == 'reflectance') 

t2 <- new_spectra$meta$LS.sample_id #2022 samples

new_spectra #1793 samples

num_sp <- unique(new_spectra$meta$species_names) #2022 samples

#drop non myc species
remove <- c("Chamaedaphne calyculata", "Claytonia perfoliata", "Eriophorum vaginatum", "Kalmia angustifolia", "Populus deltoides", "Rhododendron groenlandicum", 
            "Salix alba", "Salix interior", "Typha angustifolia", "Polystichum munitum", "Tsuga canadensis", "Pinus strobus", "Abies balsamea", "Pinus resinosa", 
            "Thuja occidentalis", "Picea mariana", "Larix laricina", "Picea glauca", "Picea abies", "Pinus banksiana", "Pinus rigida")

refls_sub <- new_spectra %>% 
  subset(!(new_spectra$meta$species_names %in% remove)) #left with 1640 samples

num_sp <- unique(refls_sub$meta$species_names) #79 species

# #drop groups with very low sampling (3< indvd) and qusp - not needed
# few <- refls_sub$meta %>% 
#   group_by(species_names) %>% 
#   summarize(num_samp = n()) %>% 
#   arrange(num_samp) %>% 
#   filter(num_samp < 4)
# 
# fewer_samp <- refls_sub %>% 
#   subset(!(refls_sub$meta$species_names %in% few$species_names)) 
# 
# fewer_samp #1626 samples left

data_frame_refl <- as.data.frame(refls_sub)
colnames(data_frame_refl)

#add myc data
for (i in 1:nrow(data_frame_refl)){
  data_frame_refl$myc[i] <- myc_data$myc[which(myc_data$species == data_frame_refl$species_names[i])]  
}

data_frame_refl$myc

#Isolate wvl+ categories
data <- data_frame_refl %>% 
  dplyr::select(myc, '400':'2400') %>%
  mutate(myc = as.factor(myc))

dim(data) #91 species and 1640 obs

#data_frame_refl[,c(1,2003)] 
dim(data)

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
##### single run to find ncomp ####
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

# Pick the row indices for the training data
training_idxs = caret::createDataPartition(y     = data$myc,
                                           times = 1,
                                           p     = 0.75,
                                           list  = TRUE)[[1]]

# Split the data for training and for testing
train_data <- data[ training_idxs, ]
test_data  <- data[ - training_idxs, ]

r <- table(data$myc) %>% as.data.frame()
write.csv(r, file = "./analysis/plsda/categories_count_only_ang.csv")


# FB:if all trainmod on equal number, would be 14, Prni and Piba removed

# Du: We're using very few bootstrap replicates, but that should be enough to figure
# out the number of components given the variability of the data.

# Comment AS: here I would run a couple of iterations (50 or 100) such that we
# can pick ncomp based on n.s. differences in Kappa between components: i.e. by running
# Tukey tests. Often this procedure gives a smaller no of comps, one that is not
# worse than a higher number of ncomps

ncomp_train_ctrl <- caret::trainControl(method            = "boot",
                                       number            = 20,
                                       repeats           = NA,
                                       p                 = 0.5,
                                       returnData        = FALSE,
                                       returnResamp      = "final",
                                       savePredictions   = "final",
                                       sampling          = "down",
                                       allowParallel     = TRUE)

# Initial fit. This is used to find out the number of components
#Kappa simulations
nsims <- 50
ncomp_train_fit<-list()
for (i in seq(nsims)){
  print(i)
  ncomp_train_fit[[i]]  = caret::train(myc ~ .,
                                       train_data,
                                       method     = "pls",
                                       metric     = "Kappa",
                                       tuneLength = 100,
                                       trControl  = ncomp_train_ctrl)
}

##### pick ncomp ####
## no adjustments ncomps
ncomps_finalmodel <- vector(length = nsims)
for (i in 1:nsims){
  ncomps_finalmodel[i]<-ncomp_train_fit[[i]]$finalModel$ncomp
}
table(ncomps_finalmodel)#tells you the number of times the computer choose number of components
#most: around 26 (17) then 27 (9)

## pick_pls function
# Du: There is probably a smarter way of doing this
# Du: Proceed with caution!
# AS: This is where I would run Tukey tests on Kappa

ncomp<-list()
for (i in seq(nsims)){
  ncomp[[i]] = pick_pls_ncomp_caret(ncomp_train_fit[[i]], SE = TRUE)
}

names(ncomp)<-as.character(seq(nsims))
ncomps_function<-matrix(unlist(ncomp))[,1]
table(ncomps_function) #tells you the number of times the computer choose number of components
#most: 20 is top

## Tuckey
### Kappa statistics
maxcomp<-100 #needs to be equal or higher than tunelenght

kappas <- data.frame(ncomps= 1:maxcomp,matrix(NA, nrow = maxcomp, ncol = length(ncomp_train_fit)))
for (i in 1:length(ncomp_train_fit)){
  kappas[,i+1] <- ncomp_train_fit[[i]]$results$Kappa
}

### Tukey test
kapp <- as.data.frame(as.numeric(t(kappas[,-1])))
kapp <- cbind(kapp, rep(1:maxcomp, each=length(ncomp_train_fit)))
names(kapp) <- c("Kappa", "ncomps")

kapp$ncomps <- as.factor(kapp$ncomps)

modi <- lm (Kappa~ncomps, kapp) #linear model between kappa and the number of components
tuk <- HSD.test (modi,"ncomps",alpha = 0.05)  #test it
#we are interested in the groups
#you can choose the significance level with alpha argument

#data frame out of everything
tuk_dat <- as.data.frame(tuk$groups)
tuk_dat$var <- as.numeric(row.names(tuk_dat))
tuk_dat <- tuk_dat[order(tuk_dat$var,decreasing = F),]
letters <- as.character(tuk_dat$groups)

### Kappa plot

par(bty="l")
boxplot(kapp$Kappa~kapp$ncomps,ylim=c(0,1),
        xlab="Number of components",ylab="Kappa")
text(x=1:maxcomp, y=rep(1,maxcomp),letters, cex = 0.4)
#choose the number of components based on plot :)
#plot is too large to see numbers: reduce plot between comps 10-60

keep<-12:32
kapp_dim<-kapp %>%
  filter(ncomps %in% keep) %>%
  droplevels()

maxi<-length(keep)
modi <- lm (Kappa~ncomps, kapp_dim) #linear model between kappa and the number of components
tuk <- HSD.test (modi,"ncomps",alpha = 0.05)  #test it
#we are interested in the groups
#you can choose the significant level with alpha argument

#data frame out of everything
tuk_dat <- as.data.frame(tuk$groups)
tuk_dat$var <- as.numeric(row.names(tuk_dat))
tuk_dat <- tuk_dat[order(tuk_dat$var,decreasing = F),]
letters <- as.character(tuk_dat$groups)

pdf("./output/plsda_kappa_comparison_finding_ncomp_ang_only.pdf", height = 8, width = 10)
boxplot(kapp_dim$Kappa~kapp_dim$ncomps,ylim=c(0,1),
        xlab="Number of components",ylab="Kappa")
text(x=1:maxi, y=rep(1,maxi),letters, cex = 0.6)
dev.off()

### if considering all models with shared letter equivalent: ncomp <30.
### if consifering last model sharing a : ncomp = 45 (already lower than previous methods)

####JJ: going to proceed with 31 (most reps in first test, and lowest num comp with only a)
####JJ: alt, can run with 24, lowest num comp including a, and 36, highest num comp including a

####JJ: second round, chose 20 because lowest num comp which includes a and is best from dudu's code

#_#_#_#_#_#_#_#_#_#_#_#
#### Iterate PLSDA ####
#_#_#_#_#_#_#_#_#_#_#_#

## things and dummies
ncomp<-20
nsims <-100
mods<-list()
train_fit<-list()
test_probs<-list()
confus<-list()
probis<-list()

## iterate analysis
## for more manageability, training, fitting, testing, confusion and probability matrix could all be computed separately

for (nsim in seq(nsims)){
  print(nsim)
  flush.console()
  set.seed(nsim)
  # Pick the row indices for the training data
  training_idxs = caret::createDataPartition(y     = data$myc,
                                             p     = 0.75,
                                             list  = TRUE)[[1]]
  #separate train and test data
  train_data = data[ training_idxs, ]
  test_data  = data[ - training_idxs, ]
  
  testclass <- as.factor(data$myc[(-training_idxs)])
  
  #train control
  train_ctrl = caret::trainControl(method = "none")
  tune_grid  = data.frame(ncomp = ncomp)
  
  
  #fit model
  train_fit[[nsim]]  = caret::train(myc ~ .,
                                    train_data,
                                    method     = "pls",
                                    tuneGrid   = tune_grid,
                                    trControl  = train_ctrl)
  # Predict the test data
  mods[[nsim]] <- predict(train_fit[[nsim]], newdata = test_data)
  
  # Build confusion matrix iteratively
  confus[[nsim]] <- caret::confusionMatrix(data = mods[[nsim]], reference = testclass)
  
  # Build probability matrix
  test_probs <- predict(train_fit[[nsim]], newdata = test_data, type = "prob")
  probs <- as.data.frame(test_probs)
  names(probs) <- sapply(strsplit(names(probs),split = ".nc"),"[",1)
  probs <- cbind(testclass, probs)
  probis[[nsim]] <- probs 
}


saveRDS(train_fit, "./analysis/plsda/PLSDA_fast_train_T6_ang_only_20comps.rds")
saveRDS(probis, "./analysis/plsda/PLSDA_fast_T6_probis_ang_only_20comps.rds")
saveRDS(confus, "./analysis/plsda/PLSDA_fast_T6_confus_ang_only_20comps.rds")

# #read files
# train_fit <- readRDS("./analysis/plsda/PLSDA_fast_train_T6_31comps.rds")
# probis <- readRDS("./analysis/plsda/PLSDA_fast_T6_probis_31comps.rds")
# confus <- readRDS("./analysis/plsda/PLSDA_fast_T6_confus_31comps.rds")


### Probability plot based on mean per class for a number of iteration
arr <- array(unlist(probis), dim = c(dim(probis[[1]]),nsims))
prob_mean <- apply(arr, 1:2, mean)
prob_mean <- as.data.frame(prob_mean)
prob_mean$V1 <- probis[[1]]$testclass
colnames(prob_mean) <- colnames(probis[[1]])

write.csv(prob_mean,"./analysis/plsda/PLSDA_T6_probmean_ang_only_20comps.csv")

pp <- melt(prob_mean, id="testclass")
pp$position <- ifelse (pp$testclass == pp$variable, 2,1) 

#create colour scale
# No margin
par(mar=c(0,0,1,0))

# Classic palette BuPu, with 4 colors
coul <- brewer.pal(11, "Spectral") 

# Add more colors to this palette :
coul <- colorRampPalette(coul)(2)

# Plot it
pie(rep(1, length(coul)), col = coul , main="") 

jpeg("./output/plsda_id_probability_by_myc_plot_ang_only_20comp.jpg", height = 5, width = 8, res = 500, units = "in")
ggplot(pp, aes(x=testclass, y=value, fill=variable, group=position))+
  geom_bar(stat="identity",position="fill")+
  theme(legend.position = "none")+
  coord_flip()+
  scale_fill_manual(values= coul)+
  theme_bw()+
  theme(legend.title=element_blank(), panel.border = element_blank(),
        panel.grid.major.x = element_line(colour = 'darkgrey'),
        panel.grid.minor.x = element_line(colour = 'darkgrey'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_blank())+
  scale_y_continuous(labels=c('0','25','50','75','100'))+
  labs(x = "Type", y= "Identification probability (%)")+
  theme(legend.position = "none")+
  coord_flip()
dev.off()

#_#_#_#_#_#_#_#_#_#_#_#_#_#
#### Confusion matrix ####
#_#_#_#_#_#_#_#_#_#_#_#_#_#

tabs <- list()
for(i in 1:length(confus)){
  tabs[[i]] <- confus[[i]]$table
}

tabsi <- Reduce('+', tabs)
tab_mean <- as.data.frame.matrix(tabsi/length(confus))

write.csv(tab_mean,"./analysis/plsda/PLSDA_all_confumean_ang_only_20comps.csv")

#with percentage values
sums <- colSums(tab_mean)
tabs_perc <- matrix(NA, length(sums),length(sums))
for (i in 1:length(sums)){
  tabs_perc[,i] <- tab_mean[,i]/sums[i]
}

colnames(tabs_perc) <- colnames(confus[[1]]$table)
rownames(tabs_perc) <- rownames(confus[[1]]$table)

write.csv(tabs_perc,"./analysis/plsda/PLSDA_all_confuperc_ang_only_20comps.csv")

#plot it
tab_mean <- read.csv("./analysis/plsda/PLSDA_all_confumean_ang_only_20comps.csv", row.names = 1)
tabs_perc <- read.csv("./analysis/plsda/PLSDA_all_confuperc_ang_only_20comps.csv", row.names = 1)

#make plot
jpeg("./output/plsda_confusion_matrix_AM_EM_ang_only_20_comps.jpg", height = 5, width = 5, res = 500, units = "in")
par(fig=c(0,1,0,1), new=F)
corrplot::corrplot(tabs_perc, is.corr = T, col = c(coul, coul), tl.col = 1, cl.pos = "n",
                   tl.offset =1, tl.cex = 0.9, tl.srt = 70,
                   addCoef.col ='black', number.font = 1, number.cex = 0.9, addCoefasPercent = T, na.label = " ") 

mtext("Prediction",2, line=2, cex=1.2)
mtext("Reference",3, line = 2.2, cex=1.2)
mtext("Classification accuracy (%)", line = -21, cex=0.9)

dev.off()


#_#_#_#_#_#_#_#_#_#_#_#
#### get loadings ##### 
#_#_#_#_#_#_#_#_#_#_#_#

library(caret)
library(pls)

##extract final model's loadings (in train_fit)
##this extracts all loadings and sums them across all components for each wavelengths 
lls <- list()
for(i in seq(nsim)){
  lls[[i]] <- abs(loadings(train_fit[[i]]$finalModel)[1:dim(loadings(train_fit[[i]]$finalModel))[1],1:ncomp])
  sumis <- lapply(lls,rowSums)
}

mm <- apply(simplify2array(sumis), 1, mean)
ss <- apply(simplify2array(sumis), 1, sd)

mm <- as.data.frame(mm) #mean loading for each band across all simulations
mm <- cbind(mm,ss) #idem but standard deviation ''

mm2<-mm %>% rownames_to_column(var ="wvl") %>% 
  mutate(wvl = as.character(wvl))

#delete quotes
mm2$wvl <-str_replace_all(mm2$wvl, '\`', "")
mm2<-mm2 %>% 
  mutate(wvl =  as.numeric(wvl))

write.csv(mm2, file = "./analysis/plsda/plsda_loadings_ang_only_20comps.csv")

### check plot loadings ####
ggplot(mm2, aes(x= as.numeric(wvl), y= mm))+
  geom_ribbon(aes(ymin=mm-(ss), ymax = mm+(ss)), fill = "grey70", col = "grey70")+
  geom_line(aes(y = mm))+
  scale_x_continuous(name = "Wavelengths (nm)",
                     breaks = c(400,900,1400,1900,2400), 
                     expand = c(0.001, 0.001))+
  scale_y_continuous(name = "loading value")

##extract local maximums 
maxis<-mm2[which(mm2$mm>0.75),]
### almost all wvl higher than 1900
#don't know if noise or.. 
local.max<-c(400, 520, 530, 550, 560, 700, 710, 750, 980, 1000, 1390, 1880, 2300, 2390)

#when two, choose highest
local.max<-c(400, 520, 550, 700, 750, 980, 1000, 1390, 1880, 2300, 2390)

local.max<-c(400, 709, 716, 1391, 1881, 1905, 2397)

#includes:
#green peak of chl reflectance (550) and 400: chl + carotenoids
#anthocyanin (max: 520)
#aromatic esters: 700-760
#water: ~980 (weak feature), ~1450(?))strong feature), 
#dry constituants: 1100-2500: especialy after 1900
#¢ose: 2300


##contruct plot
## the dashed lines are there for the purpose of defining the local maximas
#(e.g: most of 1900+ > 3, so local max defined by those 4+)
#local maxima: even if > threshold, the higest of two consecutive 10 wvl. 
load<-ggplot(mm2, aes(x= as.numeric(wvl), y= mm))+
  geom_ribbon(aes(ymin=mm-(ss), ymax = mm+(ss)), fill = "grey70", col = "grey70")+
  geom_line(aes(y = mm))+
  scale_x_continuous(name = "Wavelengths (nm)",
                     breaks = c(400,900,1400,1900,2400) 
                     # expand = c(0.001, 0.001)
  )+
  scale_y_continuous(name = "sum of loading value (absolute)")+
  #geom_vline(xintercept = maxis$wvl)+
  geom_vline(xintercept = local.max,  color = "grey50")+
  geom_hline(yintercept = c(0.5,1), linetype = 3)+
  theme_classic()+
  theme(plot.margin=unit(c(.5,.5,.5,.5),"cm"))+
  labs(title = "PLS-DA wavelength's absolute loadings")

ggsave("./output/loadings_20comp_ang_only.pdf", width = 12, height =8, dpi = 600 )

##haven't got below here to work JJ


###### comparing loadings to coefficents of variation #####
#cv= sd/mean (dispersion relative des données autour de la moyenne )

### overall cv
# cv.mean<-data %>% 
#   summarise(across(.cols = 2:length(data), mean)) %>% 
#   pivot_longer(cols = everything(), names_to = "wvl", values_to = "wvl.mean")
# cv.sd<-cv.mean<-data %>% 
#   summarise(across(.cols = 2:length(data), sd)) %>% 
#   pivot_longer(cols = everything(), names_to = "wvl", values_to = "wvl.sd")
# 
# cv.df<-left_join(cv.mean, cv.sd, by = "wvl")
# 
# cv.df$cv <- (cv.df$wvl.sd/cv.df$wvl.mean) *100


cv.df<-sapply(data[,2:length(data)], function(x) sd(x) / mean(x) * 100) %>% tibble()
cv.df$wvl<-as.numeric(seq(400,2400))#, by = 10))
names(cv.df)[1] = "cv"

###plot####
ggplot(cv.df, aes(x= as.numeric(wvl), y= cv))+
  geom_line()+
  scale_x_continuous(name = "Wavelengths (nm)",
                     breaks = c(400,900,1400,1900,2400) 
                     # expand = c(0.001, 0.001)
  )+
  scale_y_continuous(name = "coefficient of variation (%)")


##per branch
data2<-left_join(data, unique(metadata[,2:5]), by = "myc") %>% 
  relocate(names(metadata)[2:5], .before = "400") 

cv.branch<-data2 %>% 
  group_by(branch) %>% 
  group_modify(~summarize(.x, across(where(is.numeric), function(x) sd(x) / mean(x) * 100))) %>% 
  pivot_longer(where(is.numeric), names_to = "wvl", values_to = "branch.cv") %>% 
  mutate(wvl = as.numeric(wvl)) %>% 
  left_join(cv.df)

cvs<-ggplot(cv.branch, aes(x= wvl, y= cv, color = "mean"))+
  geom_line(linetype = 2)+
  geom_line(aes(x= wvl, y= branch.cv, color = branch))+
  scale_x_continuous(name = "Wavelengths (nm)",
                     breaks = c(400,900,1400,1900,2400) 
                     # expand = c(0.001, 0.001)
  )+
  scale_y_continuous(name = "coefficient of variation (%)")+
  scale_colour_manual(values = c("#1DACD6", "#F28E2B", "black"))+
  labs(title = "Coefficient of variation of reflectance spectra")
cvs


library(patchwork)

load/cvs

### per sp cv

#most of these represent area between 2 wvl increments (eg. 520-530, 970-980. draw polygons??)

ggplot(mm2, aes(x= as.numeric(wvl), y= mm))+
  geom_rect(aes(xmin=400,xmax=410,ymin=1,ymax=Inf),fill = "lightblue",alpha=0.02)+
  geom_rect(aes(xmin=520,xmax=530,ymin=1,ymax=Inf),fill = "lightblue",alpha=0.02)+
  geom_rect(aes(xmin=550,xmax=560,ymin=1,ymax=Inf),fill = "lightblue",alpha=0.02)+
  geom_rect(aes(xmin=690,xmax=720,ymin=1,ymax=Inf),fill = "lightblue",alpha=0.02)+
  geom_vline(xintercept = 750, color = "lightblue")+
  geom_rect(aes(xmin=970,xmax=980,ymin=1,ymax=Inf),fill = "lightblue",alpha=0.02)+
  geom_rect(aes(xmin=1000,xmax=1010,ymin=1,ymax=Inf),fill = "lightblue",alpha=0.02)+
  geom_vline(xintercept = 1390, color = "lightblue")+
  geom_rect(aes(xmin=1880,xmax=1940,ymin=1,ymax=Inf),fill = "lightblue",alpha=0.02)+
  geom_rect(aes(xmin=1970,xmax=2320,ymin=1,ymax=Inf),fill = "lightblue",alpha=0.02)+
  geom_rect(aes(xmin=2340,xmax=2360,ymin=1,ymax=Inf),fill = "lightblue",alpha=0.02)+
  geom_rect(aes(xmin=2380,xmax=2400,ymin=1,ymax=Inf),fill = "lightblue",alpha=0.02)+
  geom_ribbon(aes(ymin=mm-(ss)*5, ymax = mm+(ss)*5), fill = "grey70", col = "grey70")+
  geom_line(aes(y = mm))+
  scale_x_continuous(name = "Wavelengths (nm)",
                     breaks = c(400,900,1400,1900,2400), 
                     expand = c(0.001, 0.001))+
  scale_y_continuous(name = "abs (loadings)", expand = c(0.001, 0.001))+
  theme_classic()+
  theme(plot.margin=unit(c(.5,.5,.5,.5),"cm"))

text((seq(1, 14, by=1)+0.3), y=1.1,labels = mm$labi,srt=90,xpd = T, pos=2)
#dev.off()

ggsave("./loadings_rangeblue.pdf", width = 12, height =8, dpi = 600 )


###### test ####

#PLSDA_probmean_malpighiales_10comps)
#prob_mal

class_test<-prob_mal$testclass  

class(prob_mal)
prob_mal<-as.data.frame(prob_mal)

pp_test <- melt(prob_mal, id="testclass")
pp_test$position <- ifelse (pp_test$testclass == pp_test$variable, 2,1) 
pp_test <-pp_test %>% 
  mutate(testclass = as.factor(testclass))
levels(pp_test$testclass)

pp_test$testclass <- factor(pp_test$testclass, levels=rev(levels(pp_test$testclass)))  
levels(pp_test$testclass)

ggplot(pp_test, aes(x=testclass, y=value, fill=variable, group=position))+
  geom_bar(stat="identity",position="fill")+
  #scale_fill_manual(values= alpha(coli,1))+
  theme_bw()+
  theme(legend.title=element_blank(), panel.border = element_blank(),
        panel.grid.major.x = element_line(colour = 'darkgrey'),
        panel.grid.minor.x = element_line(colour = 'darkgrey'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_blank())+
  scale_y_continuous(labels=c('0','25','50','75','100'))+
  labs(x = "Sites", y= "Identification probability (%)")+
  coord_flip()



