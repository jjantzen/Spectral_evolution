#trait only models - for most general exploratory model predictors
library(ggplot2)
library(tidyr)

#read data - still missing myc data
complete_data <- readRDS("./data/for_analysis/final_data.rds")

#####VIF assessment
colnames(complete_data) <- c("Species", "Herb", "Shrub", "Tree", "Woody", "Subclass", "Superorder", "Order", "Family", "Shade", "Drought", "Coarse_soil", "Fine_soil", "Medium_soil", "min_pH", "max_pH", "leaf_persistence", "Raunk_lf", "Habitat", "Raunk_broad")

#list of predictors: Woody, Order, Shade, Drought, Coarse_soil, Fine_soil, Medium_soil 
#min_pH, max_pH, leaf_persistence, Raunk_broad

#make trait not factors
complete_data$Woody <- as.character(complete_data$Woody)

#get reduced dataset for only relevant traist
complete_data <- complete_data[,c(1:5,8:17,20)]
colnames(complete_data)

#get only complete cases
complete_data <- complete_data[complete.cases(complete_data),]

colnames(complete_data)

#model with predictor as response - only works with numerical variables - is this a flawed approach?

#woody predictor
model_woody <- lm(as.numeric(as.factor(Woody)) ~ Shade + Drought + leaf_persistence + Raunk_broad, data = complete_data)

model_woody2 <- lm(as.numeric(as.factor(Woody)) ~ Shade + Drought + leaf_persistence + Coarse_soil + Fine_soil, data = complete_data)

summary(model_woody) #strong correlation between woodiness and p of Raunk_broad
summary(model_woody2) #strong correlation with intercept (Shade factor), and Drought (none), and coarse soil yes
rs.woody <- summary(model_woody)$r.squared #0.9321943
rs.woody2 <- summary(model_woody2)$r.squared #0.2297181

#leaf persistence
model_lp <- lm(as.numeric(as.factor(leaf_persistence)) ~ Woody + Shade + Drought + Coarse_soil + Fine_soil, data = complete_data)

summary(model_lp) #strong correlation with intercept (nonwoody), some correlation with medium soil
rs.lp <- summary(model_lp)$r.squared #0.107465
# complete_data[which(complete_data$leaf_persistence == "e"),]

#shade
model_shade <- lm(as.numeric(as.factor(Shade)) ~ Woody + leaf_persistence + Drought + Coarse_soil + Fine_soil, data = complete_data)

summary(model_shade) #strong correlation with intercept (nonwoody) and fine soil
rs.shade <- summary(model_shade)$r.squared #0.1358456

#Drought
model_drought <- lm(as.numeric(as.factor(Drought)) ~ Woody + leaf_persistence + Shade + Coarse_soil + Fine_soil, data = complete_data)

summary(model_drought) #strong correlation with intercept (nonwoody)
rs.drought <- summary(model_drought)$r.squared #0.09937482

# #Raunk broad
# model_raunk <- lm(as.numeric(as.factor(Raunk_broad)) ~ Woody + leaf_persistence + Shade + Drought, data = complete_data)
# 
# summary(model_raunk) #strong correlation with intercept (nonwoody) and woody
# rs.raunk <- summary(model_raunk)$r.squared #0.5259306

#Coarse soil
model_coarse <- lm(as.numeric(as.factor(Coarse_soil)) ~ Woody + leaf_persistence + Shade + Drought + Fine_soil, data = complete_data)

summary(model_coarse) #strong correlation with intercept (nonwoody) and woody and medium soil
rs.coarse <- summary(model_coarse)$r.squared #0.1274963

#Fine soil
model_fine <- lm(as.numeric(as.factor(Fine_soil)) ~ Woody + leaf_persistence + Shade + Drought + Coarse_soil, data = complete_data)

summary(model_fine) #strong correlation with intercept (nonwoody) and shade tolerant
rs.fine <- summary(model_fine)$r.squared #0.1708713

# #Medium soil
# model_med <- lm(as.numeric(as.factor(Medium_soil)) ~ Woody + leaf_persistence + Shade + Drought + Coarse_soil + Fine_soil, data = complete_data)
# 
# summary(model_med) #strong correlation with intercept (nonwoody) and woody, shade intolerant, coarse soil
# rs.med <- summary(model_med)$r.squared #0.2973518

#maybe exclude medium soil and Raunk from predictors

as.numeric(as.factor(complete_data$Raunk_broad))

#get VIFs from R squareds

vif_raunk <- 1/(1-rs.raunk) #2.109
vif_lp <- 1/(1-rs.lp) #1.120404
vif_shade <- 1/(1-rs.shade) #1.157201
vif_drought <- 1/(1-rs.drought) #1.11034
vif_woody <- 1/(1-rs.woody) #14.748
vif_woody2 <- 1/(1-rs.woody2) #1.298226 

vif_coarse <- 1/(1-rs.coarse) #1.146127 
vif_fine <- 1/(1-rs.fine) #1.206085 
vif_med <- 1/(1-rs.med) #1.423187

#excluded raunk because of correlation with woody and inflated VIF when included in woody model
#excluded medium soil because of correlation with coarse soil and woody 
#woody still has slightly higher VIF cause of correlation with coarse soil - but not too bad

#####Instead of doing VIF for categorical variables, need to do chi-squared tests


#chi squared tests
chisq.test(complete_data$Woody, complete_data$Shade) #p=0.677
chisq.test(complete_data$Woody, complete_data$Drought) #p=0.03458 **
chisq.test(complete_data$Woody, complete_data$Coarse_soil) #p=0.01563 **
chisq.test(complete_data$Woody, complete_data$Fine_soil) #p=0.5132
chisq.test(complete_data$Woody, complete_data$Medium_soil) #p=0.6477
chisq.test(complete_data$Woody, complete_data$leaf_persistence) #p=0.2746
chisq.test(complete_data$Woody, complete_data$Raunk_broad) #p=<0.0001 ***

chisq.test(complete_data$Coarse_soil, complete_data$Shade) #p=0.7031
chisq.test(complete_data$Coarse_soil, complete_data$Drought) #p=0.5177
chisq.test(complete_data$Coarse_soil, complete_data$Fine_soil) #p=0.6159
chisq.test(complete_data$Coarse_soil, complete_data$Medium_soil) #p=0.003357 ***
chisq.test(complete_data$Coarse_soil, complete_data$leaf_persistence) #p=1
chisq.test(complete_data$Coarse_soil, complete_data$Raunk_broad) #p=0.01462 **
chisq.test(complete_data$Coarse_soil, complete_data$Woody) #p=0.01563 **

chisq.test(complete_data$Shade, complete_data$Coarse_soil) #p=0.7031
chisq.test(complete_data$Shade, complete_data$Drought) #p=0.3589
chisq.test(complete_data$Shade, complete_data$Fine_soil) #p=0.009578 ***
chisq.test(complete_data$Shade, complete_data$Medium_soil) #p=0.6975
chisq.test(complete_data$Shade, complete_data$leaf_persistence) #p=0.2646
chisq.test(complete_data$Shade, complete_data$Raunk_broad) #p=0.1181
chisq.test(complete_data$Shade, complete_data$Woody) #p=0.677

chisq.test(complete_data$Drought, complete_data$Coarse_soil) #p=0.5177
chisq.test(complete_data$Drought, complete_data$Shade) #p=0.3589
chisq.test(complete_data$Drought, complete_data$Fine_soil) #p=0.9633
chisq.test(complete_data$Drought, complete_data$Medium_soil) #p=0.4937
chisq.test(complete_data$Drought, complete_data$leaf_persistence) #p=0.9119
chisq.test(complete_data$Drought, complete_data$Raunk_broad) #p=0.09901 
chisq.test(complete_data$Drought, complete_data$Woody) #p=0.03258 **

chisq.test(complete_data$Fine_soil, complete_data$Coarse_soil) #p=0.6159
chisq.test(complete_data$Fine_soil, complete_data$Shade) #p=0.009578 ***
chisq.test(complete_data$Fine_soil, complete_data$Drought) #p=0.9633
chisq.test(complete_data$Fine_soil, complete_data$Medium_soil) #p=0.7849
chisq.test(complete_data$Fine_soil, complete_data$leaf_persistence) #p=0.1815
chisq.test(complete_data$Fine_soil, complete_data$Raunk_broad) #p=0.6242 
chisq.test(complete_data$Fine_soil, complete_data$Woody) #p=0.5132

chisq.test(complete_data$Medium_soil, complete_data$Coarse_soil) #p=0.003357 ***
chisq.test(complete_data$Medium_soil, complete_data$Shade) #p=0.6975
chisq.test(complete_data$Medium_soil, complete_data$Drought) #p=0.4937
chisq.test(complete_data$Medium_soil, complete_data$Fine_soil) #p=0.7849
chisq.test(complete_data$Medium_soil, complete_data$leaf_persistence) #p=0.811
chisq.test(complete_data$Medium_soil, complete_data$Raunk_broad) #p=0.007805 *** 
chisq.test(complete_data$Medium_soil, complete_data$Woody) #p=0.6477

chisq.test(complete_data$leaf_persistence, complete_data$Coarse_soil) #p=1
chisq.test(complete_data$leaf_persistence, complete_data$Shade) #p=0.2646
chisq.test(complete_data$leaf_persistence, complete_data$Drought) #p=0.9119
chisq.test(complete_data$leaf_persistence, complete_data$Fine_soil) #p=0.1815
chisq.test(complete_data$leaf_persistence, complete_data$Medium_soil) #p=0.811
chisq.test(complete_data$leaf_persistence, complete_data$Raunk_broad) #p=0.1201 
chisq.test(complete_data$leaf_persistence, complete_data$Woody) #p=0.2746

chisq.test(complete_data$Raunk_broad, complete_data$Coarse_soil) #p=0.01462 **
chisq.test(complete_data$Raunk_broad, complete_data$Shade) #p=0.1181
chisq.test(complete_data$Raunk_broad, complete_data$Drought) #p=0.09901
chisq.test(complete_data$Raunk_broad, complete_data$Fine_soil) #p=0.6242
chisq.test(complete_data$Raunk_broad, complete_data$Medium_soil) #p=0.007805 ***
chisq.test(complete_data$Raunk_broad, complete_data$leaf_persistence) #p=0.1201 
chisq.test(complete_data$Raunk_broad, complete_data$Woody) #p=1.781e ***

#plan to still exclude Medium soil and Raunk - keeping woody and other variables
#can justify reversing and keeping Raunk and excluding woody

#####plot distribution of quantitative variables (histogram) - exclude because not including pH for now

#####plot distribution of qualitative variables (frequency distribution - look for rare combos)
#for interactions, need combinations occuring at least twice

#contingency tables

table(complete_data$Woody, complete_data$Shade) #all combos >2x
table(complete_data$Woody, complete_data$Drought) #only one high and nonwoody, only 2 none and woody
table(complete_data$Woody, complete_data$Coarse_soil) #all combos >2x
table(complete_data$Woody, complete_data$Fine_soil) #all combos >2x
table(complete_data$Woody, complete_data$leaf_persistence) #only one evergreen nonwoody (fern)
table(complete_data$Woody, complete_data$Raunk_broad) #lots of combos with 0-2 occurrences

table(complete_data$Coarse_soil, complete_data$Shade) #all combos >2x
table(complete_data$Coarse_soil, complete_data$Drought) #only one high drought non-coarse, only two no drought tol non-coarse
table(complete_data$Coarse_soil, complete_data$Fine_soil) #all combos >2x
table(complete_data$Coarse_soil, complete_data$leaf_persistence) #all combos >2x
table(complete_data$Coarse_soil, complete_data$Raunk_broad) #lots of combos with 0-2 occurrences

table(complete_data$Shade, complete_data$Drought) #no drought tolerance - 1 intolerant and 2 intermediate shade
table(complete_data$Shade, complete_data$Fine_soil) #all combos >2x
table(complete_data$Shade, complete_data$leaf_persistence) #all combos >2x
table(complete_data$Shade, complete_data$Raunk_broad) #lots of combos with 0-2 occurrences

table(complete_data$Drought, complete_data$Fine_soil) #only 2 with no drought tolerance and non fine soil
table(complete_data$Drought, complete_data$leaf_persistence) #one evergreen with none, two evergreen with high
table(complete_data$Drought, complete_data$Raunk_broad) #lots of combos with 0-2 occurrences

table(complete_data$Fine_soil, complete_data$leaf_persistence) #all combos >2x
table(complete_data$Fine_soil, complete_data$Raunk_broad) #lots of combos with 0-2 occurrences

table(complete_data$leaf_persistence, complete_data$Raunk_broad) #lots of combos with 0-2 occurrences

#because of combos - need to deal with drought tolerance groupings
#switch to high/med vs low/none
#also, exclude fern from analysis?

#####final decision on variables to use (and structure)

#combine highiand medium drought, and low and none drought

complete_data$Drought_bin <- complete_data$Drought
complete_data$Drought_bin[which(complete_data$Drought_bin == "High")] <- "Yes"
complete_data$Drought_bin[which(complete_data$Drought_bin == "Medium")] <- "Yes"

complete_data$Drought_bin[which(complete_data$Drought_bin == "Low")] <- "No"
complete_data$Drought_bin[which(complete_data$Drought_bin == "None")] <- "No"

complete_data$Drought_bin

#exclude fern??? at later stage if needed

#####final decision on variables to use (and structure)

#make dataframe long for plotting
colnames(complete_data)

long_data <- complete_data %>% 
  dplyr::select(c(1,5,8,10,11,15,16,17)) %>% 
  pivot_longer(cols = -c(1))

jpeg("./output/trait_plots/frequency_traits.jpg", width = 15, height = 10, units = "in", res = 600)
ggplot(data = long_data) +
  geom_bar(aes(value, fill = name)) + #fill = value
  facet_wrap(vars(name), scales = "free")
dev.off()

#make variables as factors
final_data <- complete_data[,c(1,5,8,10,11,15,17)]
final_data$Woody <- as.factor(final_data$Woody)
final_data$Shade <- as.factor(final_data$Shade)
final_data$Coarse_soil <- as.factor(final_data$Coarse_soil)
final_data$Fine_soil <- as.factor(final_data$Fine_soil)
final_data$leaf_persistence <- as.factor(final_data$leaf_persistence)
final_data$Drought_bin <- as.factor(final_data$Drought_bin)

str(final_data)

gaps_db <- read.csv("./data/predictors/shade_gaps_data.csv", stringsAsFactors = FALSE)
colnames(gaps_db) <- c("Species", "Location", "Shade", "Drought_bin", "Coarse_soil", "Fine_soil", "refs", "X")

gaps_db <- gaps_db[,c(1:6)]

colnames(final_data)

colnames(complete_data)

first_data$Woody[which(first_data$Species %in% gaps_db$Species)]
topic_db$Foliage.Persistence[which(topic_db$Species %in% gaps_db$Species)]
topic_db$Raunkiaer.Life.Form[which(topic_db$Species %in% gaps_db$Species)]

first_data <- readRDS("./data/tidy/new_combo_data_matched_spectra.rds")
#filled_db <- read.csv("./data/predictors/filling_in_gaps_2.csv", stringsAsFactors = FALSE)
topic_db <- read.csv("./data/TOPIC/TOPIC_JJantzen_Filled.csv", stringsAsFactors = FALSE)

#combine data into gaps
colnames(final_data)

gaps_db_sorted <- gaps_db[order(gaps_db$Species),] #, decreasing = TRUE) #mtcars[order(mpg, cyl),]

gaps_db_sorted$Woody <- first_data$Woody[which(first_data$Species %in% gaps_db_sorted$Species)]
gaps_db_sorted$Species2 <- first_data$Species[which(first_data$Species %in% gaps_db_sorted$Species)]
gaps_db_sorted$leaf_persistence <- topic_db$Foliage.Persistence[which(topic_db$Species %in% gaps_db_sorted$Species)]
gaps_db_sorted$Species3 <- topic_db$Species[which(topic_db$Species %in% gaps_db_sorted$Species)]

gaps_db_sorted[,c(1,7,8,9, 10)]

colnames(gaps_db_sorted)

#reduce dataframe
gaps_db_final <- gaps_db_sorted[,c(1,7,3,5,6,9,4)]

colnames(gaps_db_final)

#"leaf_persistence"

full_100_data <- rbind(final_data, gaps_db_final)

full_100_data <- full_100_data[order(full_100_data$Species),]


saveRDS(full_100_data, "./data/for_analysis/trait_data_100_taxa.rds")


saveRDS(final_data, "./data/for_analysis/reduced_trait_data.rds")

#####is missing data phylogenetically random?
#how do I evaluate this given that I'll be pruning the phylogeny to match only taxa with data?
#our sampling was already biased with respect to woodiness
#and traits were slightly biased towards woody plants?

#need to prune trees and spectral dataset to match smaller trait dataset



######after modeling, assess model stability by removing collinear variables

