#assess distribution of trait values

library("MASS")
library("rcompanion")
library("lsr")
library("vcd")
library("DescTools")

#read data - still missing myc data
combo_data <- readRDS("./data/for_analysis/final_data.rds")

#contingency tables outline - expected true relationship
table(combo_data$Woody, combo_data$Tree)

#chi squared test
chisq.test(combo_data$Woody, combo_data$Tree) #significant

######assessing traits for myc hypothesis testing

#contingency tables
table(combo_data$Woody, combo_data$Myc_NM)
table(combo_data$Woody, combo_data$Subclass)
table(combo_data$Woody, combo_data$Myc_AM) #only one nonwoody nonAM
table(combo_data$Woody, combo_data$Myc_EM) #no nonwoody EM
table(combo_data$Woody, combo_data$Myc_ErM) #very few ErM, both woody
table(combo_data$Woody, combo_data$leaf_persistence) #very evergreen non woody (4)
table(combo_data$Myc_AM, combo_data$leaf_persistence)
table(combo_data$Myc_EM, combo_data$leaf_persistence)
table(combo_data$Myc_NM, combo_data$leaf_persistence) #few evergreen NM (3)
table(combo_data$Tree, combo_data$leaf_persistence)
table(combo_data$Herb, combo_data$leaf_persistence) #few herbaceous evergreen (4)
table(combo_data$Shrub, combo_data$leaf_persistence) #few shrubby evergreen (3)
table(combo_data$Tree, combo_data$Myc_AM) #only 3 nontree nonAM
table(combo_data$Tree, combo_data$Myc_EM) #only 2 nontree EM
table(combo_data$Herb, combo_data$Myc_AM) #only 1 herbaceous nonAM
table(combo_data$Herb, combo_data$Myc_EM) #no herbaceous EM
table(combo_data$Shrub, combo_data$Myc_AM) #only 3 shrubby nonAM
table(combo_data$Shrub, combo_data$Myc_EM) #only 5 shrubby EM

#chi squared tests
chisq.test(combo_data$Woody, combo_data$Myc_NM)
chisq.test(combo_data$Woody, combo_data$Subclass) #significant
chisq.test(combo_data$Woody, combo_data$Myc_AM) #marginally significant
chisq.test(combo_data$Woody, combo_data$Myc_EM) #significant
chisq.test(combo_data$Tree, combo_data$Myc_AM) #significant
chisq.test(combo_data$Tree, combo_data$Myc_EM) #significant
chisq.test(combo_data$Herb, combo_data$Myc_AM) #marginally significant
chisq.test(combo_data$Herb, combo_data$Myc_EM) #significant
chisq.test(combo_data$Shrub, combo_data$Myc_AM) #not significant
chisq.test(combo_data$Shrub, combo_data$Myc_EM) #not significant
chisq.test(combo_data$Woody, combo_data$leaf_persistence) #not significant
chisq.test(combo_data$Myc_AM, combo_data$leaf_persistence) #significant
chisq.test(combo_data$Myc_EM, combo_data$leaf_persistence) #not significant
chisq.test(combo_data$Myc_NM, combo_data$leaf_persistence) #not significant

#missing data tests
myc_categories <- combo_data[,c(1,8,16:19,43)]
myc_categories$Myc_NM[which(!is.na(myc_categories$Myc_NM))] <- "yes"
myc_categories$Myc_NM[which(is.na(myc_categories$Myc_NM))] <- "no"

myc_categories$Myc_AM[which(!is.na(myc_categories$Myc_AM))] <- "yes"
myc_categories$Myc_AM[which(is.na(myc_categories$Myc_AM))] <- "no"

myc_categories$Myc_EM[which(!is.na(myc_categories$Myc_EM))] <- "yes"
myc_categories$Myc_EM[which(is.na(myc_categories$Myc_EM))] <- "no"

table(myc_categories$Myc_NM, myc_categories$Woody)
#table(myc_categories$Myc_AM, myc_categories$Woody)
#table(myc_categories$Myc_EM, myc_categories$Woody)

table(myc_categories$Myc_NM, myc_categories$leaf_persistence) #3 evergreen species missing data for myc
#table(myc_categories$Myc_AM, myc_categories$leaf_persistence)
#table(myc_categories$Myc_EM, myc_categories$leaf_persistence)

chisq.test(myc_categories$Myc_NM, myc_categories$Woody) #not significant (missing/not for woody)
chisq.test(myc_categories$Myc_NM, myc_categories$leaf_persistence) #not significant
chisq.test(myc_categories$Myc_NM, myc_categories$Woody) 

colnames(combo_data)


######other traits for exploratory PC analyses

#contingency tables
table(combo_data$Woody, combo_data$ShadeTolerance_Qual)
table(combo_data$Woody, combo_data$DroughtTolerance_Qual)
table(combo_data$ShadeTolerance_Qual, combo_data$Myc_AM)
table(combo_data$ShadeTolerance_Qual, combo_data$Myc_EM)
table(combo_data$DroughtTolerance_Qual, combo_data$Myc_EM)
table(combo_data$DroughtTolerance_Qual, combo_data$Myc_AM)
table(combo_data$DroughtTolerance_Quant, combo_data$Woody)
table(combo_data$ShadeTolerance_Quant, combo_data$Woody)
table(combo_data$coarse_soil, combo_data$Woody)
table(combo_data$fine_soil, combo_data$Woody)
table(combo_data$medium_soil, combo_data$Woody)

#chi squared tests
chisq.test(combo_data$Woody, combo_data$ShadeTolerance_Qual)
chisq.test(combo_data$Woody, combo_data$DroughtTolerance_Qual)
chisq.test(combo_data$ShadeTolerance_Qual, combo_data$Myc_AM)
chisq.test(combo_data$ShadeTolerance_Qual, combo_data$Myc_EM)
chisq.test(combo_data$DroughtTolerance_Qual, combo_data$Myc_EM)
chisq.test(combo_data$DroughtTolerance_Qual, combo_data$Myc_AM)
chisq.test(combo_data$DroughtTolerance_Qual, combo_data$Myc_NM)
chisq.test(combo_data$ShadeTolerance_Qual, combo_data$Myc_NM)
chisq.test(combo_data$coarse_soil, combo_data$Woody) #significant
chisq.test(combo_data$fine_soil, combo_data$Woody)
chisq.test(combo_data$medium_soil, combo_data$Woody) 
chisq.test(combo_data$medium_soil, combo_data$Myc_NM)
chisq.test(combo_data$fine_soil, combo_data$Myc_NM) 
chisq.test(combo_data$coarse_soil, combo_data$Myc_NM)
chisq.test(combo_data$coarse_soil, combo_data$Myc_EM)
chisq.test(combo_data$fine_soil, combo_data$Myc_EM)
chisq.test(combo_data$medium_soil, combo_data$Myc_EM)
chisq.test(combo_data$coarse_soil, combo_data$DroughtTolerance_Qual)
chisq.test(combo_data$fine_soil, combo_data$DroughtTolerance_Qual)
chisq.test(combo_data$medium_soil, combo_data$DroughtTolerance_Qual)
chisq.test(combo_data$coarse_soil, combo_data$Subclass)
chisq.test(combo_data$fine_soil, combo_data$Subclass) #marginal, almost sig
chisq.test(combo_data$medium_soil, combo_data$Subclass)

#testing for missing data and other variables

colnames(combo_data)

shade_categories <- combo_data[,c(1,8,9)]
shade_categories$ShadeTolerance_Quant[which(!is.na(shade_categories$ShadeTolerance_Quant))] <- "yes"
shade_categories$ShadeTolerance_Quant[which(is.na(shade_categories$ShadeTolerance_Quant))] <- "no"

drought_categories <- combo_data[,c(1,8,11)]
drought_categories$DroughtTolerance_Quant[which(!is.na(drought_categories$DroughtTolerance_Quant))] <- "yes"
drought_categories$DroughtTolerance_Quant[which(is.na(drought_categories$DroughtTolerance_Quant))] <- "no"

chisq.test(shade_categories$ShadeTolerance_Quant, shade_categories$Woody) #significant (missing/not for woody)
chisq.test(drought_categories$DroughtTolerance_Quant, drought_categories$Woody) #significant (missing/not for woody)


myc_categories2 <- combo_data[,c(1,10,18)]
myc_categories2$Myc_NM[which(!is.na(myc_categories2$Myc_NM))] <- "yes"
myc_categories2$Myc_NM[which(is.na(myc_categories2$Myc_NM))] <- "no"

chisq.test(myc_categories2$Myc_NM, myc_categories2$ShadeTolerance_Qual) #not significant (missing/not for shade tolerance)


shade_qual_categoriers <- combo_data[,c(1,8,10)]
shade_qual_categoriers$ShadeTolerance_Qual[which(!is.na(shade_qual_categoriers$ShadeTolerance_Qual))] <- "yes"
shade_qual_categoriers$ShadeTolerance_Qual[which(is.na(shade_qual_categoriers$ShadeTolerance_Qual))] <- "no"

chisq.test(shade_qual_categoriers$ShadeTolerance_Qual, shade_qual_categoriers$Woody) #significant


drought_qual_categoriers <- combo_data[,c(1,8,12)]
drought_qual_categoriers$DroughtTolerance_Qual[which(!is.na(drought_qual_categoriers$DroughtTolerance_Qual))] <- "yes"
drought_qual_categoriers$DroughtTolerance_Qual[which(is.na(drought_qual_categoriers$DroughtTolerance_Qual))] <- "no"

chisq.test(drought_qual_categoriers$DroughtTolerance_Qual, drought_qual_categoriers$Woody) #significant relationship between missing data and woodiness

colnames(combo_data)

soil_coarse_cat <- combo_data[,c(1,8,38)]
soil_coarse_cat$coarse_soil[which(!is.na(soil_coarse_cat$coarse_soil))] <- "Yes"
soil_coarse_cat$coarse_soil[which(is.na(soil_coarse_cat$coarse_soil))] <- "No"

chisq.test(soil_coarse_cat$coarse_soil, soil_coarse_cat$Woody) #significant relationship between missing data and woodiness


soil_fine_cat <- combo_data[,c(1,8,39)]
soil_fine_cat$fine_soil[which(!is.na(soil_fine_cat$fine_soil))] <- "Yes"
soil_fine_cat$fine_soil[which(is.na(soil_fine_cat$fine_soil))] <- "No"

chisq.test(soil_fine_cat$fine_soil, soil_fine_cat$Woody) #significant relationship between missing data and woodiness


soil_med_cat <- combo_data[,c(1,8,40)]
soil_med_cat$medium_soil[which(!is.na(soil_med_cat$medium_soil))] <- "Yes"
soil_med_cat$medium_soil[which(is.na(soil_med_cat$medium_soil))] <- "No"

chisq.test(soil_med_cat$medium_soil, soil_med_cat$Woody) #significant relationship between missing data and woodiness