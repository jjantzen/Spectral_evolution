#retrieving TRY database data for species

library(dplyr)
library(stringr)

#read species list from CABO
plants <- read.csv("./data/metadata/plants.csv", stringsAsFactors = FALSE)

#read species list from TRY
try <- read.delim("../Data/Predictor data/TryAccSpecies_full.txt", sep = "\t", stringsAsFactors = FALSE)

#edit names to match
CABO_names <- unique(plants$scientific_name)
CABO_short <- word(CABO_names, 1,2, sep=" ")

#match names
keep_try <- try[which(try$AccSpeciesName %in% CABO_short),]
str(keep_try)

#make column into list
output <- paste0(keep_try$AccSpeciesID, collapse=", ")

#write file for TRY query
write.table(output, "./data/predictors/try_name_query.txt")
