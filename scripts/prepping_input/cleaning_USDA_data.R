#cleaning USDA data directly downloaded from https://plants.usda.gov/home/downloads

library(stringr)


usda <- read.table('./data/raw/09-02-02PLANTSdata.csv', header = TRUE, sep = ' ', colClasses = 'character')

colnames(usda)

#want columns 3, 5, 6, 10, 13, 22, 25, 26, 27, 28, 30, 31, 37, 38, 39, 36, 43, 44, 48, 49, 50, 53, 54, 56, 57, 58

cols_to_choose <- c(3, 5, 6, 10, 13, 22, 25, 26, 27, 28, 30, 31, 36, 37, 38, 39, 43, 44, 48, 49, 50, 53, 54, 56, 57, 58)

#want rows matching species in species list

#read species list
species <- read.csv("./data/list_of_species_100.csv", header = 1, stringsAsFactors = FALSE)

#select columns
usda_sm <- usda[,cols_to_choose]

colnames(usda_sm)

#select rows

usda_sm$species <- word(usda_sm$Scientific.Name, 1,2, sep=" ")

usda_sm_sp <- usda_sm[which(usda_sm$species %in% species$species),]

nrow(usda_sm_sp)

#write.table(data.out, './data/tidy/USDA_data.out.txt')

species$species[which(species$species %in% usda_sm_sp$species)]

species$species

usda_sm_sp$species[which(usda_sm_sp$Leaf.Retention == "No")]

usda_sm_sp$Growth.Rate
