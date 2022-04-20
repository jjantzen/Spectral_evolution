#Parsing TRY data

library(tidyr)
library(dplyr)

#read trait data
raw <- read.csv("./data/predictors/try_data.csv", stringsAsFactors = FALSE)

#remove data with no trait name
reduced <- raw[-which(is.na(raw$TraitID)),]

#remove extra columns
reduced <- reduced[,-c(1,2,4,5,6,9,10,14,17,19,22:25)]

colnames(reduced)

#keep select traits

keep <- unique(reduced$DataName)[c(3:10,27)]

reduced_soil <- reduced %>% 
  subset(DataName %in% keep)

unique(reduced_soil$DataName)

#look at data for Mycorrhizal type

myc_type <- reduced_soil %>% 
  select(-Comment) %>% 
  filter(DataName == "Mycorrhizal type") %>% 
  distinct(AccSpeciesName)


myc_type <- reduced_soil %>% 
  select(-c(Comment, DatasetID, ObservationID, DataID)) %>% 
  filter(DataName == "Mycorrhizal type") %>% 
  select(-c(TraitName, Reference)) %>% 
  arrange(AccSpeciesName) 

unique(myc_type$OrigValueStr)

#think this is duplicated from other script (above)
unique(reduced_soil$DataName)

coarse <- reduced_soil %>% 
  select(-c(Comment, DatasetID, ObservationID, DataID)) %>% 
  filter(DataName == "Adapted to Coarse Textured Soils") %>% 
  select(-c(TraitName, Reference, Replicates, StdValue, OrigUncertaintyStr, OrigUnitStr)) %>% 
  arrange(AccSpeciesName) 

fine <- reduced_soil %>% 
  select(-c(Comment, DatasetID, ObservationID, DataID)) %>% 
  filter(DataName == "Adapted to Fine Textured Soils") %>% 
  select(-c(TraitName, Reference, Replicates, StdValue, OrigUncertaintyStr, OrigUnitStr)) %>% 
  arrange(AccSpeciesName) 

medium <- reduced_soil %>% 
  select(-c(Comment, DatasetID, ObservationID, DataID)) %>% 
  filter(DataName == "Adapted to Medium Textured Soils") %>% 
  select(-c(TraitName, Reference, Replicates, StdValue, OrigUncertaintyStr, OrigUnitStr)) %>% 
  arrange(AccSpeciesName) 

ph_max <- reduced_soil %>% 
  select(-c(Comment, DatasetID, ObservationID, DataID)) %>% 
  filter(DataName == "pH, Maximum") %>% 
  select(-c(TraitName, Reference, Replicates, StdValue, OrigUncertaintyStr, OrigUnitStr)) %>% 
  arrange(AccSpeciesName) 

ph_min <- reduced_soil %>% 
  select(-c(Comment, DatasetID, ObservationID, DataID)) %>% 
  filter(DataName == "pH, Minimum") %>% 
  select(-c(TraitName, Reference, Replicates, StdValue, OrigUncertaintyStr, OrigUnitStr)) %>% 
  arrange(AccSpeciesName) 

#don't include these one
soil_type <- reduced_soil %>% 
  select(-c(Comment, DatasetID, ObservationID, DataID)) %>% 
  filter(DataName == "Species adapted to soil type") %>% 
  select(-c(TraitName, Reference, Replicates, StdValue, OrigUncertaintyStr, OrigUnitStr)) %>% 
  arrange(AccSpeciesName) 

microbes <- reduced_soil %>% 
  select(-c(Comment, DatasetID, ObservationID, DataID)) %>% 
  filter(DataName == "Microbial interactions and mycorrhizal status") %>% 
  select(-c(TraitName, Reference, Replicates, StdValue, OrigUncertaintyStr, OrigUnitStr)) %>% 
  arrange(AccSpeciesName) 

#subset to only keeping the above
keep2 <- unique(reduced_soil$DataName)[c(3:7)]

reduced_soil2 <- reduced_soil %>% 
  subset(DataName %in% keep2)

#get refs for soil traits
refs <- reduced_soil2 %>% 
  select(-c(ObservationID, DataID, OrigUnitStr, OrigUncertaintyStr, Replicates, StdValue,Comment))

#all traits come from USDA database
unique(refs$Reference)

#############
#spread data by DataName column - aggregates multiple records by species 

colnames(reduced_soil2)

names_data <- reduced_soil2%>%
  dplyr::select(-c(Comment, TraitName,DatasetID, ObservationID, DataID, OrigUnitStr, OrigUncertaintyStr, Replicates, StdValue, Reference)) %>% 
  filter(AccSpeciesName != "Aronia melanocarpa") %>% 
  group_by(DataName) %>%
  mutate(row = row_number()) #%>% 
  #dplyr::select(-c(DatasetID, ObservationID, DataID, OrigUnitStr, OrigUncertaintyStr, Replicates, StdValue))

spread_data <- names_data %>% 
  group_by(AccSpeciesName) %>% 
  tidyr::pivot_wider(names_from = DataName, values_from = OrigValueStr, values_fill = NA) %>% 
  dplyr::select(-row) #%>% 
  #summarise_at(vars('Adapted to Coarse Textured Soils':'pH, Maximum'),list(~ coalesce(!!! .)))

class(names_data$OrigValueStr)

unique(names_data$AccSpeciesName)
colnames(spread_data)


#save spread data
saveRDS(spread_data, "./data/predictors/spread_try_data_soil.rds")

spread_data <- readRDS("./data/predictors/spread_try_data_soil.rds")

###############
#get list of species to keep 
pruned_ang_con_tree <- readRDS("./data/ang_conifer/rds/pruned_ang_con_tree.rds")
species_names <- pruned_ang_con_tree$tip.label

#species missing from data
species_names[-which(species_names %in% spread_data$AccSpeciesName)]

#prune spread data to match species in tree
spread_data_sp <- spread_data[which(spread_data$AccSpeciesName %in% species_names),]

unique(spread_data_sp$AccSpeciesName)

#save latest data for matching with tree and spectra
saveRDS(spread_data_sp, "./data/predictors/spread_try_data_soil_sp.rds")
