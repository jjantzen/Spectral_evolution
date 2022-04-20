#Parsing TRY data

library(tidyr)
library(dplyr)

#read trait data
raw <- read.csv("./data/predictors/try_data.csv", stringsAsFactors = FALSE)

#remove data with no trait name
reduced <- raw[-which(is.na(raw$TraitID)),]

#remove extra columns
reduced <- reduced[,-c(1,2,4,5,6,9,10,14,17,19,22:25)]

# #get info
# colnames(reduced)
# head(reduced)
# unique(reduced$TraitName)
# unique(reduced$DataName)
# head(reduced)[1:13]

##############
#detailed (trait by trait) summary of datavalues into "final set"

#first shade tolerance:
sp_seedling <- reduced %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species tolerance to shade", DataName == "Shade tolerance seedlings (percent of full sunlight)") %>% 
  distinct(AccSpeciesName)

sp_full <- reduced %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species tolerance to shade", DataName == "Tolerance to shade") %>% 
  distinct(AccSpeciesName)

sp_adult <- reduced %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species tolerance to shade", DataName == "Shade tolerance adults (1 high, 9 low shade tolerance)") %>% 
  distinct(AccSpeciesName)

#shade tolerance: exclude adult and seedling datasets because species in the full dataset
#but do the data match?

shade_tolerance <- reduced %>% 
  select(-c(Comment, DatasetID, ObservationID, DataID)) %>% 
  filter(TraitName == "Species tolerance to shade" & AccSpeciesName %in% sp_adult$AccSpeciesName) %>% 
  select(-TraitName) %>% 
  arrange(AccSpeciesName)

#remove alternative metrics (keep only "Tolerance to shade")
reduced_shade <- reduced %>% 
  subset(DataName!="Shade tolerance seedlings (percent of full sunlight)" & DataName!="Shade tolerance adults (1 high, 9 low shade tolerance)")
nrow(reduced)
nrow(reduced_shade)

#next: drought tolerance
drought <- reduced_shade %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species tolerance to drought") %>% #, DataName == "Physiological drought tolerance (Psicrit)") #%>% 
  distinct(DataName)

drought_sp_phys <- reduced_shade %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species tolerance to drought", DataName == "Physiological drought tolerance (Psicrit)") %>% 
  distinct(AccSpeciesName)

drought_sp_tol <- reduced_shade %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species tolerance to drought", DataName == "Tolerance to drought") %>% 
  distinct(AccSpeciesName)

drought_sp_phys$AccSpeciesName %in% drought_sp_tol$AccSpeciesName

drought_tolerance <- reduced_shade %>% 
  select(-c(DatasetID, ObservationID, DataID)) %>% 
  filter(TraitName == "Species tolerance to drought" & AccSpeciesName %in% drought_sp_phys$AccSpeciesName) %>% 
  select(-TraitName) %>% 
  arrange(AccSpeciesName)

drought_tolerance <- reduced_shade %>% 
  select(-c(DatasetID, ObservationID, DataID)) %>% 
  filter(TraitName == "Species tolerance to drought", DataName == "Tolerance to drought") %>% 
  select(-TraitName) %>% 
  arrange(AccSpeciesName) %>% 
  distinct(OrigValueStr)

#reduce dataset - omit "Physiological drought tolerance (Psicrit)"
reduced_drought <- reduced_shade %>% 
  subset(DataName != "Physiological drought tolerance (Psicrit)") %>% 
  arrange(AccSpeciesName)

nrow(reduced)
nrow(reduced_shade)
nrow(reduced_drought)

#next: Mycorrhiza type

myc <- reduced_drought %>% 
  select(-Comment) %>% 
  filter(TraitName == "Mycorrhiza type", DataName == "Mycorrhiza type according to Maherali") %>% 
  distinct(AccSpeciesName)

var <- reduced_drought %>% 
  select(-Comment) %>% 
  filter(TraitName == "Mycorrhiza type") %>% 
  distinct(DataName)

myc <- reduced_drought %>% 
  select(-Comment) %>% 
  filter(TraitName == "Mycorrhiza status and microbial interactions") #%>% 
distinct(DataName)

colnames(myc)

#get rid of  Mycorrhiza: Likelihood the species was in the 'Stable AM Loss' state inferred under our best HRM-model (SI Figure 1, Extended Methods)
#get rid of Mycorrhiza: Likelihood the species was in the 'Labile' state inferred under our best HRM-model (SI Figure 1, Extended Methods)
#get rid of Mycorrhiza status and microbial interactions in TraitName

reduced_myc <- reduced_drought %>% 
  subset(DataName != "Mycorrhiza: Likelihood the species had retained AM interactions under our best HRM-model" & DataName != "Mycorrhiza: Likelihood the species was in the 'Stable AM' state inferred under our best HRM-model (SI Figure 1, Extended Methods)" & DataName != "Mycorrhiza: Likelihood the species had lost AM interactions under our best HRM-model" & DataName != "Mycorrhiza: Likelihood the species was in the 'Stable AM Loss' state inferred under our best HRM-model (SI Figure 1, Extended Methods)" & DataName != "Mycorrhiza: Likelihood the species was in the 'Labile' state inferred under our best HRM-model (SI Figure 1, Extended Methods)" & TraitName != "Mycorrhiza status and microbial interactions")

nrow(reduced)
nrow(reduced_shade)
nrow(reduced_drought)
nrow(reduced_myc)

#next: Plant nitrogen(N) fixation capacity

var <- reduced_myc %>% 
  select(-Comment) %>% 
  filter(TraitName == "Plant nitrogen(N) fixation capacity") %>% 
  distinct(DataName)

nit <- reduced_myc %>% 
  select(-Comment) %>% 
  filter(TraitName == "Plant nitrogen(N) fixation capacity", DataName == "Nitrogen-fixation capacity") %>% 
  distinct(AccSpeciesName)

#keep all nit fix columns for now

#next: Leaf venation type
var <- reduced_myc %>% 
  select(-Comment) %>% 
  filter(TraitName == "Leaf venation type") %>% 
  distinct(DataName)

vasc <- reduced_myc %>% 
  select(-Comment) %>% 
  filter(TraitName == "Leaf venation type") %>% #, DataName == "Nitrogen-fixation capacity") %>% 
  distinct(AccSpeciesName)

#remove leaf venation type because no data
reduced_vasc <- reduced_myc %>% 
  subset(TraitName != "Leaf venation type")

#next: Soil
#Species habitat characterisation / Plant requirement: soil texture
var <- reduced_vasc %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species habitat characterisation / Plant requirement: soil texture") %>% 
  distinct(DataName)

#Species habitat characterisation / Plant requirement: soil pH
var <- reduced_vasc %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species habitat characterisation / Plant requirement: soil pH") %>% 
  distinct(DataName)

#Species habitat characterisation / Plant requirement: soil type
var <- reduced_vasc %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species habitat characterisation / Plant requirement: soil type")# %>% 
distinct(DataName)

#Species habitat characterization: soil hydrology
var <- reduced_vasc %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species habitat characterization: soil hydrology")# %>% 
distinct(DataName)

#Species soil pH requirement
var <- reduced_vasc %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species soil pH requirement") %>% 
  distinct(DataName)

#Species soil moisture requirements
var <- reduced_vasc %>% 
  select(-Comment) %>% 
  filter(TraitName == "Species soil moisture requirements") %>% 
  distinct(AccSpeciesName)

#remove soil type because too subjective and descriptive (for now)
#remove soil hydrology and soil moisture requirements because data limited
reduced_soil <- reduced_vasc %>% 
  subset(TraitName != "Species soil moisture requirements" & TraitName != "Species habitat characterisation / Plant requirement: soil type" & TraitName != "Species habitat characterization: soil hydrology")

nrow(reduced)
nrow(reduced_shade)
nrow(reduced_drought)
nrow(reduced_myc)
nrow(reduced_vasc)
nrow(reduced_soil)

#now, with 7 variables (with sub-variables), spread to bigger dataframe

#############
#spread data by DataName column - aggregates multiple records by species 
#how to get a single consensus value instead of list of values?

colnames(reduced_soil)

spread_data <- reduced_soil %>% 
  dplyr::select(-c(DatasetID, ObservationID, DataID, TraitName, OrigUnitStr, OrigUncertaintyStr, Replicates, StdValue, Comment)) %>% 
  pivot_wider(id_cols = AccSpeciesName, names_from = DataName, values_from = c(OrigValueStr), values_fill = NA)#, OrigUnitStr, ValueKindName))

spread_data_obs <- reduced_soil %>% 
  dplyr::select(-c(DatasetID, DataID, TraitName, OrigUnitStr, OrigUncertaintyStr, Replicates, StdValue, Comment)) %>% #ObservationID, 
  pivot_wider(id_cols = c(AccSpeciesName,ObservationID, Reference), names_from = DataName, values_from = OrigValueStr, values_fill = NA)#, OrigUnitStr, ValueKindName))

spread_data_dup <- reduced_soil %>% 
  dplyr::select(-c(DatasetID, ObservationID, DataID, TraitName, OrigUnitStr, OrigUncertaintyStr, Replicates, StdValue, Comment)) %>% 
  pivot_wider(id_cols = AccSpeciesName, names_from = DataName, values_from = c(OrigValueStr), values_fn = length, values_fill = NA)#, OrigUnitStr, ValueKindName))

# colnames(spread_data_dup)

reduced_soil %>% 
  dplyr::select(AccSpeciesName, DataName, OrigValueStr)

###first, look at references for different variables

#Green 2009 - text
#Niinemets etc - numbers
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Tolerance to shade`))])
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Tolerance to drought`))])

#Werner et al 2018 for mycorrhiza
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Mycorrhiza: Arbuscular mycorrhizal (AM) fungi`))])
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Mycorrhiza: Presence of any non-AM mycorrhizal fungus (binary)`))])

#mark brundrett - look for him
#topic database - canada - contact isabelle aubin

#Iversen 2017 for mycorrhiza type
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Mycorrhiza type according to Maherali`))])

#Asem 2012 for mycorrhiza type
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Original term for mycorrhizal type given by Selivanov`))])

#janet sprent - nit fix
#Werner et al 2014 - only 140 species
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Nitrogen fixing clade`))])
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Does the species engage in symbiotic N2-fixation with bacterial N2-fixing symbionts (nodulation)?`))])
length(unique(spread_data_obs$AccSpeciesName[which(!is.na(spread_data_obs$`Does the species engage in symbiotic N2-fixation with bacterial N2-fixing symbionts (nodulation)?`))]))

#from 18 different pubs - 189 species
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Nitrogen-fixation capacity`))])
length(unique(spread_data_obs$AccSpeciesName[which(!is.na(spread_data_obs$`Nitrogen-fixation capacity`))]))

#Green 2009 and pH
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Adapted to Coarse Textured Soils`))])
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`pH, Minimum`))])

#Fitter et al 1994 for pH
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Species pH requirement ( soil, typical minimum)`))])
unique(spread_data_obs$Reference[which(!is.na(spread_data_obs$`Species pH requirement (water, extreme maximum)`))])

# colnames(spread_data_obs)
# head(spread_data_obs)

######
#Options for summarizing across obs

#create new columns with summaries (eg search grep for word or strings)
colnames(spread_data_dup)

#look at original dataframe
unique(reduced_soil$OrigValueStr[which(reduced_soil$DataName == "Tolerance to shade")])
unique(reduced_soil$OrigValueStr[which(reduced_soil$DataName == "Tolerance to drought")])
unique(reduced_soil$OrigValueStr[which(reduced_soil$DataName == "Mycorrhiza: Arbutoid mycorrhizal (ARB) fungi")])
unique(reduced_soil$OrigValueStr[which(reduced_soil$DataName == "Mycorrhiza: Ericoid Mycorrhiza (ER) fungi")])
#unique(reduced_soil[which(reduced_soil$DataName == "Mycorrhiza: Ericoid Mycorrhiza (ER) fungi"),])

#look at spread data

#overall nit fix
spread_data[,4][[1]]

#species-level summary of shade tolerance
spread_data[which(spread_data$AccSpeciesName == "Abies balsamea"),2][[1]]
spread_data[which(spread_data$AccSpeciesName == "Acer platanoides"),2][[1]]

#species-level summary of nit fix: columns 4, 12, 13
spread_data[which(spread_data$AccSpeciesName == "Abies balsamea"),4][[1]]
spread_data[which(spread_data$AccSpeciesName == "Abies balsamea"),12][[1]]
spread_data[which(spread_data$AccSpeciesName == "Abies balsamea"),13][[1]]

spread_data[which(spread_data$AccSpeciesName == "Acer pensylvanicum"),4][[1]]
spread_data[which(spread_data$AccSpeciesName == "Acer pensylvanicum"),13][[1]]
#reduced_soil[which(reduced_soil$AccSpeciesName == "Acer pensylvanicum" & reduced_soil$DataName == "Nitrogen-fixation capacity"),][]

#species-level of Original term for mycorrhizal type given by Selivanov
spread_data[which(spread_data$AccSpeciesName == "Abies balsamea"),11][[1]]

#get summary of duplicate obs per species
col_dup_num <- spread_data_dup %>% 
  #summarise(counts = sum(`Tolerance to shade` > 0, na.rm = TRUE))
  summarise_all(funs(sum(. > 1 , na.rm = TRUE)))

col_dup_num <- as.data.frame(col_dup_num)

####################
#get summaries
colnames(spread_data)[1:20]
head(spread_data)[1:20]
nrow(spread_data)

#get info 
colnames(reduced_soil)
reduced_soil[which(reduced_soil$DataName == "Tolerance to shade"),]
reduced_soil[c(270,301),]

head(reduced_soil)
unique(reduced_soil$TraitName)
unique(reduced_soil$DataName)
head(reduced)[1:13]

#save output