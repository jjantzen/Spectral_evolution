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

#remove alternative metrics (keep only "Tolerance to shade")
reduced_shade <- reduced %>% 
  subset(DataName!="Shade tolerance seedlings (percent of full sunlight)" & DataName!="Shade tolerance adults (1 high, 9 low shade tolerance)")

#check dimensions
nrow(reduced)
nrow(reduced_shade)

#next: drought tolerance

#reduce dataset - omit "Physiological drought tolerance (Psicrit)"
reduced_drought <- reduced_shade %>% 
  subset(DataName != "Physiological drought tolerance (Psicrit)") %>% 
  arrange(AccSpeciesName)

#check dimensions
nrow(reduced)
nrow(reduced_shade)
nrow(reduced_drought)

#next: Mycorrhiza type

#get rid of  Mycorrhiza: Likelihood the species was in the 'Stable AM Loss' state inferred under our best HRM-model (SI Figure 1, Extended Methods)
#get rid of Mycorrhiza: Likelihood the species was in the 'Labile' state inferred under our best HRM-model (SI Figure 1, Extended Methods)
#get rid of Mycorrhiza status and microbial interactions in TraitName

reduced_myc <- reduced_drought %>% 
  subset(DataName != "Mycorrhiza: Likelihood the species had retained AM interactions under our best HRM-model" & DataName != "Mycorrhiza: Likelihood the species was in the 'Stable AM' state inferred under our best HRM-model (SI Figure 1, Extended Methods)" & DataName != "Mycorrhiza: Likelihood the species had lost AM interactions under our best HRM-model" & DataName != "Mycorrhiza: Likelihood the species was in the 'Stable AM Loss' state inferred under our best HRM-model (SI Figure 1, Extended Methods)" & DataName != "Mycorrhiza: Likelihood the species was in the 'Labile' state inferred under our best HRM-model (SI Figure 1, Extended Methods)" & TraitName != "Mycorrhiza status and microbial interactions")

#check dimensions
nrow(reduced)
nrow(reduced_shade)
nrow(reduced_drought)
nrow(reduced_myc)

#next: Plant nitrogen(N) fixation capacity

#keep all nit fix columns for now

#next: Leaf venation type

#remove leaf venation type because no data
reduced_vasc <- reduced_myc %>% 
  subset(TraitName != "Leaf venation type")

#next: Soil

#remove soil type because too subjective and descriptive (for now)
#remove soil hydrology and soil moisture requirements because data limited
reduced_soil <- reduced_vasc %>% 
    subset(TraitName != "Species soil moisture requirements" & TraitName != "Species habitat characterisation / Plant requirement: soil type" & TraitName != "Species habitat characterization: soil hydrology")

#check dimensions
nrow(reduced)
nrow(reduced_shade)
nrow(reduced_drought)
nrow(reduced_myc)
nrow(reduced_vasc)
nrow(reduced_soil)

#############
#spread data by DataName column - aggregates multiple records by species 

colnames(reduced_soil)

names_data <- reduced_soil%>%
  group_by(DataName) %>%
  mutate(row = row_number()) %>%
  dplyr::select(-c(DatasetID, ObservationID, DataID, OrigUnitStr, TraitName, OrigUncertaintyStr, Replicates, StdValue))

spread_data <- names_data %>% 
  tidyr::pivot_wider(names_from = DataName, values_from = OrigValueStr, values_fill = NA) %>% 
  dplyr::select(-row)
  
# spread_data <- reduced_soil %>% 
#   dplyr::select(-c(DatasetID, ObservationID, DataID, TraitName, OrigUnitStr, OrigUncertaintyStr, Replicates, StdValue, Comment)) %>% 
#   pivot_wider(id_cols = AccSpeciesName, names_from = DataName, values_from = c(OrigValueStr), values_fn = list, values_fill = NA)#, OrigUnitStr, ValueKindName))

# spread_data_obs <- names_data %>% 
#   dplyr::select(-c(DatasetID, DataID, TraitName, OrigUnitStr, OrigUncertaintyStr, Replicates, StdValue, Comment)) %>% #ObservationID, 
#   pivot_wider(id_cols = c(AccSpeciesName,ObservationID, Reference), names_from = DataName, values_from = OrigValueStr, values_fill = NA)#, OrigUnitStr, ValueKindName))

# spread_data_dup <- reduced_soil %>% 
#   dplyr::select(-c(DatasetID, ObservationID, DataID, TraitName, OrigUnitStr, OrigUncertaintyStr, Replicates, StdValue, Comment)) %>% 
#   pivot_wider(id_cols = AccSpeciesName, names_from = DataName, values_from = c(OrigValueStr), values_fn = length, values_fill = NA)#, OrigUnitStr, ValueKindName))

# colnames(spread_data_dup)

#with spread data, figure out how to select the best traits/reconcile conflicts

#################
#save spread data
saveRDS(spread_data, "./data/predictors/spread_try_data.rds")

spread_data <- readRDS("./data/predictors/spread_try_data.rds")

###############
#get list of species to keep 
pruned_ang_con_tree <- readRDS("./data/ang_conifer/rds/pruned_ang_con_tree.rds")
species_names <- pruned_ang_con_tree$tip.label

#species missing from data
species_names[-which(species_names %in% spread_data$AccSpeciesName)]

#prune spread data to match species in tree
spread_data_sp <- spread_data[which(spread_data$AccSpeciesName %in% species_names),]

unique(spread_data_sp$AccSpeciesName)

#make new dataframe to enter selected data into
fillable_df <- as.data.frame(matrix(ncol = 1, nrow = length(unique(spread_data_sp$AccSpeciesName))))
colnames(fillable_df) <- c("Species")

fillable_df$Species <- unique(spread_data_sp$AccSpeciesName)

###############
#figure out traits

#shade

#get refs
refs <- unique(spread_data_sp$Reference[which(!is.na(spread_data_sp$`Tolerance to shade`))])

refs[3]

#Quantitative shade from Niinemets et al
values1 <- unique(spread_data_sp$`Tolerance to shade`[which(!is.na(spread_data_sp$`Tolerance to shade`) & spread_data_sp$Reference == refs[1])])
shade1 <- spread_data_sp[which(!is.na(spread_data_sp$`Tolerance to shade`) & spread_data_sp$Reference == refs[1]),]

length(unique(shade1$AccSpeciesName))


#Qualitative shade from Green et al
values2 <- unique(spread_data_sp$`Tolerance to shade`[which(!is.na(spread_data_sp$`Tolerance to shade`) & spread_data_sp$Reference == refs[2])])
shade2 <- spread_data_sp[which(!is.na(spread_data_sp$`Tolerance to shade`) & spread_data_sp$Reference == refs[2]),]

length(unique(shade2$AccSpeciesName))

#Quantitative shade from Wirth et al
values3 <- unique(spread_data_sp$`Tolerance to shade`[which(!is.na(spread_data_sp$`Tolerance to shade`) & spread_data_sp$Reference == refs[3])])
shade3 <- spread_data_sp[which(!is.na(spread_data_sp$`Tolerance to shade`) & spread_data_sp$Reference == refs[3]),]

length(unique(shade3$AccSpeciesName))

# #Qualitative shade from unpubl
# unique(spread_data$`Tolerance to shade`[which(!is.na(spread_data$`Tolerance to shade`) & spread_data$Reference == refs[4])])

#Qualitative shade from Fitter et al
values5 <- unique(spread_data_sp$`Tolerance to shade`[which(!is.na(spread_data_sp$`Tolerance to shade`) & spread_data_sp$Reference == refs[5])])
shade5 <- spread_data_sp[which(!is.na(spread_data_sp$`Tolerance to shade`) & spread_data_sp$Reference == refs[5]),]

length(unique(shade5$AccSpeciesName))

#use 1 and 2 as two separate metrics (one quant, one qual)
nrow(shade1)
length(unique(shade1$AccSpeciesName))

colnames(shade1)

#add qual shade variable with comment and reference
fillable_df2 <- fillable_df %>% 
  left_join(shade1[,c(1,4,2,3)], by = c("Species" = "AccSpeciesName"))

colnames(fillable_df2) <- c("Species", "ShadeTolerance_Quant", "ReferenceShade_Quant", "ShadeComment_Quant")

#add quant shade variable with comment and reference
fillable_df2 <- fillable_df2 %>% 
  left_join(shade2[,c(1,4,2)], by = c("Species" = "AccSpeciesName"))

colnames(fillable_df2) <- c("Species", "ShadeTolerance_Quant", "ReferenceShade_Quant", "ShadeComment_Quant", "ShadeTolerance_Qual", "ReferenceShade_Qual")

head(fillable_df2)

#Drought

#get refs
refs <- unique(spread_data_sp$Reference[which(!is.na(spread_data_sp$`Tolerance to drought`))])

refs[2]

#Quantitative drought from Niinemets et al
values1 <- unique(spread_data_sp$`Tolerance to drought`[which(!is.na(spread_data_sp$`Tolerance to drought`) & spread_data_sp$Reference == refs[1])])
drought1 <- spread_data_sp[which(!is.na(spread_data_sp$`Tolerance to drought`) & spread_data_sp$Reference == refs[1]),]

length(unique(drought1$AccSpeciesName))

#Quantitative drought from Green 2009
values2 <- unique(spread_data_sp$`Tolerance to drought`[which(!is.na(spread_data_sp$`Tolerance to drought`) & spread_data_sp$Reference == refs[2])])
drought2 <- spread_data_sp[which(!is.na(spread_data_sp$`Tolerance to drought`) & spread_data_sp$Reference == refs[2]),]

length(unique(drought2$AccSpeciesName))

#add qual drought variable with comment and reference
fillable_df3 <- fillable_df2 %>% 
  left_join(drought1[,c(1,5,2,3)], by = c("Species" = "AccSpeciesName"))

colnames(fillable_df3) <- c("Species", "ShadeTolerance_Quant", "ReferenceShade_Quant", 
                            "ShadeComment_Quant", "ShadeTolerance_Qual", "ReferenceShade_Qual",
                            "DroughtTolerance_Quant", "ReferenceDrought_Quant", "DroughtComment_Quant")

#add quant drought variable with comment and reference
fillable_df3 <- fillable_df3 %>% 
  left_join(drought2[,c(1,5,2)], by = c("Species" = "AccSpeciesName"))

colnames(fillable_df3) <- c("Species", "ShadeTolerance_Quant", "ReferenceShade_Quant", 
                            "ShadeComment_Quant", "ShadeTolerance_Qual", "ReferenceShade_Qual",
                            "DroughtTolerance_Quant", "ReferenceDrought_Quant", "DroughtComment_Quant", 
                            "DroughtTolerance_Qual", "ReferenceDrought_Qual")

head(fillable_df3)

#Nitfix

#get refs
refs <- unique(spread_data_sp$Reference[which(!is.na(spread_data_sp$`Nitrogen-fixation capacity`))])

refs[2]

colnames(spread_data_sp)

#get feel for data columns
values1 <- unique(spread_data_sp$`Nitrogen-fixation capacity`[which(!is.na(spread_data_sp$`Nitrogen-fixation capacity`))]) # & spread_data_sp$Reference == refs[1]
nitfix1 <- spread_data_sp[which(!is.na(spread_data_sp$`Nitrogen-fixation capacity`) ),] #& spread_data_sp$Reference == refs[1]

length(unique(nitfix1$AccSpeciesName))

#get feel for data columns
values2 <- unique(spread_data_sp$`Does the species engage in symbiotic N2-fixation with bacterial N2-fixing symbionts (nodulation)?`[which(!is.na(spread_data_sp$`Does the species engage in symbiotic N2-fixation with bacterial N2-fixing symbionts (nodulation)?`))]) # & spread_data_sp$Reference == refs[1]
nitfix2 <- spread_data_sp[which(!is.na(spread_data_sp$`Does the species engage in symbiotic N2-fixation with bacterial N2-fixing symbionts (nodulation)?`) ),] #& spread_data_sp$Reference == refs[1]

length(unique(nitfix2$AccSpeciesName))

unique(nitfix1$`Nitrogen-fixation capacity`)

#replace synonymous repsonses

nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "no, not an N fixer")] <- "No"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "None")] <- "No"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "no")] <- "No"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "n")] <- "No"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "N")] <- "No"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "0")] <- "No"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "not N2 fixing")] <- "No"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "NO-N-fixer")] <- "No"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "non Fixer")] <- "No"

nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "yes")] <- "Yes"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "2")] <- "Yes"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "Y")] <- "Yes"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "y")] <- "Yes"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "yes, an N fixer")] <- "Yes"
nitfix1$`Nitrogen-fixation capacity`[which(nitfix1$`Nitrogen-fixation capacity` == "N-FIXER")] <- "Yes"

nitfix1$Comment[which(nitfix1$`Nitrogen-fixation capacity` == "Low")]
#0,1,2 (Capacity for nitrogen fixation; 0=none; 1=non-symbiotic rhizosphere N-fixation reported; 2=forms symbiotic nodules with N-fixing bacteria)

unique(nitfix1$`Nitrogen-fixation capacity`)

nitfix_grouped <- nitfix1 %>% 
  dplyr::select(AccSpeciesName, `Nitrogen-fixation capacity`, Reference) %>% 
  group_by(AccSpeciesName, `Nitrogen-fixation capacity`) %>% 
  summarize(Reference = list(unique(Reference)))

nitfix2_grouped <- nitfix2 %>% 
  dplyr::select(AccSpeciesName, `Does the species engage in symbiotic N2-fixation with bacterial N2-fixing symbionts (nodulation)?`, Reference) %>% 
  group_by(AccSpeciesName) %>% 
  summarize(NitFix_v2 = unique(`Does the species engage in symbiotic N2-fixation with bacterial N2-fixing symbionts (nodulation)?`), Reference = Reference)

nitfix_grouped$Reference[2]
nitfix_grouped$Reference[3]

#Need to figure out what to do about conflicts

nitfix1[which(nitfix1$AccSpeciesName == "Acer negundo"),c(1,2,6)]

#all entries from the paper Werner et al 2014 are "yes" regardless of species identity or nit fix capacity (conflicts with their excel supp 2 sheet)
unique(spread_data_sp$`Nitrogen-fixation capacity`[which(spread_data_sp$Reference == "Werner, GDA, WK Cornwell, JI Sprent, J Kattge, ET Kiers (2014) A single evolutionary innovation drives the deep evolution of symbiotic N2-fixation in angiosperms. Nature Communications. DOI: 10.1038/ncomms5087" & !is.na(spread_data_sp$`Nitrogen-fixation capacity`))])

#get nitfix data without bad ref
values3 <- unique(spread_data_sp$`Nitrogen-fixation capacity`[which(!is.na(spread_data_sp$`Nitrogen-fixation capacity`)& spread_data_sp$Reference != "Werner, GDA, WK Cornwell, JI Sprent, J Kattge, ET Kiers (2014) A single evolutionary innovation drives the deep evolution of symbiotic N2-fixation in angiosperms. Nature Communications. DOI: 10.1038/ncomms5087")]) # & spread_data_sp$Reference == refs[1]
nitfix3 <- spread_data_sp[which(!is.na(spread_data_sp$`Nitrogen-fixation capacity`) & spread_data_sp$Reference != "Werner, GDA, WK Cornwell, JI Sprent, J Kattge, ET Kiers (2014) A single evolutionary innovation drives the deep evolution of symbiotic N2-fixation in angiosperms. Nature Communications. DOI: 10.1038/ncomms5087"),] #& spread_data_sp$Reference == refs[1]

#standardize responses again
#replace synonymous repsonses

nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "no, not an N fixer")] <- "No"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "None")] <- "No"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "no")] <- "No"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "n")] <- "No"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "N")] <- "No"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "0")] <- "No"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "not N2 fixing")] <- "No"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "NO-N-fixer")] <- "No"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "non Fixer")] <- "No"

nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "yes")] <- "Yes"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "2")] <- "Yes"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "Y")] <- "Yes"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "y")] <- "Yes"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "yes, an N fixer")] <- "Yes"
nitfix3$`Nitrogen-fixation capacity`[which(nitfix3$`Nitrogen-fixation capacity` == "N-FIXER")] <- "Yes"

#get rid of entries with non comparable entries
nitfix4 <- nitfix3[which(nitfix3$`Nitrogen-fixation capacity` == "Yes" | nitfix3$`Nitrogen-fixation capacity` == "No"),]

nitfix_grouped <- nitfix4 %>% 
  dplyr::select(AccSpeciesName, `Nitrogen-fixation capacity`, Reference) %>% 
  group_by(AccSpeciesName, `Nitrogen-fixation capacity`) %>% 
  summarize(Reference = list(unique(Reference)))

unique(nitfix_grouped$`Nitrogen-fixation capacity`)

#what about second variable? nfc?
length(unique(nitfix2_grouped$AccSpeciesName))
nrow(nitfix2_grouped)

dup_species <- nitfix2_grouped %>% 
  group_by(AccSpeciesName) %>% 
  summarize(n=n()) %>% 
  dplyr::filter(n>1)

#dups are consistent for second variable - remove one of them each
nitfix2_grouped <- unique(nitfix2_grouped)

#add nitfix column to dataframe with references
colnames(fillable_df5)

fillable_df4 <- fillable_df3 %>% 
  left_join(nitfix_grouped, by = c("Species" = "AccSpeciesName"))

fillable_df5 <- fillable_df4 %>% 
  left_join(nitfix2_grouped, by = c("Species" = "AccSpeciesName"))


colnames(fillable_df5) <- c("Species", "ShadeTolerance_Quant", "ReferenceShade_Quant", 
                            "ShadeComment_Quant", "ShadeTolerance_Qual", "ReferenceShade_Qual",
                            "DroughtTolerance_Quant", "ReferenceDrought_Quant", "DroughtComment_Quant", 
                            "DroughtTolerance_Qual", "ReferenceDrought_Qual", "NitFix_Capacity", 
                            "ReferenceNitFix_Capacity", "NitFix_Clade", "ReferenceNitFix_Clade")

head(fillable_df5)

#save intermediate object
saveRDS(fillable_df5, "./data/predictors/intermediate_fillable_df5.rds")

#Do mycorrhizae columns

colnames(spread_data_sp)

#want columns 1:3, 7, 13, 16:25
working_data <- spread_data_sp[,c(1:3,7,13,16:25)]

#include this one
nrow(working_data[which(!is.na(working_data$`Mycorrhizal type`)),])

#looking at different types
length(unique(working_data$AccSpeciesName[which(!is.na(working_data$`Mycorrhizal type`))]))

#ignore this one for now
length(unique(working_data$AccSpeciesName[which(!is.na(working_data$`Original term for mycorrhizal type given by Selivanov`))]))

#incorporate these columns
length(unique(working_data$AccSpeciesName[which(!is.na(working_data$`Mycorrhiza: Arbuscular mycorrhizal (AM) fungi`))]))

#ignore this one for now
length(unique(working_data$AccSpeciesName[which(!is.na(working_data$`Mycorrhiza type according to Maherali`))]))

#standardize responses to myc type
unique(working_data$`Mycorrhizal type`)

working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "ecto")] <- "EM"
working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "ECTO")] <- "EM"
working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "Ecto")] <- "EM"
working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "Ectomycorrhiza")] <- "EM"
working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "ectomycorrhizal")] <- "EM"

working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "arbuscular")] <- "AM"

working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "Ectendo")] <- "EeM"
working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "ecto, ectendo")] <- "EM+EeM"

working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "Ericoid")] <- "ErM"

working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "0")] <- "NM"
working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "no")] <- "NM"
working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "Non")] <- "NM"

working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "EC/AM")] <- "EM+AM"
working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "AM + EM")] <- "EM+AM"

working_data$`Mycorrhizal type`[which(working_data$`Mycorrhizal type` == "Gram-AM")] <- "AM"

#and remove others
#get rid of entries with non comparable entries
#remove entries where ec?, 2 (mycorrhizal), VA, non-ectomycorrhizal, 
working_data2 <- working_data[-which(working_data$`Mycorrhizal type` == "2" | 
                                      working_data$`Mycorrhizal type` == "ec?"|
                                      working_data$`Mycorrhizal type` == "non-ectomycorrhizal"|
                                      working_data$`Mycorrhizal type` == "VA"),]


unique(working_data2$`Mycorrhizal type`)


#Type of mycorrhizae formed. 
#“AM” = arbuscular mycorrhizae; 
#“EM” = ectomycorrhizae; 
#“EeM” = ectendomycorrhizae; 
#“ErM” = ericoid mycorrhizae”; 
#“mycorrhizal” means mycorrhizae are present but type is unknown; 
#“NM” = non-mycorrhizal."

#summarize by species 
myc_data <- working_data2 %>% 
  dplyr::select(AccSpeciesName, `Mycorrhizal type`, Reference) %>% 
  dplyr::filter(!is.na(`Mycorrhizal type`)) %>% 
  group_by(AccSpeciesName) %>% #`Mycorrhizal type`
  summarize(MycType = toString(unique(`Mycorrhizal type`)), Reference = list(unique(Reference)))

###second variable
#get other myc columns
working_data3 <- spread_data_sp[,c(1:3,16:19)]
colnames(spread_data_sp)

colnames(working_data3) <- c("AccSpeciesName", "Reference", "Comment", "AM", "EM", "NM", "ErM")

working_data3 <- working_data3[rowSums(is.na(working_data3[,4:7]))!=4,]

combined_myc_cats <- working_data3 %>% 
  group_by(AccSpeciesName) %>% 
  mutate(row = row_number()) %>% 
  tidyr::pivot_longer(c(AM, EM, NM, ErM), names_to = "MycType_2", values_drop_na = TRUE) %>% 
  dplyr::select(-row) 
  #tidyr::pivot_wider(names_from = MycType_2, values_from = value) 

AM <- working_data3[,c(1:4)] %>% 
  dplyr::filter(!is.na(AM))

EM <- working_data3[,c(1:3,5)] %>% 
  dplyr::filter(!is.na(EM))

NM <- working_data3[,c(1:3,6)] %>% 
  dplyr::filter(!is.na(NM))

ErM <- working_data3[,c(1:3,7)] %>% 
  dplyr::filter(!is.na(ErM))

#merge dataframes back together
myc_2_data <- AM[,-3] %>% 
  bind_cols(EM[,4]) %>% 
  bind_cols(NM[,4]) %>% 
  bind_cols(ErM[,4])

#add myc_2_data and myc_data to dataframe

colnames(fillable_df5)

fillable_df6 <- fillable_df5 %>% 
  left_join(myc_data, by = c("Species" = "AccSpeciesName"))

fillable_df7 <- fillable_df6 %>% 
  left_join(myc_2_data, by = c("Species" = "AccSpeciesName"))


colnames(fillable_df7) <- c("Species", "ShadeTolerance_Quant", "ReferenceShade_Quant", 
                            "ShadeComment_Quant", "ShadeTolerance_Qual", "ReferenceShade_Qual",
                            "DroughtTolerance_Quant", "ReferenceDrought_Quant", "DroughtComment_Quant", 
                            "DroughtTolerance_Qual", "ReferenceDrought_Qual", "NitFix_Capacity", 
                            "ReferenceNitFix_Capacity", "NitFix_Clade", "ReferenceNitFix_Clade",
                            "MycType", "ReferenceMycType", "ReferenceMycBinaries", "Myc_AM", "Myc_EM",
                            "Myc_NM", "Myc_ErM")

head(fillable_df7)

#save intermediate object
saveRDS(fillable_df7, "./data/predictors/intermediate_fillable_df7.rds")


#reorganize dataframe so data first and refs later
colnames(fillable_df7)

trait_data <- fillable_df7[,c(1,2,5,7,10,12,14,16,19,20,21,22,4,9,3,6,8,11,13,15,17,18)]
colnames(trait_data)
head(trait_data)


saveRDS(trait_data, "./data/predictors/try_trait_data_cleaned.rds")

write.csv(trait_data[,c(1:12)], "./data/predictors/try_trait_data_cleaned.csv", row.names = FALSE)
