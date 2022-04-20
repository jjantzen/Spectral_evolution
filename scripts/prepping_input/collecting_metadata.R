#Early version - see Summarizing_data.rmd for updated version (as of Nov 23, 2020)

#get metadata and sampling schemes for all projects
#including dates, species names, number of individuals per species, locations and ???

#using the spectral objects to get sample_ids and what other dataframes to get additional metadata?

##########SET UP################
#load libraries
library(spectrolab)
library(tidyr)
library(dplyr)
library(lubridate)
library(reshape2)

#make presence_absence function
presence_absence <-function(data){
  tnt <- as.data.frame(data)
  # Create an ID column to hold the row index
  tnt$ID <- seq.int(nrow(tnt))
  out=dcast(tnt, ID ~data, length)
  
  # Get rid of the ID column from the matrix
  drops=c("ID")
  out=out[ , !(names(out) %in% drops)]
  
  return(out)
}

#make splitting species names function
splitting_species_names <- function(species_names_dataframe){
  sample_meta_names <- species_names_dataframe %>% 
    separate(species, c("genus", "specific_epithet", "extra0", "extra1", "extra2", "extra3", "extra4", "extra5", "extra6", "extra7", "extra8", "extra9", "extra10"), sep = " ", fill = "right", remove = FALSE)
  
  #deal with subspecific names
  sample_meta_names$subspecific <- NA
  sample_meta_names$subspecific[which(sample_meta_names$extra5 == "australis")] <- paste0(sample_meta_names$extra4[which(sample_meta_names$extra5 == "australis")], "_", sample_meta_names$extra5[which(sample_meta_names$extra5 == "australis")])
  
  sample_meta_names$subspecific[which(sample_meta_names$extra2 == "angustifolia")] <- paste0(sample_meta_names$extra1[which(sample_meta_names$extra2 == "angustifolia")], "_", sample_meta_names$extra2[which(sample_meta_names$extra2 == "angustifolia")])
  sample_meta_names$subspecific[which(sample_meta_names$extra2 == "vaginatum")] <- paste0(sample_meta_names$extra1[which(sample_meta_names$extra2 == "vaginatum")], "_", sample_meta_names$extra2[which(sample_meta_names$extra2 == "vaginatum")]) 
  
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "spissum")] <- paste0(sample_meta_names$extra0[which(sample_meta_names$extra1 == "spissum")], "_", sample_meta_names$extra1[which(sample_meta_names$extra1 == "spissum")])
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "strigosus")] <- paste0(sample_meta_names$extra0[which(sample_meta_names$extra1 == "strigosus")], "_", sample_meta_names$extra1[which(sample_meta_names$extra1 == "strigosus")]) 
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "rugosa")] <- paste0(sample_meta_names$extra0[which(sample_meta_names$extra1 == "rugosa")], "_", sample_meta_names$extra1[which(sample_meta_names$extra1 == "rugosa")])
  return(sample_meta_names)
}

##########READ DATA################
#read in spectral data object

spectra <- read.csv("https://web.fulcrumapp.com/shares/5c5290b884ab3734.csv") # Fulcrum data share link

spectral_data <- readRDS("./data/all_spectra.rds")

#get metadata object from spectra
metadata <- meta(spectral_data)
colnames(metadata)

#get species list from spectral data
s <- meta(spectral_data, "species", simplify = TRUE)

#get metadata from additional metadata file (downloaded separately from fulcrum)
meta_fulcrum <- spectra
meta_full <- read.csv("./data/metadata/leaf_spectra.csv", stringsAsFactors = FALSE)
pardo_names <- read.csv("./data/raw/Downloads_index_folders/2019-Pardo-MSc-UdeM-spectra-processed-project_all_combined.csv", stringsAsFactors = FALSE)
warren_names <- read.csv("./data/raw/Downloads_index_folders/SWA-Warren-spectra-processed-project_all_combined.csv", stringsAsFactors = FALSE)
phrag_names <- read.csv("./data/raw/Downloads_index_folders/2019-Phragmites-temporal-spectra-processed-project_all_combined.csv", stringsAsFactors = FALSE)

str(meta_fulcrum)
str(meta_full)
#check which sample IDS are not in my sample IDs
missing_my_meta <- meta_fulcrum[-which(meta_fulcrum$sample_id %in% meta_full$sample_id),]

#it includes samples from four categories
unique(missing_my_meta$status)

#there are no NA scientific names in Etienne's data
meta_fulcrum[which(is.na(meta_fulcrum$scientific_name)),]

#there are 126 names in Etienne's data
unique(meta_fulcrum$scientific_name)

nrow(missing_my_meta)
#there are samples of 71 species missing from my metadata that are in Etienne's (maybe not unique samples)
unique(missing_my_meta$scientific_name)

#there are samples from 8 projects missing in my metadata that are in Etienne's
unique(missing_my_meta$project)

#which project are the deleted ones missing from?
unique(missing_my_meta$project[which(missing_my_meta$status == "deleted")])

#which project are the NOT deleted ones missing from?
unique(missing_my_meta$project[-which(missing_my_meta$status == "deleted")])

unique(missing_my_meta$project[which(missing_my_meta$status == "verified")])

unique(missing_my_meta$project[which(missing_my_meta$status == "verified")])

unique(meta_full$project)

##########MATCH DATASETS################
#get metadata from meta_full for the sample_IDs in the spectra dataset
meta_for_spectra_EL <- meta_fulcrum[which(meta_fulcrum$sample_id %in% metadata$sample_id),]
meta_for_spectra <- meta_full[which(meta_full$sample_id %in% metadata$sample_id),]
meta_for_pardo <- pardo_names[which(pardo_names$sample_id %in% metadata$sample_id),]
meta_for_warren <- warren_names[which(warren_names$sample_id %in% metadata$sample_id),]
meta_for_phragmites <- phrag_names[which(phrag_names$sample_id %in% metadata$sample_id),]

str(meta_for_spectra)
nrow(meta_for_spectra_EL)
nrow(metadata)

#add species names for SWA-Warren and 2019-Pardo-Msc-UdeM from downloaded spectral files
missing_sample_ids <- metadata$sample_id[which(is.na(metadata$species))]

#there are scientific names in the metadata form
meta_for_spectra$scientific_name[which(meta_for_spectra$sample_id %in% missing_sample_ids)]

#make as character to allow filling NAs
metadata$species <- as.character(metadata$species)

for (i in 1:length(missing_sample_ids)){
  id <- missing_sample_ids[i]
  if(id %in% meta_for_pardo$sample_id){
    metadata$species[which(metadata$sample_id == id)] <- unique(meta_for_pardo$scientific_name[which(meta_for_pardo$sample_id == id)])
  } else if (id %in% meta_for_warren$sample_id){
    metadata$species[which(metadata$sample_id == id)] <- unique(meta_for_warren$scientific_name[which(meta_for_warren$sample_id == id)])
  } else if (id %in% meta_for_phragmites$sample_id){
    metadata$species[which(metadata$sample_id == id)] <- unique(meta_for_phragmites$scientific_name[which(meta_for_phragmites$sample_id == id)])
  }
}

#none still missing
still_missing_species_IDS <- metadata$sample_id[which(is.na(metadata$species))]

# to_split <- metadata
# class(to_split$species) <- "character"

split_names <- splitting_species_names(metadata)

sort(unique(split_names$specific_epithet))

split_names[which(split_names$specific_epithet == "Linnaeus"),]

#43021621 from csv is Quercus Linnaeus
#spelling errors in Pensylvanica and alleghaniensis
#only Linnaeus as authority

metadata$species[which("Solidago Linnaeus" %in% metadata$species)]
unique(metadata$species)


#99 is Salix L.

#add metadata back to spectra object
meta(spectral_data)
head(meta(spectral_data))
head(metadata)

nrow(metadata)
nrow(meta(spectral_data))
tail(metadata)
tail(meta(spectral_data))

#replace original metadata with the new (species names)
meta(spectral_data) <- metadata

saveRDS(spectral_data, "./data/tidy/spectra_spnames.rds")

head(metadata)


##############GET SPECIES NAMES FOR TREE###################
#get list of names for which there are spectra - after fixing meta of spectra
species_names_meta <- data.frame(species = unique(meta(spectral_data, "species", simplify = TRUE)))

#run for species names from spectra meta
sample_meta_names_spectra <- splitting_species_names(species_names_meta)
colnames(sample_meta_names_spectra)

#get rid of extra columns
species_names_spectral <- sample_meta_names_spectra[,-c(4:14)]
colnames(species_names_spectral)

species_names_spectral <- species_names_spectral %>% 
  mutate(genus_code = toupper(substr(genus, start = 1, stop = 2)), species_code = toupper(substr(specific_epithet, start = 1, stop = 2))) %>% 
  mutate(code = paste0(genus_code, species_code)) %>% 
  unite(species_name, genus, specific_epithet, subspecific, sep = " ", na.rm = TRUE)

#compare completeness
nrow(meta_for_spectra)
nrow(metadata)
nrow(meta_full)

#overview of data
str(meta_for_spectra)

#get list of samples/projects which need metadata or names (not in fulcrum metadata)
samples_with_metadata <- meta_full[which(meta_full$sample_id %in% metadata$sample_id),]
samples_without_metadata <- metadata[-which(metadata$sample_id %in% meta_full$sample_id),]

# #now all samples have species names - not necessary after adding from individual projects
# samples_without_speciesnames <- metadata[which(is.na(metadata$species)),]
# nrow(samples_without_speciesnames)
# unique(samples_without_speciesnames$project)
# 
# #check if names in metadata full if not in metadata from spectra 
# samples_nonames_in_metafull <- meta_full[which(meta_full$sample_id %in% samples_without_speciesnames$sample_id),]
# nrow(samples_nonames_in_metafull)
# unique(samples_nonames_in_metafull$project)
# 
# #those without names in other metadata either
# samples_nonames_notin_metafull <- samples_without_speciesnames[-which(samples_without_speciesnames$sample_id %in% meta_full$sample_id),]
# nrow(samples_nonames_notin_metafull)
# unique(samples_nonames_notin_metafull$project)
# 
# #those with species names (projects)
# samples_with_speciesnames <- metadata[which(!is.na(metadata$species)),]
# unique(samples_with_speciesnames$project)
# 
# #All 2019-Phragmites-temporal have species names in metadata_full but not metadata (from spectra)
# #Some SWA-Warren have species names in full but not all (none in spectral metadata)
# #No 2019-Pardo-MSc-UdeM have species names in either

#there are still samples without full metadata (even if they now have species names)
samples_without_metadata <- meta_full[which(meta_full$sample_id %in% metadata$sample_id),]


##########TIDY DATA################
#make dataframes with structure:
#sample_ID, proj_name, sp_name, date, siteID

sample_meta <- data.frame(sample_ID = meta_for_spectra$sample_id, 
                          project = meta_for_spectra$project, species = meta_for_spectra$scientific_name, 
                          date = meta_for_spectra$date_measured, siteID = meta_for_spectra$site_id, stringsAsFactors = FALSE)
str(sample_meta)
head(sample_meta)

#convert date to ymd from mdy
sample_meta$YMD <- mdy(sample_meta$date)

#separate species names into genus, species, authority

#run for metadata species names (from index page)
sample_meta_names_data <- splitting_species_names(sample_meta)


#check for subspecific names through extra columns
unique(sample_meta_names_data$extra1)
colnames(sample_meta_names_data)

#also need to add way to get just genus names for Quercus and Solidago

#get rid of extra columns due to subspecific names  (edit as needed)
sample_meta_names <- sample_meta_names_data[,-c(6:16)]

#get separate columns for year month day
sample_meta_dates <- sample_meta_names %>%
  mutate(year = lubridate::year(YMD), month = lubridate::month(YMD),
         day = lubridate::day(YMD))

#use presence_absence function to convert locations into separate columns of P/A matrix
locations_matrix <- as.data.frame(presence_absence(sample_meta_dates$siteID))

#join locations to date/names meta dataframe
sample_locations_dates <- cbind(sample_meta_dates, locations_matrix)
colnames(locations_matrix)

##########SUMMARIZE DATA################
#group by project and species to get number of individuals per species per project with P/A of locations per species
aggregated_project_data <- sample_locations_dates %>% 
  group_by(project, genus, specific_epithet, subspecific) %>% 
  summarise(n_ind = n(), month_range = paste0(min(month), "-", max(month)),  n_locations = length(unique(siteID)), across(c(9:47), max)) #, n_species = n(unique(species))

aggregated_species_data <- sample_locations_dates %>% 
  group_by(genus, specific_epithet, subspecific) %>% 
  summarise(n_ind = n(), month_range = paste0(min(month), "-", max(month)), n_locations = length(unique(siteID)), across(c(10:49), max)) #, n_species = n(unique(species))

#get sample-specific metadata
sample_detailed_metadata <- subset(sample_locations_dates, select=-c(6,13:51))
colnames(sample_locations_dates)


#get list of all species
species_all <- aggregated_species_data %>% 
  unite(species_names, genus, specific_epithet, subspecific, sep = "_", na.rm = TRUE) %>% 
  select(species_names)


##########IDENTIFY BEST PROJECTS/SPECIES/SITES################

#get stats on various thresholds of species records
#get list of locations where species are found in >2 locations
n2_locations <- aggregated_species_data %>% 
  melt(id = c(1:6)) %>% 
  filter(n_locations > 1 & value == 1) %>% 
  select(c(1:7)) 

n2_locations_matrix <- as.data.frame(presence_absence(n2_locations$variable))

n2_dataframe <- cbind(n2_locations, n2_locations_matrix)
ncol(n2_dataframe)

#get locations with the most taxa
num_sp_per_location <- summarize(n2_dataframe, across(c(8:43), sum))

loc_4sp <- num_sp_per_location[which(num_sp_per_location >3)]

loc_3sp <- num_sp_per_location[which(num_sp_per_location >2)]

############GET DATA FOR SAVING##############

#save list of all species names
write.csv(aggregated_species_data, "./data/tidy/species_data_list.csv", row.names = FALSE)

#save list of projects with data
write.csv(aggregated_project_data, "./data/tidy/project_data_list.csv", row.names = FALSE)

#save list of samples with data
write.csv(sample_detailed_metadata, "./data/tidy/sample_data_list.csv", row.names = FALSE)

#save list of taxa/samples missing metadata
write.csv(samples_without_metadata, "./data/tidy/missing_data_list.csv", row.names = FALSE)

#save list of locations with more than 2 taxa
write.csv(loc_3sp, "./data/tidy/locations_3sp_list.csv", row.names = FALSE)

#save list of all taxa for phylogeny
write.csv(species_names_spectral, "./data/tidy/species_names.csv", row.names = FALSE)


# check out data
# #for example, get range of months for Acer saccharum Marshall
# unique(sample_meta_dates$month[which(sample_meta_dates$species == "Acer saccharum Marshall")])
# 
# unique(sample_meta_dates$YMD[which(sample_meta_dates$species == "Acer saccharum Marshall")])


######################################
#metadata from spectra themselves (limited info)

#get number of species
n_sp <- length(unique(s))

#get number of individuals per species
n_indXsp <- metadata %>% 
  count(as.character(species)) %>% 
  arrange()

#count individuals by project
n_indXproj <- metadata %>% 
  count(as.character(project))

#count individuals by project
n2_indXproj <- metadata %>% 
  group_by(project) %>% 
  summarise(n_ind = n())

#ind per species by project
n_spXproj <- metadata %>% 
  group_by(project, species) %>% 
  summarise(n_ind = n())

#group individuals by project - reshape dataframe to be project names as columns??
indXproj <- metadata %>% 
  mutate()

metadata[which(metadata$project == "2019-Pardo-MSc-UdeM"),]

           