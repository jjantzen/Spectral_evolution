#libraries
library(tidyr)
library(dplyr)

#read data

data <- read.csv("./data/leaf_spectra.csv", stringsAsFactors = FALSE)

#summary
head(data)
colnames(data)
rownames(data)

#remove bad quality leaves, those with only a single measurement, and those with only reflectance measured
#group by species and summarize by date measured or site
dates_sampled_species <- data %>% 
  filter(quality_leafs == "yes" & measurements > 1 & properties_measured == "both") %>% 
  group_by(scientific_name) %>% summarise(unique(date_measured))

sites_sampled_species <- data %>% 
  filter(quality_leafs == "yes" & measurements > 1 & properties_measured == "both") %>% 
  group_by(scientific_name) %>% summarise(unique(site_id))

#get unique date and unique sites for species into dataframe
dates_sampled_species_df <- as.data.frame(dates_sampled_species)
sites_sampled_species_df <- as.data.frame(sites_sampled_species)

saveRDS(dates_sampled_species_df, file = "./data/dates_sampled_species_df.rds")

#view comments on data
unique(data$sample_remarks)
unique(data$quality_leaves_comments)
unique(data$event_remarks)
unique(data$leaves_without_good_quality)
unique(data$properties_measured)

#get number of species measured total
length(unique(dates_sampled_species_df$scientific_name))
nspecies <- length(unique(sites_sampled_species_df$scientific_name))

#get number of species with multiple samples and with single samples
sp_samp_counts <- as.data.frame(dates_sampled_species %>% 
  count(scientific_name))

single_samp_species <- sp_samp_counts[which(sp_samp_counts$n == 1),]
mult_samp_species <- sp_samp_counts[which(sp_samp_counts$n > 1),]

count_mult <- nrow(mult_samp_species)
count_single <- nrow(single_samp_species)

unique(sites_sampled_species_df$`unique(site_id)`)
