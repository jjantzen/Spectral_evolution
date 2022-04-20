#sorting out growth form data

library(taxize)
library(stringr)

#get data
growth_form_tidy <- readRDS("./data/predictors/growth_form_tidy.rds")

#get names in right format
growth_form_tidy$scientific_name <- word(growth_form_tidy$scientific_name, 1,2, sep=" ")

#read in taxonomy df
taxonomy <- readRDS("./data/predictors/vascan_taxonomy.rds")

#change colnames
colnames(taxonomy)[1] <- "species"
colnames(growth_form_tidy)[1] <- "species"

#merge dataframes
combo_df <- left_join(taxonomy, growth_form_tidy, by = "species")

#fill in missing data for Mahonia
combo_df$Herb[52] <- 0
combo_df$Shrub[52] <- 1
combo_df$Tree[52] <- 0
combo_df$`Tree/Shrub`[52] <- 1
combo_df$No_info[52] <- 0
combo_df$Woody[52] <- 1
combo_df[52,15] <- 0

#remove NA column
combo_df <- combo_df[,-15]

#save object
saveRDS(combo_df, "./data/predictors/combined_data.rds")

#save as dataframe
write.csv(combo_df, "./data/predictors/growth_form_taxonomy.csv")
