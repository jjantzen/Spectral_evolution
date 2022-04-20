#getting taxonomy sorted out

library(taxize)
library(stringr)

#get list of names
pruned_ang_con_tree <- readRDS("./data/ang_conifer/rds/pruned_ang_con_tree.rds")

species_names <- pruned_ang_con_tree$tip.label

#retrieve taxonomy from vascan
vascan_output <- vascan_search(species_names, format = "json", raw = FALSE)
str(vascan_output)

#make dataframe for output
taxonomy_df <- as.data.frame(matrix(ncol = 2, nrow = 98))
colnames(taxonomy_df) <- c("species", "highertaxonomy")

for (i in 1:(length(vascan_output))){
  species <- vascan_output[i][[1]]$matches[[1]][[1]]$acceptednameusage
  taxonomy <- vascan_output[i][[1]]$matches[[1]][[1]]$higherclassification
  if (vascan_output[i][[1]]$matches[[1]][[1]]$taxonomicstatus == "synonym"){
    taxonomy_df$species[i] <- species
    taxonomy_df$highertaxonomy[i] <- "synonym"
  } else {
    taxonomy_df$species[i] <- species
    taxonomy_df$highertaxonomy[i] <- taxonomy
  }
}

#search  mahonia separately 
mahonia <- vascan_search("Mahonia aquifolium", format = "json", raw = FALSE)

taxonomy_df$highertaxonomy[which(taxonomy_df$species == "Mahonia aquifolium")] <- mahonia[[1]]$matches[[1]][[1]]$higherclassification

#split names and taxonomy on semicolon
taxonomy_df$species <- word(taxonomy_df$species, 1,2, sep=" ")

split_taxonomy <- str_split(taxonomy_df$highertaxonomy, pattern = ";", simplify = TRUE)

split_taxonomy <- as.data.frame(split_taxonomy, stringsAsFactors = FALSE)

split_taxonomy <- cbind(taxonomy_df$species, split_taxonomy)

#rename columns with taxonomic group
unique(split_taxonomy[,3])
tail(split_taxonomy)
#spread column

#select order values
split_taxonomy_orders <- filter(split_taxonomy, grepl("ales", V3))

#select subclass values
split_taxonomy_superorder <- filter(split_taxonomy, grepl("anae", V3))

#reformat
split_taxonomy_superorder_renamed <- split_taxonomy_superorder %>%
  select(-c(V10, V11)) %>%
  rename("Class" = V1, "Subclass" = V2, "Superorder" = V3, "Order" = V4, "Family" = V5)

split_taxonomy_orders_renamed <- split_taxonomy_orders %>%
  select(-c(V6, V7, V8, V9, V10)) %>%
  rename("Class" = V1, "Subclass" = V2, "Order" = V3, "Family" = V4, "Genus" = V5)

combo_split_taxonomy <- split_taxonomy_superorder_renamed %>% 
  bind_rows(split_taxonomy_orders_renamed)

combo_split_taxonomy[c(80:100),]

#select subfamily values
split_taxonomy_subfamily <- filter(combo_split_taxonomy, grepl("oideae", V6))

remainder <- filter(combo_split_taxonomy, !grepl("oideae", V6))

#reformat
split_taxonomy_subfamily_renamed <- split_taxonomy_subfamily %>%
  select(-V11) %>%
  rename("Subfamily" = V6)

#reformat
remainder_renamed <- remainder %>%
  select(-c(V9, V11)) 

combo_split_taxonomy2 <- split_taxonomy_subfamily_renamed %>% 
  bind_rows(remainder_renamed)

#select tribe values
split_taxonomy_tribe1 <- filter(combo_split_taxonomy2, grepl("eae", V6))
remainder2 <- filter(combo_split_taxonomy2, !grepl("eae", V6))

split_taxonomy_tribe1_renamed <- split_taxonomy_tribe1 %>%
  select(-V9) %>%
  rename("Tribe" = V6)

combo_split_taxonomy3 <- split_taxonomy_tribe1_renamed %>% 
  bind_rows(remainder2)

split_taxonomy_tribe2 <- filter(combo_split_taxonomy3, grepl("eae", V7))
remainder3 <- filter(combo_split_taxonomy3, !grepl("eae", V7))

split_taxonomy_tribe2_renamed <- split_taxonomy_tribe2 %>%
  select(-Tribe) %>%
  rename("Tribe" = V7)

combo_split_taxonomy4 <- split_taxonomy_tribe2_renamed %>% 
  bind_rows(remainder3)

#select subtribe values
split_taxonomy_subtribe <- filter(combo_split_taxonomy4, grepl("inae", V8))

remainder3 <- filter(combo_split_taxonomy4, !grepl("inae", V8))

split_taxonomy_subtribe_renamed <- split_taxonomy_subtribe %>%
  select(-c(V6, V7, V9, Genus)) %>%
  rename("Subtribe" = V8)

remainder3_renamed <- remainder3 %>% 
  select(-c(Genus, V9, V7, V6, V8))

combo_split_taxonomy5 <- split_taxonomy_subtribe_renamed %>% 
  bind_rows(remainder3_renamed)

#final taxonomy dataset - to merge back with species names
final_taxonomy <- combo_split_taxonomy5

saveRDS(final_taxonomy, "./data/predictors/vascan_taxonomy.rds")

nrow(final_taxonomy)
