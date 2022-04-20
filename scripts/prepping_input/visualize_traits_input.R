#plot characters onto tree to visualize data we have

#libraries
library(dplyr)
library(ape)
library(phytools)
library(tidyr)
library(ggplot2)
library(ggtree)
#library(ggnewscale)
library(pals)
library(randomcoloR)
library(arules)
library(plyr)

#write function

splitting_names <- function(species_names_dataframe){
  sample_meta_names <- species_names_dataframe %>% 
    separate(tip.names, c("genus", "specific_epithet", "OTT_ID", "extra1", "extra2", "extra3", "extra4", "extra5", "extra6"), sep = "_", fill = "right", remove = FALSE)
  
  #deal with subspecific names
  sample_meta_names$subspecific <- NA
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "australis")] <- paste0(sample_meta_names$OTT_ID[which(sample_meta_names$extra5 == "australis")], " ", sample_meta_names$extra1[which(sample_meta_names$extra1 == "australis")])
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "angustifolia")] <- paste0(sample_meta_names$OTT_ID[which(sample_meta_names$extra1 == "angustifolia")], " ", sample_meta_names$extra1[which(sample_meta_names$extra1 == "angustifolia")])
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "vaginatum")] <- paste0(sample_meta_names$OTT_ID[which(sample_meta_names$extra1 == "vaginatum")], " ", sample_meta_names$extra1[which(sample_meta_names$extra1 == "vaginatum")]) 
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "spissum")] <- paste0(sample_meta_names$OTT_ID[which(sample_meta_names$extra1 == "spissum")], " ", sample_meta_names$extra1[which(sample_meta_names$extra1 == "spissum")])
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "strigosus")] <- paste0(sample_meta_names$OTT_ID[which(sample_meta_names$extra1 == "strigosus")], " ", sample_meta_names$extra1[which(sample_meta_names$extra1 == "strigosus")]) 
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "rugosa")] <- paste0(sample_meta_names$OTT_ID[which(sample_meta_names$extra1 == "rugosa")], " ", sample_meta_names$extra1[which(sample_meta_names$extra1 == "rugosa")])
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "nigrum")] <- paste0(sample_meta_names$OTT_ID[which(sample_meta_names$extra1 == "nigrum")], " ", sample_meta_names$extra1[which(sample_meta_names$extra1 == "nigrum")])
  
  sample_meta_names$OTT_ID[which(sample_meta_names$OTT_ID == "-species")] <- sample_meta_names$extra5[which(sample_meta_names$OTT_ID == "-species")]
  return(sample_meta_names)
}

splitting_species_names <- function(species_names_dataframe){
  sample_meta_names <- species_names_dataframe %>% 
    #either species or scientific_name
    separate(species, c("genus", "specific_epithet", "extra0", "extra1", "extra2", "extra3", "extra4", "extra5", "extra6", "extra7", "extra8", "extra9", "extra10"), sep = " ", fill = "right", remove = FALSE)
  
  #deal with subspecific names
  sample_meta_names$subspecific <- NA
  sample_meta_names$subspecific[which(sample_meta_names$extra5 == "australis")] <- paste0(sample_meta_names$extra4[which(sample_meta_names$extra5 == "australis")], " ", sample_meta_names$extra5[which(sample_meta_names$extra5 == "australis")])
  
  sample_meta_names$subspecific[which(sample_meta_names$extra2 == "angustifolia")] <- paste0(sample_meta_names$extra1[which(sample_meta_names$extra2 == "angustifolia")], " ", sample_meta_names$extra2[which(sample_meta_names$extra2 == "angustifolia")])
  sample_meta_names$subspecific[which(sample_meta_names$extra2 == "vaginatum")] <- paste0(sample_meta_names$extra1[which(sample_meta_names$extra2 == "vaginatum")], " ", sample_meta_names$extra2[which(sample_meta_names$extra2 == "vaginatum")]) 
  
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "spissum")] <- paste0(sample_meta_names$extra0[which(sample_meta_names$extra1 == "spissum")], " ", sample_meta_names$extra1[which(sample_meta_names$extra1 == "spissum")])
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "strigosus")] <- paste0(sample_meta_names$extra0[which(sample_meta_names$extra1 == "strigosus")], " ", sample_meta_names$extra1[which(sample_meta_names$extra1 == "strigosus")]) 
  sample_meta_names$subspecific[which(sample_meta_names$extra1 == "rugosa")] <- paste0(sample_meta_names$extra0[which(sample_meta_names$extra1 == "rugosa")], " ", sample_meta_names$extra1[which(sample_meta_names$extra1 == "rugosa")])
  return(sample_meta_names)
}


#read tree
spectra_tree <- read.tree("./data/tol_tree.tre")
spectra_tree

#read data
trait_data <- read.csv("./data/tidy/merged_meta_spectra_trimmed.csv", stringsAsFactors = FALSE)
growth_form <- readRDS("./data/tidy/growth_form_tidy.rds")
taxonomy <- read.csv("./data/tidy/taxonomy_table_higher_class.csv", stringsAsFactors = FALSE)

#match labels
trait_data <- trait_data %>% 
  unite(species_names, genus, specific_epithet, subspecific, sep = " ", na.rm = TRUE)

growth_form <- splitting_species_names(growth_form)
growth_form <- growth_form %>% 
  unite(species_names, genus, specific_epithet, subspecific, sep = " ", na.rm = TRUE)
growth_form <- growth_form[,c(2,14:19)]
colnames(growth_form)
growth_form

tip_names <- data.frame(tip.names = spectra_tree$tip.label, stringsAsFactors = FALSE)
tip_names_fixed <- splitting_names(tip_names)
tip_names_fixed <- tip_names_fixed %>% 
  unite(species_names, genus, specific_epithet, subspecific, sep = " ", na.rm = TRUE)
spectra_tree$tip.label <- tip_names_fixed$species_names

#need to change function in between line
taxonomy <- splitting_species_names(taxonomy)
taxonomy <- taxonomy %>% 
  unite(species_names, genus, specific_epithet, subspecific, sep = " ", na.rm = TRUE)
taxonomy <- taxonomy[,c(2,14:17)]
colnames(taxonomy)

trait_data$species_names
tip_names_fixed$species_names

tip_names_fixed[which(tip_names_fixed$species_names %in% trait_data$species_names == FALSE),]
unique(trait_data[which(trait_data$species_names %in% tip_names_fixed$species_names == FALSE),28])

#synonyms:
#Kalmia angustifolia vs Kalmia angustifolia var. angustifolia
#Acer saccharum subsp. nigrum vs Acer nigrum
#Phragmites australis vs Phragmites australis subsp. australis
#Solidago fistulosa vs Solidago altissima

#change the synonyms in the tree because fewer changes
spectra_tree$tip.label[which(spectra_tree$tip.label == "Kalmia angustifolia")] <- "Kalmia angustifolia var. angustifolia"
spectra_tree$tip.label[which(spectra_tree$tip.label == "Acer saccharum subsp. nigrum")] <- "Acer nigrum"
spectra_tree$tip.label[which(spectra_tree$tip.label == "Phragmites australis")] <- "Phragmites australis subsp. australis"
spectra_tree$tip.label[which(spectra_tree$tip.label == "Solidago fistulosa")] <- "Solidago altissima"

#get smaller trait dataset that matches tree
traits_for_mapping <- trait_data[which(trait_data$species_names %in% spectra_tree$tip.label),]
nrow(traits_for_mapping)

traits_summary <- traits_for_mapping %>% 
  group_by(species_names) %>% 
  select(project, species_names, PL.site_id, day_measured, month_measured, year_measured) %>% 
  summarise(n_sample_dates = length(unique(paste(day_measured, month_measured))), earliest_month = min(month_measured), latest_month = max(month_measured), month_range = paste0(min(month_measured), "-", max(month_measured)), duration_sampling = (max(month_measured)-min(month_measured)), n_locations = length(unique(PL.site_id)))

#compare names for taxonomy table and tree
taxonomy[which(taxonomy$species_names %in% spectra_tree$tip.label == FALSE),]
spectra_tree$tip.label[which(spectra_tree$tip.label %in% taxonomy$species_names == FALSE)]

#in taxonomy:
#duplicate Eriophorum vaginatum subsp. vaginatum and change subsp to spissum
#add column for Berberis aquifolium
taxonomy <- rbind(taxonomy, c("Berberis aquifolium", "Equisetopsida", "Ranunculanae", "Ranunculales", "Berberidaceae"))
taxonomy <- rbind(taxonomy, taxonomy[34,])
taxonomy$species_names[105] <- "Eriophorum vaginatum subsp. spissum"
taxonomy <- rbind(taxonomy, c("Tsuga canadensis", "Equisetopsida", "Pinidae", "Pinales", "Pinaceae"))
#get rid of duplicate solidago
taxonomy <- taxonomy[-88,]

#compare names for growth form table and tree
growth_form[which(growth_form$species_names %in% spectra_tree$tip.label == FALSE),]
spectra_tree$tip.lable[which(spectra_tree$tip.label %in% growth_form$species_names == FALSE)]
#get rid of duplicate solidago
growth_form <- growth_form[-101,]

#compare names for trait summary and tree
traits_summary[which(traits_summary$species_names %in% spectra_tree$tip.label == FALSE),]
spectra_tree$tip.label[which(spectra_tree$tip.label %in% traits_summary$species_names == FALSE)]

#merge data together for easier plotting
combo_data <- merge(tip_names_fixed, traits_summary, by = "species_names", all = TRUE)
nrow(combo_data)
combo_data <- merge(combo_data, growth_form, by = "species_names", all  = TRUE)
nrow(combo_data)
combo_data <- merge(combo_data, taxonomy, by = "species_names", all = TRUE)
nrow(combo_data)
colnames(combo_data)

combo_data <- combo_data[,c(1,10:25)]
combo_data$species_names

#give rownames for 
rownames(combo_data) <- combo_data$species_names

#trim to tree
combo_data <- combo_data[which(combo_data$species_names %in% spectra_tree$tip.label),]

saveRDS(combo_data, "./data/tidy/combo_data_for_plotting.rds")

combo_data <- readRDS("./data/tidy/combo_data_for_plotting.rds")

#make columns right classes
combo_data$n_sample_dates <- as.factor(combo_data$n_sample_dates)
levels(combo_data$n_sample_dates) <- list("1-3"=c("1", "2", "3"), "4-7"=c("4", "5", "6", "7"), "8-11" = c("8", "9", "11"), "13-15" = c("13", "14", "15"), "17-20" = c("17", "20"))

combo_data$duration_sampling <- as.factor(combo_data$duration_sampling)

combo_data$n_locations <- as.factor(combo_data$n_locations)
levels(combo_data$n_locations) <- list("1"=c("1"), "2-4" =c("2", "3", "4"), "5-7"=c("5", "6", "7"), "8-10" = c("8", "10"))

combo_data$Herb <- as.factor(combo_data$Herb)
levels(combo_data$Herb) <- list("No"=c("0"), "Yes" =c("1"))

combo_data$Shrub <- as.factor(combo_data$Shrub)
levels(combo_data$Shrub) <- list("No"=c("0"), "Yes" =c("1"))

combo_data$Tree <- as.factor(combo_data$Tree)
levels(combo_data$Tree) <- list("No"=c("0"), "Yes" =c("1"))

combo_data$Vine <- as.factor(combo_data$Vine)
levels(combo_data$Vine) <- list("No"=c("0"), "Yes" =c("1"))

combo_data$No_info <- as.factor(combo_data$No_info)
levels(combo_data$No_info) <- list("No"=c("0"), "Yes" =c("1"))

combo_data$Woody <- as.factor(combo_data$Woody)
levels(combo_data$Woody) <- list("No"=c("0"), "Yes" =c("1"))

combo_data$clade <- as.factor(combo_data$clade)

combo_data$class <- as.factor(combo_data$class)

combo_data$order <- as.factor(combo_data$order)
levels(combo_data$order) <- list("A_orders"=c("Apiales", "Asparagales", "Asterales"), 
                                 "C_orders" =c("Caryophyllales", "Cornales", "Cupressales"), "DE_orders" = c("Dipsacales", "Ericales"), "FG_orders" = c("Fabales", "Fagales", "Gentianales"), "LM_orders" = c("Lamiales", "Malpighiales", "Malvales", "Myrtales"), "P_orders" = c("Pinales", "Poales"), "R_orders" = c("Ranunculales", "Rosales"), "SV_orders" = c("Sapindales", "Vitales"))

combo_data$family <- as.factor(combo_data$family)
levels(combo_data$family) <- list("A_families"=c("Anacardiaceae", "Apiaceae", "Apocynaceae", "Asparagaceae", "Asteraceae"), 
                                 "B_families" =c("Berberidaceae", "Betulaceae"), "C_families" = c("Cannabaceae", "Caprifoliaceae", "Cornaceae", "Cupressaceae", "Cyperaceae"),
                                 "EF_families" = c("Ericaceae", "Fabaceae", "Fagaceae"), "JLM_families" = c("Juglandaceae", "Lythraceae", "Malvaceae", "Montiaceae"), 
                                 "OP_families" = c("Oleaceae", "Pinaceae", "Poaceae"), "RS_families" = c("Rhamnaceae", "Rosaceae", "Salicaceae", "Sapindaceae"), 
                                 "TUV_families" = c("Typhaceae", "Ulmaceae", "Vitaceae"))


#give tree branch lengths
spectra_tree_br <- compute.brlen(spectra_tree, method = "Grafen", power = 1)

##set the colours for the different variables so that colours match
#get list of 8 colours (max num states)
myColor <- randomcoloR::distinctColorPalette(k = 8)
plot(myColor)

#get subset of columns for plotting
data_subset <- combo_data[,c(2,6:13,15:17)]
#colnames(combo_data)

#create empty list of right length
empty_list <- vector(mode = "list", length = length(colnames(data_subset)))

#get levels of variables as names and get associated length of colours as values
for (i in 1:ncol(data_subset)){
  value <- myColor[1:length(levels(data_subset[,i]))]
  names <- levels(data_subset[,i])
  combined <- setNames(value, names)
  empty_list[i]<- list(combined)
  names(empty_list)[i] <- colnames(data_subset[i])
}

colours <- empty_list

#first rescale size of tree to fit better
scale <- 0.5
spectra_tree_br$edge.length<-
  spectra_tree_br$edge.length/max(nodeHeights(spectra_tree_br)[,2])*scale

#make plot
jpeg("./figures/exploratory_data_tree.jpg", width = 2000, height = 1500)
object <- plotTree.datamatrix(spectra_tree_br, data_subset, length=5,ftype="i", fsize = 1.2, colors = colours, legend = 2, x.lim = 40, y.lim = 100)#args.plotTree = list(), # col = myColor, 
#par(mar=c(100,100,100,100))
x<-object$end.x+diff(par()$usr[1:2])*0.01
y<-Ntip(spectra_tree_br)
for(i in 1:ncol(data_subset)){
  text(x,y,colnames(data_subset)[i],pos=4,cex=1.2,offset=0)
  add.simmap.legend(colors=object$colors[[i]],shape="square",
                    prompt=FALSE,x=x,y=y-2*strheight("W")*1.2,fsize=1.2)
  y<-y-1.5*0.9*strheight("W")*(length(object$colors[[i]])-1)-6
}
dev.off()

#maybe want to switch around yes and no colours for more intuitive scheme

#plot time sampling

time_to_plot <- readRDS("./data/tidy/time_to_plot.rds")

time_to_plot$discrete_times_cluster <- discretize(time_to_plot$day_of_year, method = "cluster")
time_to_plot$discrete_times_freq <- discretize(time_to_plot$day_of_year, method = "frequency")
time_to_plot$discrete_times_interval <- discretize(time_to_plot$day_of_year, method = "interval")

unique(time_to_plot$discrete_times_cluster)
unique(time_to_plot$discrete_times_freq)

time_points <- c(150, 191, 215, 281)

as.Date(time_points, origin = "2018-01-01")

time_to_plot$discrete_times_cluster <- revalue(time_to_plot$discrete_times_cluster, c("[150,185)"="May 31 to July 4", "[185,223)" = "July 5 to August 11", "[223,281]"="August 12 to October 9"))
time_to_plot$discrete_times_interval <- revalue(time_to_plot$discrete_times_interval, c("[150,194)"="May 31 to July 13", "[194,237)" = "July 14 to August 25", "[237,281]"="August 26 to October 9"))
time_to_plot$discrete_times_freq <- revalue(time_to_plot$discrete_times_freq, c("[150,191)"="May 31 to July 10", "[191,215)" = "July 11 to August 3", "[215,281]"="August 4 to October 9"))

#For clustered: May 31 to July 4; July 5 to August 11, August 12 to Oct 9
#For frequency: May 31 to July 10; July 11 to August 3, August 4 to Oct 9
#For interval: May 31 to July 13; July 14 to August 25, August 26 to Oct 9



jpeg("./figures/time_sampling.jpg")
ggplot(time_to_plot)+
  geom_point(aes(x = as.factor(day_measured), y = species_name))+
  facet_wrap(vars(month_measured), scales = "free_x", ncol = 2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

jpeg("./figures/time_discrete.jpg", 20, 15, unit = "cm", res = 400)
time_to_plot %>% 
  pivot_longer(discrete_times_cluster:discrete_times_interval, names_to = "groupings", values_to = "response") %>%
  ggplot(aes(x = response, y = species_name, colour = groupings)) +
  facet_wrap(vars(groupings), scales = "free_x", ncol = 3) +
  geom_point() +
  labs(x = "Discrete time periods", y = "Species")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


jpeg("./figures/time_sampling.jpg", 20, 15, unit = "cm", res = 400)
ggplot(time_to_plot)+
  geom_point(aes(x = day_of_year, y = species_name))+
  #geom_point(aes(x = discrete_times_interval, y = species_name))+
  #facet_wrap(vars(month_measured), scales = "free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

unique(time_to_plot$discrete_times_interval)

time_points <- c(150, 194, 237, 281)

as.Date(time_points, origin = "2018-01-01")

#For clustered: May 31 to July 4; July 5 to August 11, August 12 to Oct 9
#For frequency: May 31 to July 10; July 11 to August 3, August 4 to Oct 9
#For interval: May 31 to July 13; July 14 to August 25, August 26 to Oct 9

