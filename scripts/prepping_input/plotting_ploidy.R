#ploidy data

library("ggplot2")

ploidy <- read.csv("./data/predictors/ploidy_data_cleaned.csv")

colnames(ploidy)

ploidy


conifers <- c("Abies balsamea", "Larix laricina", "Picea abies", "Picea glauca", "Picea mariana", "Pinus banksiana", 
              "Pinus resinosa", "Pinus rigida", "Pinus strobus", "Polystichum munitum", "Phragmites australis", "Thuja occidentalis", "Tsuga canadensis")

ggplot(data = ploidy[-which(ploidy$ï..resolved_name %in% conifers),])+
  geom_point(aes(x = chromosome_count, y = cvalue))

