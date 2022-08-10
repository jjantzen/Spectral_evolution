#play with spectral data - spectrolab version 0.0.8 needed
library(spectrolab)
library(ggplot2)

#read in spectral data
CABO_spectra <- readRDS("./data/all_spectra.rds")

spectra_df <- as.data.frame(CABO_spectra)

#save as csv file to import with version 0.0.10
write.csv(spectra_df, "./data/raw/Shan_spectral_data_saved.csv")

dim(CABO_spectra)

plot(CABO_spectra)

#get structure of spectra object
class(CABO_spectra)

CABO_spectra

#plot spectra all together

dim(CABO_spectra)

plot(CABO_spectra)


metadata <- meta(CABO_spectra)

spectra_df <- as.data.frame(CABO_spectra)
str(spectra_df)


unique(metadata$species)

unique(metadata$project)


CABO_spectra[which(CABO_spectra)]