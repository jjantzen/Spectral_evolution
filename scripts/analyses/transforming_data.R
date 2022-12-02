#transform data for modeling (to see if it improves the distribution of residuals)


#read spectral data
data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

#get only spectra
spectra_only <- data_spectra$spectra

#transform spectra

#multiply by 100 and square root

spectra_only_trans <- sqrt(spectra_only*100)


spectra_only[1:10,1:10]
spectra_only_trans[1:10, 1:10]

#make new dataframe
data_spectra_trans <- data_spectra

data_spectra_trans$spectra <- spectra_only_trans

str(data_spectra_trans)

str(spectra_only_trans)

#save output
saveRDS(data_spectra_trans, "./data/for_analysis/myc_data_list_92sp_for_analysis_transformed.rds")
