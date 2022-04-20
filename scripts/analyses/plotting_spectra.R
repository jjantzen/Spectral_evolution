#plotting real spectra
library(spectrolab)
library(ggplot2)

#improve the colour choices


#plot actual spectra by myc association

#read spectra
data_spectra <- readRDS("./data/for_analysis/myc_data_list_for_analysis.rds")

#convert to spectra object
spec_for_plot <- as_spectra(data_spectra$spectra)

#get the metadata in right format
meta(spec_for_plot)$myc <- data_spectra$myc

#get means
AM_mean <- mean(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),])

EM_mean <- mean(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),])

colours_touse <- c("#008837", "#E69F00", "#7b3294")

colours_trans <- c(make.transparent("#008837", 0.15), make.transparent("#E69F00", 0.15), make.transparent("#7b3294", 0.2))

#colours_trans <- c(make.transparent("#008837", 0.3), make.transparent("#E69F00", 0.3), make.transparent("#7b3294", 0.1))

#plot means
plot(AM_mean, lwd = 0.75, lty = 1, col = "blue", main = "Spectra by mycorrhizal association")
plot(EM_mean, lwd = 0.75, lty = 1, col = "red", add = TRUE)

jpeg("./output/spectra_by_myc.jpg", res = 600, width = 10, height = 10, units = "in")
plot(AM_mean, lwd = 0.75, lty = 1, col = rgb(1,0,0), ylim = c(0,0.6), main = "Spectra by mycorrhizal association (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance")
plot(EM_mean, lwd = 0.75, lty = 1, col = rgb(0,0,1), add = TRUE)
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], total_prob = 0.95, col = rgb(1,0,0,0.1), border = FALSE, add = TRUE) #main = "Spectra by mycorrhizal association (80% quantile)"
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], total_prob = 0.95, col = rgb(0,0,1,0.1), border = FALSE, add = TRUE)
legend(2000, 0.6, legend=c("AM", "EM"), fill=c(rgb(1,0,0,0.25), rgb(0,0,1,0.25)), cex=2)#lty=1,col
dev.off()

jpeg("./output/spectra_by_myc.jpg", res = 600, width = 10, height = 10, units = "in")
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], ylim = c(0,0.6), total_prob = 0.95, col = colours_trans[1], border = FALSE, main = "Spectra by mycorrhizal association (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association (80% quantile)"
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], total_prob = 0.95, col = colours_trans[3], border = FALSE, add = TRUE)
plot(AM_mean, lwd = 1, lty = 1.5, col = colours_touse[1], add = TRUE)
plot(EM_mean, lwd = 1, lty = 1.5, col = colours_touse[3], add = TRUE)
legend(2000, 0.6, legend=c("AM", "EM"), fill=c(colours_trans[c(1,3)]), cex=2)#lty=1,col
dev.off()

#plot actual spectra by habit
#get the metadata in right format
meta(spec_for_plot)$gf <- data_spectra$gf

#get means
Tree_mean <- mean(spec_for_plot[which(spec_for_plot$meta$gf == "Tree"),])

Herb_mean <- mean(spec_for_plot[which(spec_for_plot$meta$gf == "Herb"),])

Shrub_mean <- mean(spec_for_plot[which(spec_for_plot$meta$gf == "Shrub"),])

colours_trans
colours_touse

#plot means
plot(Tree_mean, lwd = 0.75, lty = 1, col = "blue", main = "Spectra by mycorrhizal association")
plot(Herb_mean, lwd = 0.75, lty = 1, col = "red", add = TRUE)
plot(Shrub_mean, lwd = 0.75, lty = 1, col = "green", add = TRUE)

jpeg("./output/spectra_by_gf.jpg", res = 600, width = 10, height = 10, units = "in")
plot_quantile(spec_for_plot[which(spec_for_plot$meta$gf == "Tree"),], total_prob = 0.95, col = colours_trans[1], border = FALSE, ylim = c(0,0.6), main = "Spectra by growth form (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association (80% quantile)"
plot_quantile(spec_for_plot[which(spec_for_plot$meta$gf == "Herb"),], total_prob = 0.95, col = colours_trans[3], border = FALSE, add = TRUE)
plot_quantile(spec_for_plot[which(spec_for_plot$meta$gf == "Shrub"),], total_prob = 0.95, col = colours_trans[2], border = FALSE, add = TRUE)
plot(Tree_mean, lwd = 0.75, lty = 1, col = colours_touse[1], add = TRUE)
plot(Herb_mean, lwd = 0.75, lty = 2, col = colours_touse[3], add = TRUE)
plot(Shrub_mean, lwd = 0.75, lty = 3, col = colours_touse[2], add = TRUE)
legend(1900, 0.6, legend=c("Tree", "Shrub", "Herb"), fill=c(colours_trans), cex=2)#lty=1,col
dev.off()

#plot actual spectra by leaf persistence
#get the metadata in right format
meta(spec_for_plot)$lp <- data_spectra$lp

#get means
deciduous_mean <- mean(spec_for_plot[which(spec_for_plot$meta$lp == "deciduous"),])

evergreen_mean <- mean(spec_for_plot[which(spec_for_plot$meta$lp == "evergreen"),])

#plot means
plot(deciduous_mean, lwd = 0.75, lty = 1, col = "blue", main = "Spectra by mycorrhizal association")
plot(evergreen_mean, lwd = 0.75, lty = 1, col = "red", add = TRUE)

jpeg("./output/spectra_by_lp.jpg", res = 600, width = 10, height = 10, units = "in")
plot_quantile(spec_for_plot[which(spec_for_plot$meta$lp == "deciduous"),], total_prob = 0.95, col = colours_trans[1], border = FALSE, ylim = c(0,0.6), main = "Spectra by leaf persistence (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association (80% quantile)"
plot_quantile(spec_for_plot[which(spec_for_plot$meta$lp == "evergreen"),], total_prob = 0.95, col = colours_trans[3], border = FALSE, add = TRUE)
plot(deciduous_mean, lwd = 0.75, lty = 1, col = colours_touse[1], add = TRUE)
plot(evergreen_mean, lwd = 0.75, lty = 1, col = colours_touse[3], add = TRUE)
legend(1700, 0.6, legend=c("Deciduous", "Evergreen"), fill=c(colours_trans[c(1,3)]), cex=2)#lty=1,col
dev.off()

####for larger dataset and other predictors
#read spectra
traits <- readRDS("./data/for_analysis/trait_data_100_taxa.rds")
spectra <- readRDS("./data/for_analysis/spectra_not_reordered_to_tree.rds")

#convert to spectra object
spec_for_plot_2 <- as_spectra(spectra)

colnames(traits)
#get the metadata in right format
meta(spec_for_plot_2)$drought <- traits$Drought_bin
meta(spec_for_plot_2)$woody <- traits$Woody
meta(spec_for_plot_2)$shade <- traits$Shade
meta(spec_for_plot_2)$coarse <- traits$Coarse_soil
meta(spec_for_plot_2)$fine <- traits$Fine_soil
meta(spec_for_plot_2)$lp <- traits$leaf_persistence

#get means
deciduous_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$lp == "d"),])
evergreen_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$lp == "e"),])

drought_yes_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$drought == "Yes"),])
drought_no_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$drought == "No"),])

coarse_yes_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$coarse == "Yes"),])
coarse_no_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$coarse == "No"),])

fine_yes_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$fine == "Yes"),])
fine_no_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$fine == "No"),])

woody_yes_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$woody == "1"),])
woody_no_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$woody == "0"),])

shade_yes_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$shade == "Tolerant"),])
shade_mid_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$shade == "Intermediate"),])
shade_no_mean2 <- mean(spec_for_plot_2[which(spec_for_plot_2$meta$shade == "Intolerant"),])

###make plots
jpeg("./output/spectra_by_lp_full100.jpg", res = 600, width = 10, height = 10, units = "in")
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$lp == "d"),], total_prob = 0.95, col = colours_trans[1], border = FALSE, ylim = c(0,0.6), main = "Spectra by leaf persistence (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association (80% quantile)"
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$lp == "e"),], total_prob = 0.95, col = colours_trans[3], border = FALSE, add = TRUE)
plot(deciduous_mean2, lwd = 0.75, lty = 1, col = colours_touse[1], add = TRUE)
plot(evergreen_mean2, lwd = 0.75, lty = 1, col = colours_touse[3], add = TRUE)
legend(1700, 0.6, legend=c("Deciduous", "Evergreen"), fill=c(colours_trans[c(1,3)]), cex=2)#lty=1,col
dev.off()

jpeg("./output/spectra_by_shade_full100.jpg", res = 600, width = 10, height = 10, units = "in")
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$shade == "Tolerant"),], total_prob = 0.95, col = colours_trans[1], border = FALSE, ylim = c(0,0.6), main = "Spectra by shade tolerance (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association (80% quantile)"
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$shade == "Intermediate"),], total_prob = 0.95, col = colours_trans[2], border = FALSE, add = TRUE)
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$shade == "Intolerant"),], total_prob = 0.95, col = colours_trans[3], border = FALSE, add = TRUE)
plot(shade_yes_mean2, lwd = 0.75, lty = 1, col = colours_touse[1], add = TRUE)
plot(shade_mid_mean2, lwd = 0.75, lty = 1, col = colours_touse[2], add = TRUE)
plot(shade_no_mean2, lwd = 0.75, lty = 1, col = colours_touse[3], add = TRUE)
legend(1700, 0.6, legend=c("Tolerant", "Intermediate", "Intolerant"), fill=c(colours_trans), cex=2)#lty=1,col
dev.off()

jpeg("./output/spectra_by_drought_full100.jpg", res = 600, width = 10, height = 10, units = "in")
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$drought == "Yes"),], total_prob = 0.95, col = colours_trans[1], border = FALSE, ylim = c(0,0.6), main = "Spectra by drought tolerance (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association (80% quantile)"
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$drought == "No"),], total_prob = 0.95, col = colours_trans[3], border = FALSE, add = TRUE)
plot(drought_yes_mean2, lwd = 0.75, lty = 1, col = colours_touse[1], add = TRUE)
plot(drought_no_mean2, lwd = 0.75, lty = 1, col = colours_touse[3], add = TRUE)
legend(1500, 0.6, legend=c("Drought tolerant", "Drought intolerant"), fill=c(colours_trans[c(1,3)]), cex=2)#lty=1,col
dev.off()

jpeg("./output/spectra_by_coarse_full100.jpg", res = 600, width = 10, height = 10, units = "in")
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$coarse == "Yes"),], total_prob = 0.95, col = colours_trans[1], border = FALSE, ylim = c(0,0.6), main = "Spectra by coarse soil adaptation (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association (80% quantile)"
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$coarse == "No"),], total_prob = 0.95, col = colours_trans[3], border = FALSE, add = TRUE)
plot(coarse_yes_mean2, lwd = 0.75, lty = 1, col = colours_touse[1], add = TRUE)
plot(coarse_no_mean2, lwd = 0.75, lty = 1, col = colours_touse[3], add = TRUE)
legend(1500, 0.6, legend=c("Adapted to coarse soil", "Not adapted to coarse soil"), fill=c(colours_trans[c(1,3)]), cex=1.5)#lty=1,col
dev.off()

jpeg("./output/spectra_by_fine_full100.jpg", res = 600, width = 10, height = 10, units = "in")
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$fine == "Yes"),], total_prob = 0.95, col = colours_trans[1], border = FALSE, ylim = c(0,0.6), main = "Spectra by fine soil adaptation (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association (80% quantile)"
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$fine == "No"),], total_prob = 0.95, col = colours_trans[3], border = FALSE, add = TRUE)
plot(fine_yes_mean2, lwd = 0.75, lty = 1, col = colours_touse[1], add = TRUE)
plot(fine_no_mean2, lwd = 0.75, lty = 1, col = colours_touse[3], add = TRUE)
legend(1500, 0.6, legend=c("Adapted to fine soil", "Not adapted to fine soil"), fill=c(colours_trans[c(1,3)]), cex=1.5)#lty=1,col
dev.off()

jpeg("./output/spectra_by_woody_full100.jpg", res = 600, width = 10, height = 10, units = "in")
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$woody == "1"),], total_prob = 0.95, col = colours_trans[1], border = FALSE, ylim = c(0,0.6), main = "Spectra by woodiness (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association (80% quantile)"
plot_quantile(spec_for_plot_2[which(spec_for_plot_2$meta$woody == "0"),], total_prob = 0.95, col = colours_trans[3], border = FALSE, add = TRUE)
plot(woody_yes_mean2, lwd = 0.75, lty = 1, col = colours_touse[1], add = TRUE)
plot(woody_no_mean2, lwd = 0.75, lty = 1, col = colours_touse[3], add = TRUE)
legend(1700, 0.6, legend=c("Woody", "Non-woody"), fill=c(colours_trans[c(1,3)]), cex=1.5)#lty=1,col
dev.off()

#overall mean
spectra_mean <- mean(spec_for_plot_2)

jpeg("./output/spectra_all_100.jpg", res = 600, width = 10, height = 10, units = "in")
plot_quantile(spec_for_plot_2, total_prob = 0.95, col = colours_trans[1], border = FALSE, ylim = c(0,0.6), main = "Spectra by woodiness (mean with 95% quantile)", xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association (80% quantile)"
plot(spectra_mean, lwd = 0.75, lty = 1, col = colours_touse[1], add = TRUE)
legend(1700, 0.6, legend=c("All spectra"), fill=c(colours_trans[c(1)]), cex=1.5)#lty=1,col
dev.off()