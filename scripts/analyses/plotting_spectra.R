#plotting real spectra
library(spectrolab)
library(ggplot2)
library(phytools)

#improve the colour choices


#plot actual spectra by myc association

#read spectra
data_spectra <- readRDS("./data/for_analysis/myc_data_list_92sp_binary_for_analysis.rds")

#convert to spectra object
spec_for_plot <- as_spectra(data_spectra$spectra)

#get the metadata in right format
meta(spec_for_plot)$myc <- data_spectra$myc

#get means
AM_mean <- mean(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),])

EM_mean <- mean(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),])

#colours_touse <- c("#008837", "#E69F00", "#7b3294")
colours_touse <- c("darkorange", "blue", "#7b3294")

#colours_trans <- c(make.transparent("#008837", 0.15), make.transparent("#E69F00", 0.15), make.transparent("#7b3294", 0.2))
colours_trans <- c(make.transparent("orange", 0.15), make.transparent("blue", 0.1), make.transparent("#7b3294", 0.2))

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

jpeg("./output/spectra_by_myc_orange.jpg", res = 600, width = 6, height = 6, units = "in")
#pdf("./output/spectra_by_myc_orange.pdf", width = 6, height = 6)
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], ylim = c(0,0.6), total_prob = 0.95, col = colours_trans[1], cex.lab = 1, cex.main = 1.2, cex.axis = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association\n(mean with 95% quantile)", 
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], total_prob = 0.95, col = colours_trans[2], border = FALSE, add = TRUE)
plot(AM_mean, lwd = 1, lty = 1.5, col = colours_touse[1], add = TRUE)
plot(EM_mean, lwd = 1, lty = 1.5, col = colours_touse[2], add = TRUE)
legend(2000, 0.6, legend=c("AM", "EM"), fill=c(colours_trans[c(1,2)]), border = c(colours_touse[c(1,2)]), cex=1)#lty=1,col
dev.off()


#split into ang conifer 

#get the metadata in right format
ang_trees <- readRDS("./data/for_analysis/ang_only_trees_for_myc.rds")

all_trees <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

ang_species <- ang_trees[[1]]$tip.label

all_species <- all_trees[[1]]$tip.label

#get subset of species for myc
data_spectra$clade[which(!(data_spectra$species %in% ang_species))] <- "Con"
data_spectra$clade[which(data_spectra$species %in% ang_species)] <- "Ang"

meta(spec_for_plot)$clade <- data_spectra$clade
meta(spec_for_plot)$species <- data_spectra$species

str(spec_for_plot)

#get means
Ang_mean <- mean(spec_for_plot[which(spec_for_plot$meta$clade == "Ang"),])

Con_mean <- mean(spec_for_plot[which(spec_for_plot$meta$clade == "Con"),])

Ang_mean_AM <- mean(spec_for_plot[which(spec_for_plot$meta$clade == "Ang" & spec_for_plot$meta$myc == "AM"),])
Ang_mean_EM <- mean(spec_for_plot[which(spec_for_plot$meta$clade == "Ang" & spec_for_plot$meta$myc == "EM"),])

Con_mean_AM <- mean(spec_for_plot[which(spec_for_plot$meta$clade == "Con" & spec_for_plot$meta$myc == "AM"),])
Con_mean_EM <- mean(spec_for_plot[which(spec_for_plot$meta$clade == "Con" & spec_for_plot$meta$myc == "EM"),])



colours_trans <- c(make.transparent("orange", 0.15), make.transparent("blue", 0.1), make.transparent("#7b3294", 0.2))

jpeg("./output/spectra_by_myc_3panels_individual_means.jpg", res = 600, width = 6, height = 6, units = "in")
#pdf("./output/spectra_by_myc_3panels.pdf", width = 6, height = 6)
#plot.new() ## clean up device
#par(mfrow = c(2, 2))
par(cex=0.7, mai=c(0.6,0.7,0.1,0.1))
par(fig=c(0.01,0.99,0.51,0.99))
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], ylim = c(0,0.6), total_prob = 0.95, col = colours_trans[1], cex.lab = 1, cex.main = 1.2, cex.axis = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association\n(mean with 95% quantile)", 
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], total_prob = 0.95, col = colours_trans[2], border = FALSE, add = TRUE)
#plot(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[1], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
#plot(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[2], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
plot(AM_mean, lwd = 1, lty = 1.5, col = colours_touse[1], add = TRUE)
plot(EM_mean, lwd = 1, lty = 1.5, col = colours_touse[2], add = TRUE)
legend(2000, 0.6, legend=c("AM", "EM"), fill=c(colours_trans[c(1,2)]), border = c(colours_touse[c(1,2)]), cex=1)#lty=1,col
#p <- recordPlot()


#plot.new() ## clean up device
#par(mfrow = c(1, 1))
par(fig=c(0.01,0.49,0.01,0.49), new = TRUE)
#plot_quantile(spec_for_plot[which(spec_for_plot$meta$clade == "Ang"),], ylim = c(0,0.6), total_prob = 0.95, col = "mediumpurple1", cex.lab = 1, cex.main = 1.2, cex.axis = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association\n(mean with 95% quantile)", 
plot(spec_for_plot[which(spec_for_plot$meta$clade == "Ang"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = "mediumpurple1", cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
plot(Ang_mean, lwd = 1, lty = 1.5, col = "purple4", add = TRUE)
legend(1500, 0.6, legend=c("Angiosperms"), fill=c("mediumpurple1"), border = c("purple4"), cex=0.7)#lty=1,col
#q1 <- recordPlot()


#plot.new() ## clean up device
par(fig=c(0.51,0.99,0.01,0.49), new = TRUE)
#plot_quantile(spec_for_plot[which(spec_for_plot$meta$clade == "Con"),], ylim = c(0,0.6), total_prob = 0.95, col = "lightgreen", cex.lab = 1, cex.main = 1.2, cex.axis = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association\n(mean with 95% quantile)", 
plot(spec_for_plot[which(spec_for_plot$meta$clade == "Con"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = "seagreen3", cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
plot(Con_mean, lwd = 1, lty = 1.5, col = "darkgreen", add = TRUE)
legend(1500, 0.6, legend=c("Conifers and fern"), fill=c("seagreen3"), border = c("darkgreen"), cex=0.7)#lty=1,col
#q2 <- recordPlot()

#ggarrange(p, q1, q2)
dev.off()


colours_trans <- c(make.transparent("orange", 0.4), make.transparent("blue", 0.4), make.transparent("#7b3294", 0.4))

jpeg("./output/spectra_by_myc_3panels_individual_spectra_ang_conifer.jpg", res = 600, width = 6, height = 6, units = "in")
#pdf("./output/spectra_by_myc_3panels.pdf", width = 6, height = 6)
#plot.new() ## clean up device
#par(mfrow = c(2, 2))
par(cex=0.7, mai=c(0.6,0.7,0.1,0.1))
par(fig=c(0.01,0.49,0.51,0.99))
#plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], ylim = c(0,0.6), total_prob = 0.95, col = colours_trans[1], cex.lab = 1, cex.main = 1.2, cex.axis = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association\n(mean with 95% quantile)", 
#plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], total_prob = 0.95, col = colours_trans[2], border = FALSE, add = TRUE)
plot(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[1], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
#plot(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[2], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
plot(AM_mean, lwd = 1, lty = 1.5, col = colours_touse[1], add = TRUE)
#plot(EM_mean, lwd = 1, lty = 1.5, col = colours_touse[2], add = TRUE)
legend(1500, 0.6, legend=c("AM"), fill=c(colours_trans[c(1)]), border = c(colours_touse[c(1)]), cex=1)#lty=1,col
#p <- recordPlot()

par(fig=c(0.51,0.99,0.51,0.99), new = TRUE)
#plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], ylim = c(0,0.6), total_prob = 0.95, col = colours_trans[1], cex.lab = 1, cex.main = 1.2, cex.axis = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association\n(mean with 95% quantile)", 
#plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], total_prob = 0.95, col = colours_trans[2], border = FALSE, add = TRUE)
#plot(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[1], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
plot(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[2], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
#plot(AM_mean, lwd = 1, lty = 1.5, col = colours_touse[1], add = TRUE)
plot(EM_mean, lwd = 1, lty = 1.5, col = colours_touse[2], add = TRUE)
legend(1500, 0.6, legend=c("EM"), fill=c(colours_trans[c(2)]), border = c(colours_touse[c(2)]), cex=1)#lty=1,col
#p <- recordPlot()

#plot.new() ## clean up device
#par(mfrow = c(1, 1))
par(fig=c(0.01,0.49,0.01,0.49), new = TRUE)
#plot_quantile(spec_for_plot[which(spec_for_plot$meta$clade == "Ang"),], ylim = c(0,0.6), total_prob = 0.95, col = "mediumpurple1", cex.lab = 1, cex.main = 1.2, cex.axis = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association\n(mean with 95% quantile)", 
plot(spec_for_plot[which(spec_for_plot$meta$clade == "Ang"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = "mediumpurple1", cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance", border = FALSE)
plot(Ang_mean, lwd = 1, lty = 1.5, col = "mediumpurple4", add = TRUE)
legend(1500, 0.6, legend=c("Angiosperm"), fill=c("mediumpurple1"), border = c("mediumpurple4"), cex=0.7)#lty=1,col
#q1 <- recordPlot()


#plot.new() ## clean up device
par(fig=c(0.51,0.99,0.01,0.49), new = TRUE)
#plot_quantile(spec_for_plot[which(spec_for_plot$meta$clade == "Con"),], ylim = c(0,0.6), total_prob = 0.95, col = "lightgreen", cex.lab = 1, cex.main = 1.2, cex.axis = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association\n(mean with 95% quantile)", 
plot(spec_for_plot[which(spec_for_plot$meta$clade == "Con"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = "seagreen3", cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
plot(Con_mean, lwd = 1, lty = 1.5, col = "darkgreen", add = TRUE)
legend(1500, 0.6, legend=c("Conifers and fern"), fill=c("seagreen3"), border = c("darkgreen"), cex=0.7)#lty=1,col
#q2 <- recordPlot()

#ggarrange(p, q1, q2)
dev.off()


#final figure for manuscript
colours_trans <- c(make.transparent("orange", 0.15), make.transparent("blue", 0.1), make.transparent("#7b3294", 0.2))

jpeg("./output/spectra_by_myc_3panels_ang_con_split_by_myc.jpg", res = 600, width = 6, height = 6, units = "in")
#pdf("./output/spectra_by_myc_3panels_ang_con_split_by_myc.pdf", width = 6, height = 6)
#plot.new() ## clean up device
#par(mfrow = c(2, 2))
par(cex=0.7, mai=c(0.6,0.7,0.1,0.1), xpd=TRUE)
par(fig=c(0.01,0.99,0.51,0.99))

plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], ylim = c(0,0.6), total_prob = 0.95, col = colours_trans[1], cex.lab = 1, cex.main = 1.2, cex.axis = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association\n(mean with 95% quantile)", 
plot_quantile(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], total_prob = 0.95, col = colours_trans[2], border = FALSE, add = TRUE)
#plot(spec_for_plot[which(spec_for_plot$meta$myc == "AM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[1], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
#plot(spec_for_plot[which(spec_for_plot$meta$myc == "EM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[2], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
plot(AM_mean, lwd = 1, lty = 1.5, col = colours_touse[1], add = TRUE)
plot(EM_mean, lwd = 1, lty = 1.5, col = colours_touse[2], add = TRUE)
legend(2000, 0.6, legend=c("AM", "EM"), fill=c(colours_trans[c(1,2)]), border = c(colours_touse[c(1,2)]), cex=1)#lty=1,col
text(115,0.63, labels=c("a"), cex=2)#lty=1,col
#mtext('a', side=2, line=10, at=0)
#title(main="a",line=300, font=2) # adj=1, 

#plot.new() ## clean up device
#par(mfrow = c(1, 1))
par(fig=c(0.01,0.49,0.01,0.49), new = TRUE)
#plot_quantile(spec_for_plot[which(spec_for_plot$meta$clade == "Ang"),], ylim = c(0,0.6), total_prob = 0.95, col = "mediumpurple1", cex.lab = 1, cex.main = 1.2, cex.axis = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association\n(mean with 95% quantile)", 
plot_quantile(spec_for_plot[which(spec_for_plot$meta$clade == "Ang" & spec_for_plot$meta$myc == "AM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[1], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance", border = FALSE)
plot_quantile(spec_for_plot[which(spec_for_plot$meta$clade == "Ang" & spec_for_plot$meta$myc == "EM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[2], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance", add = TRUE, border = FALSE)
plot(Ang_mean_AM, lwd = 1, lty = 1.5, col = colours_touse[1], add = TRUE)
plot(Ang_mean_EM, lwd = 1, lty = 1.5, col = colours_touse[2], add = TRUE)
legend(1500, 0.6, legend=c("Angiosperm AM", "Angiosperm EM"), fill=c(colours_trans[1], colours_trans[2]), border = c(colours_touse[1], colours_touse[2]), cex=0.7)#lty=1,col
#q1 <- recordPlot()
text(-180, 0.63, labels=c("b"), cex=2)#lty=1,col

#plot.new() ## clean up device
par(fig=c(0.51,0.99,0.01,0.49), new = TRUE)
#plot_quantile(spec_for_plot[which(spec_for_plot$meta$clade == "Con"),], ylim = c(0,0.6), total_prob = 0.95, col = "lightgreen", cex.lab = 1, cex.main = 1.2, cex.axis = 1, border = FALSE, xlab = "Wavelength (nm)", ylab = "Reflectance") #main = "Spectra by mycorrhizal association\n(mean with 95% quantile)", 
#plot(spec_for_plot[which(spec_for_plot$meta$clade == "Con"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = "seagreen3", cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance")
#plot(Con_mean, lwd = 1, lty = 1.5, col = "darkgreen", add = TRUE)
plot_quantile(spec_for_plot[which(spec_for_plot$meta$clade == "Con" & spec_for_plot$meta$myc == "AM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[1], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance", border = FALSE)
plot_quantile(spec_for_plot[which(spec_for_plot$meta$clade == "Con" & spec_for_plot$meta$myc == "EM"),], lwd = 0.25, ylim = c(0,0.6), lty = 1, col = colours_trans[2], cex.lab = 1, cex.main = 1.2, cex.axis = 1, xlab = "Wavelength (nm)", ylab = "Reflectance", add = TRUE, border = FALSE)
plot(Con_mean_AM, lwd = 1, lty = 1.5, col = colours_touse[1], add = TRUE)
plot(Con_mean_EM, lwd = 1, lty = 1.5, col = colours_touse[2], add = TRUE)
legend(1300, 0.6, legend=c("Conifer and fern AM", "Conifer and fern EM"), fill=c(colours_trans[1], colours_trans[2]), border = c(colours_touse[1], colours_touse[2]), cex=0.7)#lty=1,col
#q2 <- recordPlot()
text(-180, 0.63, labels=c("c"), cex=2)#lty=1,col
#ggarrange(p, q1, q2)
dev.off()




#convert to spectra object
spec_for_plot <- as_spectra(data_spectra$spectra)

#get the metadata in right format
meta(spec_for_plot)$myc <- data_spectra$myc








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