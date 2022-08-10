#plot spectra by ang vs conifer
library(tidyr)

#get spectra
spectra <- read.csv("./data/intermediate_objects/spectral_dataframe_ang_conifer.csv", header = TRUE, stringsAsFactors = FALSE)

colnames(spectra) <- gsub("trait.", "", colnames(spectra))

long_spectra <- spectra %>% 
  pivot_longer(cols = c(2:2002))

class(long_spectra$name) <- "integer"

#plot
jpeg("./output/ang_conifer_spectra.jpg", height = 8, width = 12, units = "in", res = 400)
ggplot(long_spectra, aes(x = name, y = value, color = lineage))+ #[which(long_spectra$X == "Abies balsamea"),]
  #geom_line(aes(x = name, y = value, group=lineage, color = lineage))+
  geom_line(alpha = 0.1, size = 0.2)+
  stat_summary(aes(group = lineage), fun = mean, geom = 'line', size=3, alpha=0.9) +
  xlab("Wavelength (nm)")+
  ylab("Reflectance (%)")+
  scale_color_manual(values = c("mediumpurple1","darkgreen"))+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  theme_bw()
dev.off()


