#plot simulated spectra

str(bm_data_sim)




i <- 1
#plot
jpeg("./output/sim_spectra_full_bm_sigma0.5.jpg", height = 8, width = 8, units = "in", res = 400)

for (i in 1:1){ #length(bm_data_sim)
  plot(NULL, ylim=c(0,0.6),xlim=c(400,2500), type = "l", col = "blue", ylab = "Reflectance", xlab = "Wavelength (nm)")
  for (j in 1:length(bm_data_sim_redo[[i]]$trait[,1])){
    #points(400:2500,bm_data_sim_redo[[i]]$trait[j,],type="l",col="blue")
    #points(400:2500,ou_data_sim[[i]]$trait[j,],type="l",col="red")
    points(400:2500,bm_data_sim_sigma0.5[[i]]$trait[j,],type="l",col="black")
  }
}

dev.off()

str(bm_data_sim[[i]]$trait)
points(400:2500,ou_data_sim[[i]]$trait[1,],type="l",col="green")


for (i in 1:1){ #length(ou_data_sim)
  for (j in 1:length(ou_data_sim[[i]]$trait[,1]))
    if (i == 1 & j == 1) {
      plot(ou_data_sim[[i]]$trait[j,],ylim=c(0,0.6),xlim=c(400,2500), type = "l", col = "red")
    }  else {
      points(400:2500,ou_data_sim[[i]]$trait[j,],type="l",col="red")
    }
}
