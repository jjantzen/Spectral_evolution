## install.packages("hsdar")
library(hsdar)

#look for different packages with prospect - cross-inversion (to try to invert the function)
#supercomputer cluster for larger dataset

##################################
## Simulating leaf spectra with PROSPECT

## PROSPECT PARAMETERS
## N: average number of air:cell wall interfaces in the mesophyll
## Cab: chlorophyll per unit area
## Car: carotenoids per unit area
## Anth: anthocyanins per unit area
## Cbrown: brown pigments per unit area
## Cw: Equivalent water thickness
## Cm: dry matter per unit area

## simulate a default spectrum using PROSPECT
orig_spec <- PROSPECT(N = 1.2, Cab = 40, Car = 12, Anth = 1.0, Cbrown = 0, 
         Cw = 0.015, Cm = 0.01, transmittance = FALSE,
         parameterList = NULL, version = "D")

## simulate a "new" spectrum for comparison
new_spec <- PROSPECT(N = 1.2, Cab = 40, Car = 12, Anth = 1.0, Cbrown = 0, 
                    Cw = 0.015, Cm = 0.015, transmittance = FALSE,
                    parameterList = NULL, version = "D")

new_spec2 <- PROSPECT(N = 1.2, Cab = 40, Car = 12, Anth = 1.0, Cbrown = 0, 
                   Cw = 0.015, Cm = 0.02, transmittance = FALSE,
                   parameterList = NULL, version = "D")


## simulate a default spectrum using PROSPECT - vary carotenoids
orig_spec <- PROSPECT(N = 1.2, Cab = 40, Car = 12, Anth = 1.0, Cbrown = 0, 
                      Cw = 0.015, Cm = 0.01, transmittance = FALSE,
                      parameterList = NULL, version = "D")

## simulate a "new" spectrum for comparison
new_spec <- PROSPECT(N = 1.2, Cab = 140, Car = 12, Anth = 1.0, Cbrown = 0, 
                     Cw = 0.015, Cm = 0.01, transmittance = FALSE,
                     parameterList = NULL, version = "D")

new_spec2 <- PROSPECT(N = 1.2, Cab = 400, Car = 12, Anth = 1.0, Cbrown = 0, 
                      Cw = 0.015, Cm = 0.01, transmittance = FALSE,
                      parameterList = NULL, version = "D")

## simulate a "new" spectrum for comparison
new_spec3 <- PROSPECT(N = 1.2, Cab = 40, Car = 100, Anth = 1.0, Cbrown = 0, 
                     Cw = 0.015, Cm = 0.01, transmittance = FALSE,
                     parameterList = NULL, version = "D")

new_spec4 <- PROSPECT(N = 1.2, Cab = 40, Car = 50, Anth = 1.0, Cbrown = 0, 
                      Cw = 0.015, Cm = 0.01, transmittance = FALSE,
                      parameterList = NULL, version = "D")

## plot all three spectra
plot(orig_spec,ylim=c(0,0.6),xlim=c(400,2500))
points(400:2500,new_spec@spectra@spectra_ma,type="l",col="red")
points(400:2500,new_spec2@spectra@spectra_ma,type="l",col="blue")
points(400:2500,new_spec3@spectra@spectra_ma,type="l",col="green")
points(400:2500,new_spec4@spectra@spectra_ma,type="l",col="purple")

## show PRI wavelengths
abline(v=531,lty="dashed")
abline(v=570,lty="dashed")

######################################
## Simulating canopy spectra with PROSAIL

## PROSAIL PARAMETERS
## psoil: is soil reflectance Lambertian?
## LAI: Leaf area index, or leaf area per ground area
## TypeLidf, lidfa, lidfb: parameters related to leaf angle distribution
## hspot, tts, tto, psi: parameters related to hotspot, solar and observer angles

canopy_spec1<-PROSAIL(N = 1.5, Cab = 40, Car = 8, Cbrown = 0.0,
                      Cw = 0.01, Cm = 0.009, psoil = 0, LAI = 0.5, 
                      TypeLidf = 1, lidfa = -0.35, lidfb = -0.15,
                      hspot = 0.01, tts = 30, tto = 10, psi = 0,
                      parameterList = NULL, rsoil = NULL)

canopy_spec2<-PROSAIL(N = 1.5, Cab = 40, Car = 8, Cbrown = 0.0,
                      Cw = 0.01, Cm = 0.009, psoil = 0, LAI = 1.5, 
                      TypeLidf = 1, lidfa = -0.35, lidfb = -0.15,
                      hspot = 0.01, tts = 30, tto = 10, psi = 0,
                      parameterList = NULL, rsoil = NULL)

canopy_spec3<-PROSAIL(N = 1.5, Cab = 40, Car = 8, Cbrown = 0.0,
                      Cw = 0.01, Cm = 0.009, psoil = 0, LAI = 2.5, 
                      TypeLidf = 1, lidfa = -0.35, lidfb = -0.15,
                      hspot = 0.01, tts = 30, tto = 10, psi = 0,
                      parameterList = NULL, rsoil = NULL)

canopy_spec4<-PROSAIL(N = 1.5, Cab = 20, Car = 8, Cbrown = 0.0,
                      Cw = 0.01, Cm = 0.009, psoil = 0, LAI = 1.5, 
                      TypeLidf = 1, lidfa = -0.35, lidfb = -0.15,
                      hspot = 0.01, tts = 30, tto = 10, psi = 0,
                      parameterList = NULL, rsoil = NULL)

canopy_spec5<-PROSAIL(N = 1.5, Cab = 60, Car = 8, Cbrown = 0.0,
                      Cw = 0.01, Cm = 0.009, psoil = 0, LAI = 1.5, 
                      TypeLidf = 1, lidfa = -0.35, lidfb = -0.15,
                      hspot = 0.01, tts = 30, tto = 10, psi = 0,
                      parameterList = NULL, rsoil = NULL)

## varying LAI
plot(canopy_spec1,ylim=c(0,0.6))
points(400:2500,canopy_spec2@spectra@spectra_ma,type="l",col="red")
points(400:2500,canopy_spec3@spectra@spectra_ma,type="l",col="blue")

## show NDVI wavelengths
abline(v=670,lty="dashed")
abline(v=800,lty="dashed")

## varying chlorophyll
plot(canopy_spec4,ylim=c(0,0.6))
points(400:2500,canopy_spec2@spectra@spectra_ma,type="l",col="red")
points(400:2500,canopy_spec5@spectra@spectra_ma,type="l",col="blue")

## show NDVI wavelengths
abline(v=670,lty="dashed")
abline(v=800,lty="dashed")
