# Baumgartner et al 1992; sediment scale deposition data for sardine and anchovy
#
# Original paper:
# Baumgartner TR, Soutar A, Ferreira-Bartrina V. 1992. Reconstruction of the history of Pacific sardine and northern anchovy populations over the past two millennia from sediments of the Santa Barbara Basin, California. CalCOFI Rep 33:24–40.

dat <- read.csv(here::here('R','DataCleanup','sardineanchovyBaumgartner.csv')) 
# units of scale dep are number/1000cm^2/yr
# These are averages across two piston cores, extracted from Fig 4 in Baumgartner et al. 1992 using the ancient software GraphClick.
# The sample timescale is 10 years.

library(dplyr)
library(ggplot2)

ggplot(dat,aes(x=year.approx,y=sdr)) + 
  geom_line() + 
  facet_wrap(~sp) + 
  theme_classic()             

# Per the paper, can convert SDR to biomass using the following regressions:
# This regression is based on 5-year averages
# sardine_B <- 0.767*SDR + 0.416
# anchovy_B <- 0.092*SDR + 0.206

dat$biomass <- NA

for(i in 1:nrow(dat)){
  if(dat$sp[i] == "sardine") {dat$biomass[i] <- 0.767*dat$sdr[i] + 0.416} 
  else {dat$biomass[i] <- 0.092*dat$sdr[i] + 0.206} # biomass units: x 10^6
}

# Plot biomass
ggplot(dat,aes(x=year.approx,y=biomass)) + 
  geom_line() + 
  facet_wrap(~sp,scales = "free_y",ncol=1) + 
  ylab("Biomass (mt * 10^6)") +
  theme_classic()         

# The covariance estimation and wavelet analysis can be applied to these data if desired.


