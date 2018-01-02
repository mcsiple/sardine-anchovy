# Baumgartner data - sardine-anchovy
# 

basedir <- "/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Datasets/Baumgartner"
dat <- read.csv(file.path(basedir,"sardineanchovyBaumgartner.csv")) # units of scale dep are number/1000cm^2/yr
# These are averages across two piston cores (extracted from Fig 4 in Baumgartner et al. 1992)
library(plyr)
library(dplyr)
library(ggplot2)

ggplot(dat,aes(x=year.approx,y=sdr)) + 
  geom_line() + 
  facet_wrap(~sp) + 
  theme_classic()             

# Per the paper, can convert SDR to biomass using the following regressions
# This regression is based on 5-year averages
# sardine_B <- 0.767*SDR + 0.416
# anchovy_B <- 0.092*SDR + 0.206
dat$biomass <- NA
for(i in 1:nrow(dat)){
  if(dat$sp[i] == "sardine") {dat$biomass[i] <- 0.767*dat$sdr[i] + 0.416} 
  else {dat$biomass[i] <- 0.092*dat$sdr[i] + 0.206} # biomass units: x 10^6
}

# Plot biomass over time
ggplot(dat,aes(x=year.approx,y=biomass)) + 
  geom_line() + 
  facet_wrap(~sp,scales = "free_y",ncol=1) + 
  ylab("Biomass (mt * 10^6)") +
  theme_classic()         

sard = subset(dat, sp=="sardine")
1e6 * quantile(sard$biomass,probs = c(0.25,0.75))
# SAMPLE TIMESCALE: 10 YEARS! - They do spectral analysis looking at periods only longer than 50 years
