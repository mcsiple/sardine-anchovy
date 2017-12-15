# Baumgartner data - sardine-anchovy
basedir <- "/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Datasets/Baumgartner"
dat <- read.csv(file.path(basedir,"sardineanchovyBaumgartner.csv")) # units of scale dep are number/1000cm^2/yr
library(plyr)
library(dplyr)
library(ggplot2)

ggplot(dat,aes(x=year.approx,y=sdr)) + geom_line() + facet_wrap(~sp) + theme_classic()             

