# Combine all data from RAM, FAO, Barange et al.
# The resulting dataset is "allsardineanchovy.RData" which is used in all subsequent analyses.
setwd("~/Dropbox/Chapter3-SardineAnchovy/R files")

load("~/Dropbox/Chapter3-SardineAnchovy/All Datasets/RAM/RAM.RData")      #RAM
load("~/Dropbox/Chapter3-SardineAnchovy/All Datasets/FAO/FAO.RData")      #FAO
load("~/Dropbox/Chapter3-SardineAnchovy/All Datasets/Barange/Barange_mystocks.RData")   #barange_noNAs
b <- which(colnames(barange_noNAs)=='b')
colnames(barange_noNAs)[b] <- 'sp'
library(plyr)
library(dplyr)

library(ggplot2)
library(reshape2)


# Combine all data so that it has all the same columns --------------------

# Stock  year  ssb   rec   landings    fishing.mortality   sp    region    subregion
# Following format of barange_noNAs

barange.new <- data.frame(datasource=rep("Barange",times=nrow(barange_noNAs)),scientificname=rep(NA,times=nrow(barange_noNAs)),barange_noNAs)

RAM.new <- data.frame(datasource=rep("RAM",times=nrow(RAM)),scientificname=RAM$scientificname,stock=RAM$stocklong,year=RAM$Years,ssb=RAM$SSB,rec=RAM$R,totalcatch=RAM$TC,fishing.mortality=RAM$F,sp=RAM$sp,region=RAM$region,subregion=rep(NA,times=nrow(RAM)))

FAO[grep("Southern African pilchard",FAO$Species), 2] <- "Africa"
FAO.new <- data.frame(datasource=rep("FAO",times=nrow(FAO)),scientificname=FAO$Scientific.name,stock=rep(NA,times=nrow(FAO)),year=FAO$year,ssb=rep(NA,times=nrow(FAO)),rec=rep(NA,times=nrow(FAO)),landings=FAO$landings,fishing.mortality=rep(NA,times=nrow(FAO)),sp=FAO$sp,region=FAO$region,subregion=FAO$Land.Area)

# Output just one species for ADMB
# European pilchard, 34 years of data
# Epilchard <- subset(RAM.new,scientificname=="Sardina pilchardus")
# Ssagax <- subset(RAM.new,scientificname=="Sardinops sagax" & region =="California")
# write.csv(x = Ssagax,file = "Sardinops_sagax.csv")
colnames(barange.new)
alldat <- rbind.fill(list(barange.new,RAM.new,FAO.new))
alldat$year <- as.numeric(alldat$year)

# Save the huge dataset!
save(alldat,file="allsardineanchovy.RData")



########################################################################
#Side note: Pacific sardine in RAM does not look like the one in Barange et al. -- is it possible that they are from different sources? Different maxima, years covered; same overall pattern.
PSRAM <- subset(alldat,sp=="Sardine" & region=="California" & datasource=="RAM")
PSBarange <- subset(alldat,sp=="Sardine" & region=="California" & datasource=="Barange")

yearrange <- range(c(PSRAM$year,PSBarange$year))
plot(PSRAM$year,PSRAM$rec/10^6,type='l',xlim=yearrange)
lines(PSBarange$year,PSBarange$rec/10^6,col="blue")

167889000/1e6

stdrec <- (Ssagax$rec - mean(Ssagax$rec)) / mean(Ssagax$rec)
plot(stdrec)
abline(h=0,col="blue")

