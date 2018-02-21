# Combine all data from RAM, FAO, Barange et al.
# The resulting dataset is "allsardineanchovy.RData" which is used in all subsequent analyses
# NOTE: One important thing this code does is turn the "NA" stocks in the FAO database into fake stocks.
# For example, there is a time series of sardine in land area "Europe" and one in land area "-"

basedir <- "~/Dropbox/Chapter3-SardineAnchovy"

load(file.path(basedir,"Datasets/RAM/RAM.RData"))      #RAM
load(file.path(basedir,"Datasets/FAO/FAO.RData"))      #FAO
load(file.path(basedir,"Datasets/Barange/Barange_mystocks.RData"))   #barange_noNAs
b <- which(colnames(barange_noNAs)=='b')
colnames(barange_noNAs)[b] <- 'sp'

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)


# Combine all data so that it has all the same columns --------------------

# Stock  year  ssb   rec   landings    fishing.mortality   sp    region    subregion
# Following format of barange_noNAs

barange.new <- data.frame(datasource=rep("Barange",times=nrow(barange_noNAs)),
                          scientificname=rep(NA,times=nrow(barange_noNAs)),
                          barange_noNAs)

RAM.new <- data.frame(datasource=rep("RAM",times=nrow(RAM)),
                      scientificname=RAM$scientificname,
                      stock=RAM$stocklong,
                      year=RAM$Years,
                      ssb=RAM$SSB,
                      rec=RAM$R,
                      totalcatch=RAM$TC,
                      fishing.mortality=RAM$F,
                      sp=RAM$sp,
                      region=RAM$region,
                      subregion=rep(NA,times=nrow(RAM)))

FAO[grep("Southern African pilchard",FAO$Species), 2] <- "Africa"
FAO.new <- data.frame(datasource="FAO",
                      scientificname=FAO$Scientific.name,
                      stock=paste(FAO$Ocean.Area,FAO$Land.Area), # This is just so that there is a stock name - needed for future comparisons
                      year=FAO$year,
                      ssb=NA,
                      rec=NA,
                      landings=FAO$landings,
                      fishing.mortality=rep(NA,times=nrow(FAO)),
                      sp=FAO$sp,
                      region=FAO$region,subregion=FAO$Land.Area)

FAO.new$year <- as.numeric(levels(FAO.new$year))[FAO.new$year]

# Output just one species for ADMB
# European pilchard, 34 years of data
# Epilchard <- subset(RAM.new,scientificname=="Sardina pilchardus")
# Ssagax <- subset(RAM.new,scientificname=="Sardinops sagax" & region =="California")
# write.csv(x = Ssagax,file = "Sardinops_sagax.csv")

colnames(barange.new)
alldat <- rbind.fill(list(barange.new,RAM.new,FAO.new))

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

