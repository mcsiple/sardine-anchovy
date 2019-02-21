# Combine all data from RAM, FAO, Barange et al.
# The resulting dataset is "allsardineanchovy.RData" which is used in all subsequent analyses
# NOTE: One important thing this code does is turn the "NA" stocks in the FAO database into fake stocks.
# For example, there is a time series of sardine in land area "Europe" and one in land area "-"

datadir <- "C:/Users/siplem/Dropbox/Chapter3-SardineAnchovy/Datasets"

# For unix OS:
# datadir <- "~/Dropbox/Chapter3-SardineAnchovy/Datasets/"

load(file.path(datadir,"RAM/RAM_2.RData"))      #RAM
load(file.path(datadir,"FAO/FAO.RData"))      #FAO
load(file.path(datadir,"Barange/BARANGE_ALL.RData"))   #ALLDAT
barange <- alldat
#load(file.path(datadir,"Barange/Barange_mystocks.RData"))   #barange_noNAs - DEPRACATED

library(reshape2)
library(tidyverse)
library(taxize) # to get common names from scientific names

# Combine all data so that it has all the same columns --------------------

# stock  year  ssb   rec   landings    fishing.mortality   sp    region    subregion
# Following format of barange_noNAs

barange.new <- data.frame(datasource = rep("Barange",times=nrow(barange)),
                          scientificname= rep(NA, times = nrow(barange)),
                          stock=barange$stock,
                          year=barange$year,
                          ssb=barange$SSB,
                          rec=barange$Rec,
                          landings=barange$landings,
                          fishing.mortality=barange$fishing.mortality,
                          sp=barange$sp,
                          region=barange$region,
                          subregion=barange$subregion) %>% filter(!is.na(year)) # had to remove entries where year is NA

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



FAO.new <- data.frame(datasource="FAO",
                      scientificname=FAO$Scientific.name,
                      stock=paste(FAO$Ocean.Area), # This is just so that there is a stock name - needed for future comparisons
                      year=FAO$year,
                      ssb=NA,
                      rec=NA,
                      landings=FAO$landings,
                      fishing.mortality=rep(NA,times=nrow(FAO)),
                      sp=FAO$sp,
                      region=FAO$region)

# We need better stock names for FAO stocks. Going to combine common name and Ocean.area to make a stock name.
stocktable <- data.frame(scientificname = unique(FAO.new$scientificname),
                         commonname = c("European anchovy",
                                        "European pilchard",
                                        "Sardinella spp",
                                        " Engraulidae",
                                        "Round sardinella",
                                        "Southern African anchovy",
                                        "Southern African pilchard",
                                        "California pilchard",
                                        "California anchovy",
                                        "Pacific anchoveta",
                                        "Japanese anchovy",
                                        "Japanese pilchard",
                                        "Japanese sardinella",
                                        "Slender raindbow sardine",
                                        "Stolephorus spp",
                                        "Peruvian anchoveta",
                                        "Longnose anchovy",
                                        "South American pilchard")  ) # These common names are not from taxize, but should be very close. 

FAO.new <- inner_join(FAO.new,stocktable) %>% mutate(stock = paste(commonname,stock,sep=" - "))                       
                         
#FAO.new$year <- as.numeric(levels(FAO.new$year))[FAO.new$year]


# Do some summary stats to double check data ------------------------------
barange.new %>% 
  group_by(sp,region,stock) %>% 
  summarize(max.ssb = max(ssb,na.rm=T), 
            max.rec = max(rec,na.rm=T),
            max.landings=max(landings,na.rm=T),
            max.f=max(fishing.mortality,na.rm=T)) %>%
            as.data.frame()



# Output just one species for ADMB
# European pilchard, 34 years of data
# Epilchard <- subset(RAM.new,scientificname=="Sardina pilchardus")
# Ssagax <- subset(RAM.new,scientificname=="Sardinops sagax" & region =="California")
# write.csv(x = Ssagax,file = "Sardinops_sagax.csv")

colnames(barange.new)
alldat <- plyr::rbind.fill(list(barange.new,RAM.new,FAO.new)) %>% select(datasource,scientificname,stock,year,ssb,rec,landings,fishing.mortality,sp,region,subregion,totalcatch)

  melt(alldat,
          id.vars = c("datasource","scientificname","stock","year","sp","region","subregion")) %>% 
  subset(variable !="rec") %>%
  ggplot(aes(x=year,y=value,linetype=datasource,colour=stock)) +
  geom_line() +
  facet_grid(variable~region,scales="free_y") +
  theme_classic() +
  theme(legend.position = "none")

  
# some time series have subregion info, so plotting by region needs to be preceded by summing by 
  
# Save the huge dataset!
save(alldat,file="~/Dropbox/Chapter3-SardineAnchovy/Datasets/allsardineanchovy_3.RData")



# Fill NAs with long-term means -- see the getMARSSstates.R file!



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

