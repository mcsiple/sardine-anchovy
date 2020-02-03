# PACF 
# Look at PACFs of the time series to see if lags are appropriate

library(tidyverse)
source(here::here("R/DataCleanup/getMARSSstates.R"))
source(here::here("R/DataCleanup/Fill_NAs_SA.R"))
load(here::here("R/DataCleanup/allsardineanchovy_3.RData")) # dataframe alldat

region <- c("Benguela","California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
variable <- c("rec","ssb","landings")
datasource <- c("Barange")

allfishies <- alldat %>% 
  melt(id.vars=c("datasource","scientificname","stock","year","sp","region","subregion"))

domfishies <- allfishies %>%
  group_by(region,datasource,variable,sp) %>% 
  filter(value==max(value,na.rm=T)) %>%
  mutate(dom_YN = 1) %>%
  select(-value,-year) %>%
  as.data.frame()

domstocks <- left_join(allfishies,domfishies) %>% 
  filter(dom_YN==1)

domstocks %>% 
  filter(datasource=="Barange") %>% 
  ggplot(aes(x=year,y=value,colour=sp)) + 
  geom_line() + 
  facet_grid(variable~region,scales="free_y")


use.these <- domstocks %>% filter(datasource=="Barange")

v ="ssb"

regions <- c("Benguela","California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
nregions=length(regions)

for(i in 1:nregions){
    r =regions[i]
    test <- use.these%>%filter(variable==v & region==r)
    ggplot(test,aes(x=year,y=value,colour=sp)) +geom_line() 
    
    # sard <- test %>% filter(sp=="Sardine") %>% select(value)
    # anch <-  test %>% filter(sp=="Anchovy") %>% select(value)
    
    dt <- subset(use.these,variable==v & region==r)
    
    sard <- dt %>% filter(sp=="Sardine")
    anch <- dt %>% filter(sp=="Anchovy")
    maxyr <- max(dt$year)
    minyr <- min(dt$year)
    sard.mars <- FillNAs.ts(ObsMat = cbind(sard$year,sard$value), # fills in missing years with NAs
                            startyear=minyr,
                            endyear=maxyr)
    anch.mars <- FillNAs.ts(cbind(anch$year,anch$value),
                            startyear=minyr,
                            endyear=maxyr)
    regioncombo <- rbind(sard.mars[2,],anch.mars[2,])
    #select only years w data for both species:
    complete.yrs <- apply(regioncombo,MARGIN = 2,FUN = function(x) all(!is.na(x)))
    rc <- regioncombo[,complete.yrs]
    tiff(paste(r,"_SSB.tiff"),width = 8,height=8,units = "in",res = 75)
    pacf(t(rc),lag.max = 5)
    dev.off()
}



