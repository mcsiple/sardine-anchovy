# Comparability between time series

# Data checks
all.medians <- alldat %>%
  group_by(datasource,region,sp,stock) %>%
  summarize(median.ssb = median(ssb,na.rm=T),
            median.landings = median(landings,na.rm=T),
            median.tc = median(totalcatch,na.rm=T)) 

dom.by.ssb <- all.medians %>%
  mutate(is.dominant = ifelse(median.ssb == max(median.ssb,na.rm = T),1,0)) %>%
  filter(is.dominant == 1) %>%
  as.data.frame() 

dom.by.landings <- all.medians %>%
  mutate(is.dominant = ifelse(median.landings == max(median.landings,na.rm = T),1,0)) %>%
  filter(is.dominant == 1) %>%
  as.data.frame() 

dom.by.landings %>% 
  reshape2::dcast(datasource+region~sp,value.var = 'median.landings') %>%
  mutate(logAS = log10(Anchovy/Sardine))

dom.by.tc <- all.medians %>%
  mutate(is.dominant = ifelse(median.tc == max(median.tc,na.rm = T),1,0)) %>%
  filter(is.dominant == 1) %>%
  as.data.frame() 

dom.by.ssb %>% 
  reshape2::dcast(datasource+region~sp,value.var = 'median.ssb') %>%
  mutate(logAS = log10(Anchovy/Sardine))

all.maxes <- alldat %>%
  group_by(datasource,region,sp,stock) %>%
  summarize(max.ssb = max(ssb,na.rm=T),
            max.landings = max(landings,na.rm=T),
            max.tc = max(totalcatch,na.rm=T))


# Code used to get comparability between dominant stocks ------------------
source(here::here("R/DominantTS/Corr_Fig.R")) # corr.fig()

# Test function
corr.fig(region_or_subregion="Benguela",variable="landings",scale = "Region",data_source="Barange")
corr.fig(region_or_subregion="California",variable="ssb",scale = "Region",data_source="RAM")
corr.fig(region_or_subregion="California",variable="landings",scale = "Region",data_source="FAO")

# Summarize peak variables for all the regions.
regions <- c("Benguela","California","Humboldt","Kuroshio-Oyashio","NE Atlantic")
variables <- c("ssb","rec","landings","fishing.mortality")
nregions <- length(regions)
nvariables <- length(variables)
megatable <- vector()

for (r in 1:nregions){
  for (v in 1: nvariables){
    xx <- try(corr.fig(region_or_subregion=regions[r],
                       variable=variables[v],
                       scale = "Region",
                       data_source="FAO")) #must rotate sources
    mt <-  xx$median.table # previously: max.table
    mt$variable=rep(variables[v],times=nrow(mt))
    mt$region=rep(regions[r],times=nrow(mt))  
    mt$domsard <- ifelse(length(xx$Dom_sard_LTmax)==0,NA,xx$Dom_sard_LTmax)
    mt$domanch <- ifelse(length(xx$Dom_anch_LTmax)==0,NA,xx$Dom_anch_LTmax)
    megatable <- bind_rows(megatable,mt)
  }
}

# Clunky- when you change source above, save each datasource like so:
RAM.summary <- megatable[-1,]
RAM.summary$datasource <- "RAM"

barange.summary <- megatable[-1,]
barange.summary$datasource <- "Barange"

FAO.summary <-megatable[-1,]
FAO.summary$datasource <- "FAO"

Both <- rbind(RAM.summary,barange.summary,FAO.summary)

#save(Both,file = "Replacement_RAM_Barange_MEDIAN2.Rdata") # For the new function that defines "dominant" species by which one has the highest median biomass
