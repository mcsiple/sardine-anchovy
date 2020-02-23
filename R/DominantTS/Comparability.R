# Comparability between sardine and anchovy 

# Data
load(here::here("R/DataCleanup/allsardineanchovy_3.RData"))# dataframe: alldat (raw data including NAs for missing years)

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

# RAM legacy
megatable <- vector()
for (r in 1:nregions){
  for (v in 1: nvariables){
    xx <- try(corr.fig(region_or_subregion=regions[r],
                       variable=variables[v],
                       scale = "Region",
                       data_source="RAM")) #must rotate sources
    mt <-  xx$median.table # previously: max.table
    mt$variable=rep(variables[v],times=nrow(mt))
    mt$region=rep(regions[r],times=nrow(mt))  
    mt$domsard <- ifelse(length(xx$Dom_sard_LTmax)==0,NA,xx$Dom_sard_LTmax)
    mt$domanch <- ifelse(length(xx$Dom_anch_LTmax)==0,NA,xx$Dom_anch_LTmax)
    megatable <- bind_rows(megatable,mt)
  }
}
RAM.summary <- megatable[-1,]
RAM.summary$datasource <- "RAM"

# Barange et al. 2009
megatable <- vector()
for (r in 1:nregions){
  for (v in 1: nvariables){
    xx <- try(corr.fig(region_or_subregion=regions[r],
                       variable=variables[v],
                       scale = "Region",
                       data_source="Barange")) #must rotate sources
    mt <-  xx$median.table # previously: max.table
    mt$variable=rep(variables[v],times=nrow(mt))
    mt$region=rep(regions[r],times=nrow(mt))  
    mt$domsard <- ifelse(length(xx$Dom_sard_LTmax)==0,NA,xx$Dom_sard_LTmax)
    mt$domanch <- ifelse(length(xx$Dom_anch_LTmax)==0,NA,xx$Dom_anch_LTmax)
    megatable <- bind_rows(megatable,mt)
  }
}
barange.summary <- megatable[-1,]
barange.summary$datasource <- "Barange"

# FAO landings
megatable <- vector()
for (r in 1:nregions){
  for (v in 1: nvariables){
    xx <- try(corr.fig(region_or_subregion=regions[r],
                       variable=variables[v],
                       scale = "Region",
                       data_source="FAO")) 
    mt <-  xx$median.table # previously: max.table
    mt$variable=rep(variables[v],times=nrow(mt))
    mt$region=rep(regions[r],times=nrow(mt))  
    mt$domsard <- ifelse(length(xx$Dom_sard_LTmax)==0,NA,xx$Dom_sard_LTmax)
    mt$domanch <- ifelse(length(xx$Dom_anch_LTmax)==0,NA,xx$Dom_anch_LTmax)
    megatable <- bind_rows(megatable,mt)
  }
}

FAO.summary <-megatable[-1,]
FAO.summary$datasource <- "FAO"

Both <- rbind(RAM.summary,barange.summary,FAO.summary)

#save(Both,file = "Replacement_RAM_Barange_MEDIAN2.Rdata") # For the new function that defines "dominant" species by which one has the highest median biomass


# Make Figure 2: comparability --------------------------------------------

# Plot log ratios of median biomass, landings, recruitment
load(here::here("ProcData/Replacement_RAM_Barange_MEDIAN2.Rdata")) #df: Both 

desired_order <- c("Northern Benguela anchovy","Northern Benguela sardine","Southern Benguela anchovy","Southern Benguela sardine","Anchovy South Africa","Sardine South Africa","California anchovy","California sardine","N Anchovy E Pacific","Pacific sardine Pacific Coast","Humboldt anchovy - Central Peru","Humboldt sardine - N Central Peru","Humboldt anchovy - South Peru N Chile","Humboldt sardine - South Peru N Chile","Chilean common sardine","Japanese anchovy","Japanese sardine","Bay of Biscay anchovy","European sardine")

summary <- filter(Both,
                  stock != "Humboldt anchovy - South Peru N Chile" &
                    stock != "Chilean common sardine" )
summary$stock2 <- factor(summary$stock,desired_order)

# Add percent differences to 'both' figure:
tp <- summary %>% 
  select(sp,stock,median.var,variable,region,datasource) # previously: max.var

tp2 <- dcast(tp,region+variable+datasource ~ sp,value.var="median.var",fun.aggregate = max,na.rm=T) # this will give errors but it's ok-- it's just bc of the cases where sardine is present but not anchovy (also: this can be changed btwn max and median depending on how you define "dominant" stocks)
tp2$Anchovy[which(is.inf(tp2$Anchovy))] <- NA  
tp2$Sardine[which(is.inf(tp2$Sardine))] <- NA  
tp3 <- tp2 %>% 
  mutate(log.diff = log10(Anchovy/Sardine)) 
tp3$percent.diff <- (tp3$Anchovy-tp3$Sardine) / rowMeans(cbind(tp3$Anchovy,tp3$Sardine),na.rm=T) * 100
tp3$zeroes <- 0
tp3$variable <- recode(tp3$variable,
                       fishing.mortality="Fishing mortality",
                       ssb="Spawning stock biomass",
                       landings="Landings",
                       rec="Recruitment")
tp3$newvar = factor(tp3$variable, levels=c("Recruitment",
                                           "Spawning stock biomass",
                                           "Landings",
                                           "Fishing mortality"))

# Check % diffs between sardine and anchovy
subset(tp3,newvar=="Landings") %>% 
  mutate(which.dom = ifelse(Sardine > Anchovy, "Sardine","Anchovy"))
subset(tp3,newvar=="Spawning stock biomass") %>% 
  mutate(which.dom = ifelse(Sardine > Anchovy, "Sardine","Anchovy"))
subset(tp3,newvar=="Recruitment") %>% 
  mutate(which.dom = ifelse(Sardine > Anchovy, "Sardine","Anchovy"))

(st <- tp3 %>% 
    filter(newvar != "Fishing mortality") %>%
    filter(!is.infinite(log.diff)) %>% 
    droplevels() %>%
    group_by(region,newvar) %>% 
    summarize(ml = median(log.diff,na.rm=T)) %>% 
    droplevels() %>%
    mutate(newvar_f=factor(newvar,levels=c("Landings",
                                           "Spawning stock biomass",
                                           "Recruitment"))) %>%
    as.data.frame() )

fig2 <- tp3 %>% 
  filter(newvar != "Fishing mortality") %>%
  droplevels() %>%
  filter(!is.inf(log.diff)) %>%
  mutate(newvar_f=factor(newvar,levels=c("Landings","Spawning stock biomass","Recruitment"))) %>%
  ggplot(aes(x=log.diff,y=region)) + 
  geom_vline(xintercept=0,lty=2) + 
  geom_point(size = 3,colour='lightgrey') +
  theme_classic(base_size = 14) %+replace% theme(strip.background  = element_blank()) +
  facet_wrap(~newvar_f,ncol = 1) + 
  geom_point(size=3,data=st,aes(x=ml,y=region,colour=region)) +
  scale_colour_manual("Region",values=pal)+
  theme(legend.position = "none") +
  xlab(expression('Log'['10']* '(Anchovy / Sardine)')) + 
  ylab("")

pdf(here::here("R/Figures/Figure2.pdf"),width = 8, height = 5,useDingbats = F)
fig2
dev.off()
