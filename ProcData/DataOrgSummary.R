#
# Annoying code for summarizing dominant species for each region and variable. 
# It's very untidy but it will have to do for now.
library(reshape2)
library(ggplot2)

source(here::here("R/DominantTS/ExtractMaxes.R"))
load(here::here("R/Data/allsardineanchovy_3.RData")) # dataframe alldat
data.all <- alldat

region.max.ssb <- matrix(NA,nrow=nrow(data.all),ncol=2)


# Test extract maxes function
extract.maxes(data = alldat,region_or_subregion = "Benguela",
              scale = "Region", # row 840 is the first FAO data point
              data_source = "Barange",
              variable = "landings" )

regions <- as.character(unique(data.all$region))
datasources <- as.character(unique(data.all$datasource))
dominant.summary <- expand.grid(region=regions,datasource=datasources)
dominant.summary$dominant.sardine <- NA
dominant.summary$dominant.anchovy <- NA
for(i in 10:nrow(dominant.summary)){
  get.max <- extract.maxes(data = data.all,
                           region_or_subregion = dominant.summary$region[i], 
                           scale="Region",
                           data_source = dominant.summary$datasource[i])  
  dominant.summary$dominant.anchovy[i]  = get.max$Dom_anch_LTmax
  dominant.summary$dominant.sardine[i]  = get.max$Dom_sard_LTmax
}
  
#  data.frame(region=NA,dominant.sardine=NA,dominant.anchovy=NA,datasource=NA)


for(i in 1:nrow(alldat)){
  region.max.ssb[i,1] <- as.character(extract.maxes(region_or_subregion = alldat$region[i], 
                                                    scale = "Region",
                                                    data_source = alldat$datasource[i],
                                                    variable = "ssb" )[1])
  region.max.ssb[i,2] <- as.character(extract.maxes(region_or_subregion = alldat$region[i], 
                                                    scale = "Region",
                                                    data_source = alldat$datasource[i],
                                                    variable = "ssb" )[2])
}

colnames(region.max.ssb) = c("domanch","domsard")
data.all <- cbind(alldat,region.max.ssb)
data.all$dom.a <- 0
data.all$dom.s <- 0
data.all$dom.a[which(data.all$stock == data.all$domanch)] = 1
data.all$dom.s[which(data.all$stock == data.all$domsard)] = 1

domstocks <- data.all %>% subset(dom.a==1 | dom.s==1)

dm.standardized <- domstocks %>% 
  group_by(datasource,scientificname,region,stock,domanch,domsard) %>% 
  mutate(ssb.st = ssb / mean(ssb,na.rm=T),
         rec.st = rec / mean(rec,na.rm=T),
         landings.st = landings / mean(landings,na.rm=T)) %>%as.data.frame()
        #select(datasource,scientificname,stock,sp,region,subregion,domanch,domsard,dom.a,dom.s,year, ssb,rec,landings)
  
md <- melt(dm.standardized, id.vars=c('datasource','scientificname','stock','sp','region','subregion','domanch','domsard','dom.a','dom.s','year'))
sa.col <- c("red","darkblue") # colors for plot

levels(md$variable)
levels(md$variable) = c("SSB","Rec","Landings","Fishing.mortality","Total.catch","SSB.st","Rec.st","Landings.st")


figwd <- "/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Figures"
setwd(figwd)
pdf("Summary_Dominant_Stocks.pdf",width = 8,height = 7,useDingbats = FALSE)
ggplot(md, aes(x=year,y=value,colour=datasource,linetype=sp)) + 
  geom_line(lwd=0.6) + 
  scale_linetype("Data source") +
  facet_grid(variable~region,scales="free_y") +
  theme_classic() +
  xlab("Year") +
  ylab("Standardized value")# +
  #scale_colour_manual("Species",values = sa.col)
dev.off()


data.summary <- data.all[,-c(ncol(data.all)-1,ncol(data.all))] %>% 
  melt(id.vars=c("datasource","scientificname","stock","year","sp","region","subregion","domanch","domsard")) %>% 
  group_by(datasource,scientificname,stock,variable,sp,region,subregion,domanch,domsard) %>% 
  summarize(nyears.data = length(which(!is.na(value))),max.value = max(value,na.rm=T)) %>% 
  as.data.frame()
unique(data.summary[,c('datasource','region','domanch','domsard')])


# Fix these later
#write.csv(data.summary,"Data-summary.csv")
#ts.summary <- data.summary %>% group_by(datasource,region,stock,nyears.data,domanch,domsard) %>% summarize()
#subset(alldat,datasource=="RAM" & scientificname == "Engraulis encrasicolus" & stock=="Anchovy ICES VIII" & sp=="Anchovy" & region=="NE Atlantic")
#write.csv(alldat,"Data-check.csv")


# TODO: Make sure dominant species function works with FAO data - should work as of 1/16/17 ---------------
load("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Datasets/FAO/FAO.Rdata")
ggplot(FAO,aes(x=year,y=landings,colour=Species,lty=sp)) + geom_line() + facet_wrap(~region)
