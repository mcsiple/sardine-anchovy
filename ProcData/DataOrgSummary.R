
# Annoying code for summarizing data. It's very untidy but it will have to do for now.

#load("~/Dropbox/Chapter3-SardineAnchovy/Datasets/allsardineanchovy.RData")
setwd("~/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/Corrs")
source("ExtractMaxes.R")
data.all <- alldat
region.max.ssb <- matrix(NA,nrow=nrow(data.all),ncol=2)
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
         landings.st = landings / mean(landings,na.rm=T)) %>%
  as.data.frame()

dm.standardized <- dm.standardized[,-(5:8)]
md <- melt(dm.standardized, id.vars=c('datasource','scientificname','stock','sp','region','subregion','domanch','domsard','dom.a','dom.s','year'))
sa.col <- c("red","darkblue")
levels(md$variable)
levels(md$variable) = c("SSB","Recruitment","Landings")


figwd <- "/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Figures"
setwd(figwd)
pdf("Summary_Dominant_Stocks.pdf",width = 8,height = 7,useDingbats = FALSE)
ggplot(md, aes(x=year,y=value,colour=sp,linetype=datasource)) + 
  geom_line(lwd=0.6) + 
  scale_linetype("Data source") +
  facet_grid(region~variable,scales="free_y") +
  theme_classic() +
  xlab("Year") +
  ylab("Standardized value") +
  scale_colour_manual("Species",values = sa.col)
dev.off()


data.summary <- data.all[,-c(ncol(data.all)-1,ncol(data.all))] %>% 
  melt(id.vars=c("datasource","scientificname","stock","year","sp","region","subregion","domanch","domsard")) %>% 
  group_by(datasource,scientificname,stock,variable,sp,region,subregion,domanch,domsard) %>% 
  summarize(nyears.data = length(which(!is.na(value))),max.value = max(value,na.rm=T)) %>% 
  as.data.frame()
unique(data.summary[,c('datasource','region','domanch','domsard')])


#write.csv(data.summary,"Data-summary.csv")
#ts.summary <- data.summary %>% group_by(datasource,region,stock,nyears.data,domanch,domsard) %>% summarize()

#subset(alldat,datasource=="RAM" & scientificname == "Engraulis encrasicolus" & stock=="Anchovy ICES VIII" & sp=="Anchovy" & region=="NE Atlantic")
#write.csv(alldat,"Data-check.csv")

