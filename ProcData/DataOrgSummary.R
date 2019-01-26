
# Annoying code for summarizing data. It's very untidy but it will have to do for now.
rootdir <- ("C:/Users/siplem/Dropbox/")
setwd(file.path(rootdir,"Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/Corrs"))
source("ExtractMaxes.R")
library(reshape2)
library(ggplot2)

load(file.path(rootdir,"/Chapter3-SardineAnchovy/Datasets/allsardineanchovy.RData")) # all the data!
data.all <- alldat 

region.max.ssb <- matrix(NA,nrow=nrow(data.all),ncol=2)
# Test extract maxes function
extract.maxes(data = alldat,region_or_subregion = alldat$region[840],scale = "Region", # row 840 is the first FAO data point
              data_source = alldat$datasource[840],
              variable = "landings" )

for(i in 1:nrow(alldat)){
  region.max.ssb[i,1] <- as.character(extract.maxes(region_or_subregion = alldat$region[i], 
                                                    scale = "Region",
                                                    data_source = alldat$datasource[i],
                                                    variable = "ssb" )[1])
  region.max.ssb[i,2] <- as.character(extract.maxes(region_or_subregion = alldat$region[i], 
                                                    scale = "Region",
                                                    data_source = alldat$datasource[i],
                                                    variable = "ssb" )[2])
  if(alldat$datasource[i] =="FAO"){
    region.max.ssb[i,1] <- as.character(extract.maxes(region_or_subregion = alldat$region[i], 
                                                      scale = "Region",
                                                      data_source = alldat$datasource[i],
                                                      variable = "landings" )[1])
    region.max.ssb[i,2] <- as.character(extract.maxes(region_or_subregion = alldat$region[i], 
                                                      scale = "Region",
                                                      data_source = alldat$datasource[i],
                                                      variable = "landings" )[2])
  }
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
