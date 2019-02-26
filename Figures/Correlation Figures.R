
load("~/Dropbox/Chapter3-SardineAnchovy/Datasets/allsardineanchovy_3.RData") # dataframe: alldat
figwd <- ("~/Dropbox/Chapter3-SardineAnchovy/Submission/Figures")
library(dplyr)
library(reshape2)
library(ggplot2)
library(beyonce)
library(RCurl)
library(RColorBrewer)

# Test data ---------------------------------------------------------------
filter(alldat,datasource=="Barange" & stock == "California anchovy")
unique(alldat$stock)

# Set region palette ------------------------------------------------------
pal <- beyonce_palette(18)[2:6]

# FIGURE 1: raw time series data ------------------------------------------
sacols <- c("#3288bd","#d53e4f") #blue = sardine, red = anchovy

mf <- melt(alldat,id.vars = c("datasource", "scientificname", "stock", "year", "sp", "region", "subregion"))
mf$variable <- recode(mf$variable,ssb="Biomass",rec="Recruitment",landings="Landings")

fig1 <- mf %>%
  group_by(datasource,scientificname,stock,sp,region,subregion,variable) %>% mutate(std.value=value/mean(value,na.rm=T)) %>% #as.data.frame() value/mean(value,na.rm=T)
  filter(variable %in% c("Biomass","Recruitment","Landings"))  %>%
    ggplot(aes(x=year,y=std.value,colour=sp,lty=datasource,group=stock)) +
      geom_line(lwd=0.6) +
      scale_color_manual(values=sacols) +
  facet_grid(region~variable,scales="free_y") +
  ylab("Standardized value") +
  theme_classic(base_size=14) %+replace% theme(strip.background  = element_blank())
     
pdf(file.path(figwd,"Figure1.pdf"),width = 10, height = 8,useDingbats = F)
fig1 
dev.off()


# FIGURE 2: schematic of spectral analysis --------------------------------
# This figure is composed in Adobe Illustrator

# Make a function to look at correlations and time series --------
# This function takes variable, region, and returns a time series of the two dominant spps, correlations, and ACFs of each ts

corr.fig <- function(data = alldat, region_or_subregion = "Benguela", scale = "Region", data_source = "Barange", variable = "landings",MARSS.cov = FALSE, plot = FALSE, ccf.calc = FALSE){
  # data = dataset including the time series you're interested in
  # Region = 1 of 5 LMEs (Benguela, California, NE Atlantic, Kuroshio-Oyashio, Humboldt)
  # variable = the variable (rec, biomass, or landings) 
  # Datasource = FAO, RAM, or Barange
  if(!region_or_subregion %in% subset(data,region==region_or_subregion | subregion==region_or_subregion)$region){stop("this dataset does not contain this region or subregion")}
  if (scale == "Region") {dataset <- filter(data, region == region_or_subregion & datasource==data_source)}
  if (scale == "Subregion") {dataset <- filter(data, subregion == region_or_subregion & datasource==data_source) }   
  if (!"Sardine" %in% dataset$sp | !"Anchovy" %in% dataset$sp){stop("missing one of the species")}
  
    #dataset <- subset(dataset,datasource == datasource)
    #print(dataset)
  # This is a little janky and old-fashioned because I wrote it a while ago. In the interest of time, I am not re-writing because it works.
  if(variable=="landings"){
    lt.maxes <- dataset %>% group_by(sp,stock) %>% summarize(max.var=round(max(landings,na.rm=T),2)) %>% as.data.frame()
    lt.medians <- dataset %>% group_by(sp,stock) %>% summarize(median.var=round(median(landings,na.rm=T),2)) %>% as.data.frame()
    } 
  if(variable=="ssb"){
    lt.maxes <- dataset %>% group_by(sp,stock) %>% summarize(max.var=round(max(ssb,na.rm=T),2)) %>% as.data.frame()
    lt.medians <- dataset %>% group_by(sp,stock) %>% summarize(median.var=round(median(ssb,na.rm=T),2)) %>% as.data.frame()
    }
  if(variable=="rec"){
    lt.maxes <- dataset %>% group_by(sp,stock) %>% summarize(max.var=round(max(rec,na.rm=T),2)) %>% as.data.frame()
    lt.medians <- dataset %>% group_by(sp,stock) %>% summarize(median.var=round(median(rec,na.rm=T),2)) %>% as.data.frame()
    }
  if(variable=="fishing.mortality"){
    lt.maxes <- dataset %>% group_by(sp,stock) %>% summarize(max.var=round(max(fishing.mortality,na.rm=T),2)) %>% as.data.frame()
    lt.medians <- dataset %>% group_by(sp,stock) %>% summarize(median.var=round(median(fishing.mortality,na.rm=T),2)) %>% as.data.frame()
    }

  print(lt.medians)
  #anchovy stats
  #if(!"Anchovy" %in% lt.maxes.sp){print("No anchovy time series")}
  a.only <- subset(lt.maxes, lt.maxes$sp=="Anchovy")
  a.max.ind <- which.max(a.only$max.var)
  lt.max.a <- a.only[a.max.ind,ncol(a.only)]
  lt.max.sp <- a.only$stock[a.max.ind]
  # Which anchovy species had biggest long term value for this time series (i.e., the "dominant anchovy species")
  
  #sardine stats
  # if(!"Sardine" %in% lt.maxes.sp){print("No sardine time series")
  #                                 lt.max.sp.sar = dom.s.ts =  NA}else{
  s.only <- subset(lt.maxes, lt.maxes$sp=="Sardine")
  s.max.ind <- which.max(s.only$max.var)
  lt.max.s <- s.only[s.max.ind,ncol(s.only)]
  lt.max.sp.sar <- s.only$stock[s.max.ind]     # "Dominant sardine species"
  
  dom.s.ts <- dataset[which(dataset$stock==lt.max.sp.sar),]
  dom.a.ts <- dataset[which(dataset$stock==lt.max.sp),]
  #   print(head(dom.a.ts))
  #   print(head(dom.s.ts))

  if(variable=="landings"){sar = dom.s.ts$landings
                           anch = dom.a.ts$landings} 
  if(variable=="ssb"){sar = dom.s.ts$ssb
                      anch = dom.a.ts$ssb} 
  if(variable=="rec"){sar = dom.s.ts$rec
                      anch = dom.a.ts$rec}
  if(variable=="fishing.mortality"){sar = dom.s.ts$fishing.mortality
                                    anch = dom.a.ts$fishing.mortality}
  #Plot the two dominant stocks/ts
  if (plot == TRUE){
  plot(dom.s.ts$year,sar,
       type="l",
       xlab="Year",ylab=variable,main=paste(c(region_or_subregion,variable)),
       xlim=c(min(dataset$year),max(dataset$year)),
       ylim=c(0,max(lt.maxes[,ncol(lt.maxes)])))
  #lines(dom.a.ts$year,anch,col='red')
  }
  
            # One correlation method: use MARSS to find covariance --------------------------------------------
            if(MARSS.cov == TRUE){
              source("Fill_NAs_SA.R")
              sard.mars <- FillNAs.ts(cbind(dom.s.ts$year,sar),
                                      startyear=min(c(dom.s.ts$year,dom.a.ts$year)),
                                      endyear=max(c(dom.s.ts$year,dom.a.ts$year)))
              anch.mars <- FillNAs.ts(cbind(dom.a.ts$year,anch),
                                      startyear=min(c(dom.s.ts$year,dom.a.ts$year)),
                                      endyear=max(c(dom.s.ts$year,dom.a.ts$year)))
              
              # sar and anch are rows, columns are years
              MAR.obj <- log(rbind(sard.mars[2,],anch.mars[2,]))    #Landings are log transformed
              colnames(MAR.obj) <- sard.mars[1,]
            
              model.sa=list()
              model.sa$Q="unconstrained"
              
              kem.sa = MARSS(MAR.obj, model=model.sa, control=list(maxit=1000)) 
              correlation = kem.sa$par$Q[2]/(sqrt(kem.sa$par$Q[3]) * sqrt(kem.sa$par$Q[1]))
              kem.sa.CIs = MARSSparamCIs(kem.sa,method="parametric",nboot=200)
              print(kem.sa.CIs)
              }
            
            else correlation = "no MARSS correlation calculated"
  
  if(ccf.calc==TRUE){
    ccf(x=sar, y=anch, na.action = na.pass,type = "correlation")
  }
  
  return(list(correlation = correlation,Dom_anch_LTmax=as.character(lt.max.sp),
             Dom_sard_LTmax=as.character(lt.max.sp.sar),
             dom.anch = dom.a.ts,dom.sard = dom.s.ts,
             max.table = lt.maxes,
             median.table = lt.medians))
  
          # print(kem.sa.CIs$par.upCI)
          # print(kem.sa.CIs$par.lowCI)
          #legend("topright",legend=c("sardine","anchovy"),lty=c(1,1),col=c('black','red'))
          #table(dataset$Land.Area,dataset$Scientific.name)
  
}  #End corr.fig function


# Test function
corr.fig(region_or_subregion="Benguela",variable="landings",scale = "Region",data_source="Barange")
corr.fig(region_or_subregion="California",variable="ssb",scale = "Region",data_source="RAM")


# Summarize peak variables for all the regions (to make comparison)

regions <- c("Benguela","California","Humboldt","Kuroshio-Oyashio","NE Atlantic")
variables <- c("ssb","rec","landings","fishing.mortality")
nregions <- length(regions)
nvariables <- length(variables)
megatable <- vector(length=6)

for (r in 1:nregions){
for (v in 1: nvariables){
  xx <- try(corr.fig(region_or_subregion=regions[r],variable=variables[v],scale = "Region",data_source="RAM"))
  mt <-  xx$max.table #xx$median.table   
  mt$variable=rep(variables[v],times=nrow(mt))
  mt$region=rep(regions[r],times=nrow(mt))  
  mt$domsard <- xx$Dom_sard_LTmax
  mt$domanch <- xx$Dom_anch_LTmax
  megatable <- rbind(megatable,mt)
}
}

RAM.summary <- megatable[-1,]
RAM.summary$datasource <- "RAM"
# RAM summary is missing some values; fill in with NA
RAM.summary <- RAM.summary[-which(is.infinite(RAM.summary$max.var)),]

barange.summary <- megatable[-1,]
barange.summary$datasource <- "Barange"

FAO.summary <-megatable[-1,]
FAO.summary$datasource <- "FAO"

Both <- rbind(RAM.summary,barange.summary,FAO.summary)
#Both <- Both[-which(is.infinite(Both$median.var)),]
Both <- Both[-which(is.infinite(Both$max.var)),]

#save(Both,file = "Replacement_RAM_Barange_MAX.Rdata")

load("~/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/Replacement_RAM_Barange_MAX.Rdata") #df: Both


# FIGURE 3: replacement; log-ratios of maximums and medians ---------------
# Plot peak biomass, landings, recruitment
#Change order of x axis tick labels so that sardines/anchovy in similar ecosystems are close together
desired_order <- c("Northern Benguela anchovy","Northern Benguela sardine","Southern Benguela anchovy","Southern Benguela sardine","Anchovy South Africa","Sardine South Africa","California anchovy","California sardine","N Anchovy E Pacific","Pacific sardine Pacific Coast","Humboldt anchovy - Central Peru","Humboldt sardine - N Central Peru","Humboldt anchovy - South Peru N Chile","Humboldt sardine - South Peru N Chile","Chilean common sardine","Japanese anchovy","Japanese sardine","Bay of Biscay anchovy","European sardine")

                summary <- filter(Both,stock != "Northern Benguela sardine" & 
                                            stock != "Humboldt sardine - South Peru N Chile" & 
                                            stock != "Humboldt anchovy - South Peru N Chile" & 
                                            stock != "Chilean common sardine" )
                                            #variable != "fishing.mortality" & 
                                            #variable != "rec")
                summary$stock2 <- factor(summary$stock,desired_order)
                # plot.barange <- filter(barange.summary,variable != "fishing.mortality")  #Take out fishing mortality
                # plot.barange$variable = recode(plot.barange$variable,ssb = "SSB",
                #                                          rec = 'Recruitment',
                #                                          landings = 'Landings')

summary$stock2 <- factor(summary$stock,desired_order)
# summary$variable <- recode(both.summary$variable,ssb = "SSB",
#                                                 rec = 'Recruitment',
#                                                 landings = 'Landings')

sa.col <- c("#ef8a62","#67a9cf")

# Add percent differences to 'both' figure:
tp <- summary %>% select(sp,stock,max.var,variable,region,datasource) #"max.var"
tp2 <- dcast(tp,region+variable+datasource ~ sp,value.var="max.var",fun.aggregate = sum,na.rm=T)
tp3 <- tp2 %>% mutate(log.diff = log(Anchovy/Sardine)) 
tp3$percent.diff <- (abs(tp3$Anchovy-tp3$Sardine) / rowMeans(cbind(tp3$Anchovy,tp3$Sardine),na.rm=T) )* 100
tp3$zeroes <- 0
tp3$variable <- recode(tp3$variable,
                       fishing.mortality="Fishing mortality",
                       ssb="Spawning stock biomass",
                       landings="Landings",
                       rec="Recruitment")
tp3$newvar = factor(tp3$variable, levels=c("Recruitment","Spawning stock biomass","Landings","Fishing mortality"))
#newpal <- c("darkgrey",beyonce_palette(78)[-1],"#F8A02E") #old palette

newblue <- c("#c6dbef","#4292c6","#08519c")
fig3 <- tp3 %>% 
        filter(variable != "Recruitment") %>%
        ggplot(aes(x=log.diff,y=region,colour=datasource)) + 
        geom_vline(xintercept=0,lty=2) + 
        geom_point(size = 3) +
        scale_colour_manual("Data source",values=rev(newblue)) +
        #scale_colour_brewer(palette = "Blues") +
        facet_wrap(~newvar,ncol = 1) + 
        theme_classic(base_size = 14) %+replace% theme(strip.background  = element_blank()) +
        xlim(-3,3) +
        xlab("Log(Anchovy / Sardine)") + 
        ylab("") 

(st <- tp3 %>% filter(!is.infinite(log.diff))%>% 
    group_by(region,newvar) %>% 
    summarise(ml = median(log.diff,na.rm=T)) %>% 
  filter(newvar != "Recruitment") %>%
  droplevels() )
  

fig3_option2 <- tp3 %>% 
      filter(variable != "Recruitment") %>%
      droplevels() %>%
      ggplot(aes(x=log.diff,y=region)) + 
      geom_vline(xintercept=0,lty=2) + 
      geom_point(size = 3,colour='lightgrey') +
      scale_colour_manual(values=rev(newblue)) +
      #scale_colour_brewer(palette = "Blues") +
      theme_classic(base_size = 14) %+replace% theme(strip.background  = element_blank()) +
      facet_wrap(~newvar,ncol = 1) + 
      geom_point(size=3,data=st,aes(x=ml,y=region,colour=region)) +
      scale_colour_manual("Region",values=pal)+
      xlim(-3,3) +
      theme(legend.position = "none") +
      xlab("Log(Anchovy / Sardine)") + 
      ylab("") 

pdf(file.path(figwd,"Figure3_option2.pdf"),width = 8, height = 5,useDingbats = F)
fig3_option2
dev.off()


# New replacement plot for defense slides
        #pdf("NewReplacementPlot_black2.pdf",width=10,height=3,useDingbats = FALSE)
        newpal <- c("lightyellow",beyonce_palette(78)[-1],"#F8A02E")
        ggplot(tp3,aes(x=log.diff,y=zeroes)) + #colour=region,shape = datasource
          geom_vline(xintercept=0,lty=2,colour="white") + 
          geom_point(size = 5,colour="white",alpha=0.7) + #
          scale_color_manual(values = newpal) +
          facet_wrap(~variable,ncol = 1) + 
          #theme_classic(base_size = 14) + 
          theme_black(base_size=12) +
          xlim(c(-3,3)) +
          xlab("Log(Anchovy / Sardine)") + 
          ylab("") +
          scale_y_discrete(breaks=c(-0.1,0,0.1),labels=c("","","")) +
          annotate("text",label="Sardine dominant",x=-1.75,y=.5,colour='white') + 
          annotate("text",label="Anchovy dominant",x=1.75,y=.5,colour='white')
        #dev.off()
              
        

# Back to bar plots -------------------------------------------------------
# Add percent differences to figure
sb <- plot.barange[,c("sp","stock","max.var","variable","region")]
sb2 <- dcast(sb,region+variable ~ sp,value.var="max.var")
sb2$Percent.diff <- abs(sb2$Anchovy - sb2$Sardine) / ((abs(sb2$Anchovy) + abs(sb2$Sardine ))/2) * 100
# Originally this was % change, so it could show directionality, but %diff is actually a different value
sb2$Absolute.diff <- sb2$Anchovy - sb2$Sardine
toplot <- merge(plot.barange,sb2)

setwd("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Figures")
pdf(file = "peaks_stocks_both.pdf",width = 6, height = 5,pointsize=16,useDingbats = FALSE)
ggplot(plot.both,aes(x=stock2,y=max.var,fill=sp,alpha=datasource)) + 
  geom_bar(stat="identity",position="dodge",colour="black") + 
  scale_fill_manual(values=sa.col) +
  scale_alpha_manual(values=c(0.9,0.3)) +
  facet_grid(variable~region,scales="free") +
  xlab("Stock") +
  theme_classic(base_size=10) + ylab("Long-term maximum") +
  theme(axis.text.x = element_text(angle = 90,hjust=1)) #+
  #geom_text(aes(x = 2, y = Inf,label = round(Percent.diff,1)),vjust = 1.1,size = 2.8)
dev.off()



# Same plot for RAM -------------------------------------------------------
plot.RAM <- RAM.summary
plot.RAM <- filter(plot.RAM,variable != "fishing.mortality")  #Take out fishing mortality
plot.RAM <- mutate(plot.RAM,variable = revalue(variable, c("ssb" = "SSB",
                                                           'rec' = 'Recruitment',
                                                           'landings' = 'Landings')))
sb3 <- dcast(plot.RAM,region+variable ~ sp,value.var="max.var")
sb3$Percent.diff <-( (sb3$Anchovy -sb3$Sardine) / sb3$Sardine )*100
  #abs(sb2$Anchovy - sb2$Sardine) / ((abs(sb2$Anchovy) + abs(sb2$Sardine ))/2) * 100
# Originally this was % change, so it could show directionality, but %diff is actually a different value
sb3$Absolute.diff <- sb3$Anchovy - sb3$Sardine
toplot <- merge(plot.RAM,sb3)

pdf("peak_stocks_RAM.pdf",width = 4,height=5,useDingbats = FALSE)
ggplot(toplot,aes(x=stock,y=max.var,fill=sp)) + 
  geom_bar(stat="identity",position="dodge",colour="black") + 
  scale_fill_manual("Species",values=sa.col) +
  facet_grid(variable~region,scales="free") +
  xlab("Stock") +
  theme_classic(base_size=10) + ylab("Long-term maximum") +
  theme(axis.text.x = element_text(angle = 90,hjust=1)) +
  geom_text(aes(x = 2, y = Inf,label = round(Percent.diff,1)),vjust = 1.1,size = 2.8)
dev.off()
# Do same plot, but scaled to max for each LME ----------------------------


#Find out max biomass, landings, etc for whole ecosystem
#Before I do maxes, I have to take out subregions w only one species
#That means Northern Benguela sardine, Chilean common sardine, and I'm only going to use the Central Peru stock from Humboldt Current, just for simplicity in the figure.
barange.summary2 <- filter(barange.summary,stock != "Northern Benguela sardine" & stock != "Humboldt sardine - South Peru N Chile" & stock != "Humboldt anchovy - South Peru N Chile" & stock != "Chilean common sardine" & variable != "fishing.mortality" & variable != "rec")
maxes <- ddply(barange.summary2, c("region","variable"),summarise, max.region =max(max.var))

# Scale to max for each ecosystem
yy <- merge(barange.summary2,maxes,by=c("region","variable"),all=TRUE)
yy$scaled.var <- yy$max.var/yy$max.region
barange.summary2 <- yy
# filter(barange.summary,stock %in% c("Southern Benguela sardine","Southern Benguela anchovy")

summary.biomass <- filter(barange.summary2,variable == "ssb")
summary.landings <- filter(barange.summary2,variable == "landings")
tiff(filename = "peaks_barange_landings_AFS145.tiff",width = 8.5, height = 3.5,units = "in", pointsize=16,res=300)

ggplot(summary.landings,aes(x=sp,y=scaled.var,fill=sp)) + 
  geom_bar(stat="identity",col='black') + 
  scale_fill_manual(values = c("#ef8a62","#67a9cf")) +
  facet_grid(variable~region,scales="free") +
  theme_classic() + ylab("Long-term peak (scaled to max)") +
  theme(axis.text.x = element_text(angle = 45,hjust=1)) 
  

dev.off()


# Sample time series for talks --------------------------------------------

# sample <- subset(alldat,datasource=="RAM" & region =="California" )
# setwd("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Figures")
# pdf("SampleCA_presentation.pdf",width=12,height = 7)
# ggplot(sample,aes(x=year,y=ssb/1000,colour=sp)) + 
#   geom_line(lwd=1.7) + 
#   scale_colour_manual(values = c("#ef8a62","#67a9cf")) +
#   theme_classic(base_size = 16) +
#   ylab("Biomass (x 1000 t)")
# dev.off()



# Distributions for null vs. observed data --------------------------------
obs <- get_obs(dat = RB,dsource = "Barange",reg = "California",var = "landings")
obstest <- get_wmr(obs$std_anchovy, obs$std_sardine)
nulltest1 = get_wmr(anchovy.ts=xx$Anchovy.surrogates[,1],sardine.ts=xx$Sardine.surrogates[,1])
nulltest2 = get_wmr(anchovy.ts=xx$Anchovy.surrogates[,2],sardine.ts=xx$Sardine.surrogates[,2])


get_large_null <- function(dat = RB,dsource = "Barange",reg = "California",var = "ssb",nsims){
  #generate as many surrogates as you need to get your full sims:
  yy <- get_surrogates(obs = get_obs(dat = dat,dsource = dsource,reg = reg,var = var),nsurrogates = nsims) 
  # Combine multiple null runs to get a good null dist:
  null.combined <- get_wmr(anchovy.ts=yy$Anchovy.surrogates[,1],sardine.ts=yy$Sardine.surrogates[,1])
  for (i in 2:nsims){
    newlist <- get_wmr(anchovy.ts=yy$Anchovy.surrogates[,i],sardine.ts=yy$Sardine.surrogates[,i])
    null.combined <- mapply(c, null.combined, newlist, SIMPLIFY=FALSE)
  }
  return(null.combined)
}

giantnull <- get_large_null(dat = RBF,dsource = "Barange",reg = "California",var = "landings",nsims=50)
str(giantnull)





distfig <- function(null, observed, return.dataframe = F){
  #' @description plots a histogram of observed WMRs vs. "null" WMR
  #' @param null is a list of three WMR matrices, each of which represents one time window (from a null surrogate)
  #' @param obs is a list of three WMR matrices, each of which represents one time window (observed time series)
  
  nl <- data.frame(ID = rep(names(null), sapply(null, length)),
                 wmrdens = unlist(null))
  obsl <- data.frame(ID = rep(names(observed), sapply(observed, length)),
                     wmrdens = unlist(observed))
  nl$ID_ord = factor(nl$ID, levels=c('less.than.5','five.ten','ten.plus'))
  obsl$ID_ord = factor(obsl$ID, levels=c('less.than.5','five.ten','ten.plus'))
  
  nl$ID_ord <- recode(nl$ID_ord, less.than.5 = "< 5 yr", five.ten = "5-10 yr",ten.plus="10+ yr" )
  obsl$ID_ord <- recode(obsl$ID_ord, less.than.5 = "< 5 yr", five.ten = "5-10 yr",ten.plus="10+ yr" )
  
  ggplot(nl,aes(x=wmrdens)) + 
    geom_density(fill="grey",colour="darkgrey",alpha=0.5) + 
    geom_density(data=obsl,fill="#41b6c4",colour="#41b6c4",alpha=0.5) +
    facet_wrap(~ID_ord,ncol = 1) +
    ylab("Density") +
    xlab("Wavelet modulus ratio (WMR)") +
    theme_classic(base_size=14)
  
  if(return.dataframe == TRUE){
    nl$nv <- "null"
    obsl$nv <- "obs"
    all <- rbind.fill(nl,obsl)
    return(all)
  }
}

test <- distfig(null = giantnull,observed = obstest,return.dataframe = T)

# Need to get data into the following format for analysis:
# Year Sardine.est Anchovy.est      region datasource variable
newmelt <- alldat %>%  
  select(datasource,stock,year,ssb,rec,landings,sp, region, subregion) %>%  
  melt(id.vars = c("datasource","stock","region", "subregion","year"))


# Get a giant dataframe of all the WMRs for each of the regions.
stocks <- alldat %>% 
  distinct(sp, stock, region, datasource) # get unique stock-region combos
n <- nrow(stocks) # there should be 15 combos of region and datasource

compares <- alldat %>% 
  distinct(region, datasource)


# LANDINGS ----------------------------------------------------------------
landings.wmrs <- vector()
for(s in 1:nrow(compares)){ 
  tryCatch ({
  obs <- get_obs(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "landings")
  obstest <- get_wmr(obs$std_anchovy, obs$std_sardine)
  giantnull <- get_large_null(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "landings",nsims=10)
  null=giantnull
  observed=obstest
  test <- distfig(null = giantnull,observed = obstest,return.dataframe = T)
  test$region <- compares$region[s]
  test$datasource <- compares$datasource[s]
  test$variable <- "landings"
  # NEED TO ADD THESE OR SIMILAR TO FIND OUT WHICH TWO STOCKS ARE BEING COMPARED
  #test$sardstock <- subset(stocks,region==compares$region[s] & datasource==compares$datasource[s] & sp =="Sardine")$stock
  #test$anchstock <- subset(stocks,region==compares$region[s] & datasource==compares$datasource[s] & sp =="Anchovy")$stock
  landings.wmrs <- rbind(landings.wmrs,test)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# BIOMASS -----------------------------------------------------------------
ssb.wmrs <- vector()
for(s in 1:nrow(compares)){ 
  tryCatch ({
    obs <- get_obs(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "ssb")
    obstest <- get_wmr(obs$std_anchovy, obs$std_sardine)
    giantnull <- get_large_null(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "ssb",nsims=10)
    null=giantnull
    observed=obstest
    test <- distfig(null = giantnull,observed = obstest,return.dataframe = T)
    test$region <- compares$region[s]
    test$datasource <- compares$datasource[s]
    test$variable <- "ssb"
    ssb.wmrs <- rbind(ssb.wmrs,test)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# RECRUITMENT -----------------------------------------------------------------
rec.wmrs <- vector()
for(s in 1:nrow(compares)){ 
  tryCatch ({
    obs <- get_obs(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "rec")
    obstest <- get_wmr(obs$std_anchovy, obs$std_sardine)
    giantnull <- get_large_null(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "rec",nsims=10)
    null=giantnull
    observed=obstest
    test <- distfig(null = giantnull,observed = obstest,return.dataframe = T)
    test$region <- compares$region[s]
    test$datasource <- compares$datasource[s]
    test$variable <- "rec"
    rec.wmrs <- rbind(rec.wmrs,test)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



# PLOTS -------------------------------------------------------------------
pal <- beyonce_palette(18)[2:6]
#tiff("WMRS_Rec.tiff",width = 8,height = 10,units = 'in',res = 150)
landings.wmrs %>% subset(datasource=="Barange" & nv=="obs") %>%
ggplot(aes(x=wmrdens,colour=region,fill=region)) + 
  geom_density(alpha=0.5,lwd=1.2,trim=F) + 
  scale_colour_manual(values=pal) +
  scale_fill_manual(values=pal) +
  facet_wrap(~ID_ord,ncol=1,scales = "free_y") +
  ylab("Density") +
  xlab("Wavelet modulus ratio (WMR)") +
  theme_classic(base_size=14) +
  ggtitle("landings")
#dev.off()



# Junk plots --------------------------------------------------------------
fs <- subset(alldat, !is.na(fishing.mortality)) %>% melt(id.vars=c("datasource","stock","year","sp","region","subregion","scientificname"))
colnum <- length(unique(fs$stock))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

fs %>% filter(variable=="fishing.mortality") %>% 
  ggplot(aes(x=year,y=value,colour=stock,lty=datasource)) + 
  geom_line(lwd=1.2) + 
  scale_colour_manual(values = getPalette(colnum)) + 
  theme_classic(base_size=14) +
  ylab("Fishing mortality")

