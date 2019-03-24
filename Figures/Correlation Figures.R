basedir <- "C:/Users/siplem"


load(file.path(basedir,"Dropbox/Chapter3-SardineAnchovy/Datasets/allsardineanchovy_3.RData")) # dataframe: alldat
figwd <- file.path(basedir,"Dropbox/Chapter3-SardineAnchovy/Submission/Figures")
fnwd <- file.path(basedir,"/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/Corrs")
source(file.path(basedir,"Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/Fill_NAs_SA.R")) # This is for corr.fig
source(file.path(basedir,"Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/Corrs/NullModel.R"))

library(dplyr)
library(reshape2)
library(ggplot2)
library(beyonce) #if needed: devtools::install_github("dill/beyonce")
library(RCurl)
library(RColorBrewer)

# Test data ---------------------------------------------------------------
filter(alldat,datasource=="Barange" & stock == "California anchovy")
unique(alldat$stock)

# Set palettes  ------------------------------------------------------
# Region palette
pal <- beyonce_palette(18)[2:6]

# Sardine-anchovy palette
sacols <- c("#3288bd","#d53e4f") #blue = sardine, red = anchovy


# FIGURE 1: raw time series data ------------------------------------------

mf <- melt(alldat,id.vars = c("datasource", "scientificname", "stock", "year", "sp", "region", "subregion"))
mf$variable <- recode(mf$variable,ssb="Biomass",rec="Recruitment",landings="Landings")

fig1 <- mf %>%
  group_by(datasource,scientificname,stock,sp,region,subregion,variable) %>% 
  mutate(std.value=value/mean(value,na.rm=T)) %>% #as.data.frame() value/mean(value,na.rm=T)
  filter(variable %in% c("Biomass","Recruitment","Landings"))  %>%
    ggplot(aes(x=year,y=std.value,colour=sp,lty=datasource,group=stock)) +
      geom_line(lwd=0.6) +
      scale_color_manual(values=sacols) +
  facet_grid(variable~region,scales="free_y") +
  ylab("Standardized value") +
  theme_classic(base_size=14) %+replace% theme(strip.background  = element_blank())
     
pdf(file.path(figwd,"Fig1_option2.pdf"),width = 12, height = 4,useDingbats = F)
fig1 
dev.off()



# FIGURE 2: schematic of spectral analysis --------------------------------
# This figure is composed in Adobe Illustrator

source(file.path(fnwd,"Corr_Fig.R")) # This includes the corr.fig function
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
  xx <- try(corr.fig(region_or_subregion=regions[r],
                     variable=variables[v],
                     scale = "Region",
                     data_source="FAO"))
  mt <-  xx$median.table #SWITCH MEDIAN/MAX
  mt$variable=rep(variables[v],times=nrow(mt))
  mt$region=rep(regions[r],times=nrow(mt))  
  mt$domsard <- xx$Dom_sard_LTmax
  mt$domanch <- xx$Dom_anch_LTmax
  megatable <- rbind(megatable,mt)
}
}

RAM.summary <- megatable[-1,]
RAM.summary$datasource <- "RAM"
if(any(is.infinite(RAM.summary$median.var))){
  RAM.summary <- RAM.summary[-which(is.infinite(RAM.summary$median.var)),]} #SWITCH MEDIAN/MAX

barange.summary <- megatable[-1,]
barange.summary$datasource <- "Barange"

FAO.summary <-megatable[-1,]
FAO.summary$datasource <- "FAO"

Both <- rbind(RAM.summary,barange.summary,FAO.summary)

#Both <- Both[-which(is.infinite(Both$median.var)),] #switch median/max here for the one you want
#Both <- Both[-which(is.infinite(Both$max.var)),]

#save(Both,file = "Replacement_RAM_Barange_MAX.Rdata")
#save(Both,file = "Replacement_RAM_Barange_MEDIAN.Rdata")

load(file.path(basedir,"Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/Replacement_RAM_Barange_MAX.Rdata")) #df: Both


# FIGURE 3: replacement; log-ratios of maximums and medians ---------------
# Plot peak biomass, landings, recruitment
#Change order of x axis tick labels so that sardines/anchovy in similar ecosystems are close together
desired_order <- c("Northern Benguela anchovy","Northern Benguela sardine","Southern Benguela anchovy","Southern Benguela sardine","Anchovy South Africa","Sardine South Africa","California anchovy","California sardine","N Anchovy E Pacific","Pacific sardine Pacific Coast","Humboldt anchovy - Central Peru","Humboldt sardine - N Central Peru","Humboldt anchovy - South Peru N Chile","Humboldt sardine - South Peru N Chile","Chilean common sardine","Japanese anchovy","Japanese sardine","Bay of Biscay anchovy","European sardine")

                summary <- filter(Both,stock != "Northern Benguela sardine" & 
                                            stock != "Humboldt sardine - South Peru N Chile" & 
                                            stock != "Humboldt anchovy - South Peru N Chile" & 
                                            stock != "Chilean common sardine" )
                
                summary$stock2 <- factor(summary$stock,desired_order)
                
summary$stock2 <- factor(summary$stock,desired_order)

# Add percent differences to 'both' figure:
tp <- summary %>% select(sp,stock,max.var,variable,region,datasource) # CHANGE THIS BTWN MAX AND MEDIAN
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

newblue <- c("#c6dbef","#4292c6","#08519c")
fig3 <- tp3 %>% 
        filter(variable != "Recruitment") %>%
        ggplot(aes(x=log.diff,y=region,colour=datasource)) + 
        geom_vline(xintercept=0,lty=2) + 
        geom_point(size = 3) +
        scale_colour_manual("Data source",values=rev(newblue)) +
        facet_wrap(~newvar,ncol = 1) + 
        theme_classic(base_size = 14) %+replace% theme(strip.background  = element_blank()) +
        xlim(-3,3) +
        xlab("Log(Anchovy / Sardine)") + 
        ylab("") 

(st <- tp3 %>% 
    filter(!is.infinite(log.diff))%>% 
    group_by(region,newvar) %>% 
    summarise(ml = median(log.diff,na.rm=T)) %>% 
    #filter(newvar != "Recruitment") %>%
    droplevels() )
  

fig3_option2 <- tp3 %>% 
      droplevels() %>%
      ggplot(aes(x=log.diff,y=region)) + 
      geom_vline(xintercept=0,lty=2) + 
      geom_point(size = 3,colour='lightgrey') +
      #scale_colour_manual(values=pal) +
      theme_classic(base_size = 14) %+replace% theme(strip.background  = element_blank()) +
      facet_wrap(~newvar,ncol = 1) + 
      geom_point(size=3,data=st,aes(x=ml,y=region,colour=region)) +
      scale_colour_manual("Region",values=pal)+
      theme(legend.position = "none") +
      xlab("Log(Anchovy / Sardine)") + 
      ylab("")

#pdf(file.path(figwd,"Figure3_max_option2.pdf"),width = 8, height = 5,useDingbats = F)
fig3_option2
#dev.off()



# FIGURE 4: MARSS covariances ---------------------------------------------
# This is the plot that requires me to estimate a bunch of covariances between dominant sardine-anchovy pairs. 
# I think they are the same as before but need to find the code.


# Distributions for null vs. observed data --------------------------------
obs <- get_obs(dat = RBF,dsource = "Barange",reg = "California",var = "landings")
obstest <- get_wmr(obs$std_anchovy, obs$std_sardine)
nulltest1 = get_wmr(anchovy.ts=xx$Anchovy.surrogates[,1],sardine.ts=xx$Sardine.surrogates[,1])
nulltest2 = get_wmr(anchovy.ts=xx$Anchovy.surrogates[,2],sardine.ts=xx$Sardine.surrogates[,2])
# Functions used here are in the NullModel.R file
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
  giantnull <- get_large_null(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "landings",nsims=50)
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
    giantnull <- get_large_null(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "ssb",nsims=50)
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
    giantnull <- get_large_null(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "rec",nsims=50)
    null=giantnull
    observed=obstest
    test <- distfig(null = giantnull,observed = obstest,return.dataframe = T)
    test$region <- compares$region[s]
    test$datasource <- compares$datasource[s]
    test$variable <- "rec"
    rec.wmrs <- rbind(rec.wmrs,test)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



# FIGURE 5A -------------------------------------------------------------------
pdf(file.path(figwd,"Fig5b_landings_FAO.pdf"),width = 8,height = 10)
landings.wmrs %>% subset(datasource=="FAO" & nv=="obs") %>%
ggplot(aes(x=wmrdens,colour=region,fill=region)) + 
  geom_density(alpha=0.5,lwd=1.2,trim=F) + 
  scale_colour_manual(values=pal) +
  scale_fill_manual(values=pal)+
  facet_wrap(~ID_ord,ncol=1,scales = "free_y") +
  ylab("Density") +
  xlab("Wavelet modulus ratio (WMR)") +
  theme_classic(base_size=14) +
  ggtitle("Landings")
dev.off()


# FIGURE 5b ---------------------------------------------------------------
#pdf(file.path(figwd,"Fig5B_landings.pdf"),width = 8,height = 10)
landings.wmrs %>% subset(datasource=="Barange" & nv=="obs") %>%
  ggplot(aes(x=wmrdens,colour=region,fill=region)) + 
  geom_density(alpha=0.5,lwd=1.2,trim=F) + 
  scale_colour_manual(values=pal) +
  scale_fill_manual(values=pal) +
  facet_wrap(~ID_ord,ncol=1,scales = "free_y") +
  ylab("Density") +
  xlab("Wavelet modulus ratio (WMR)") +
  theme_classic(base_size=14) +
  ggtitle("Landings")
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




# SUPPLEMENTAL FIGURES ----------------------------------------------------
# FIGURE S1: bar plot ---------------------------------------------------------------
load(file.path(basedir,"Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/Replacement_RAM_Barange_MAX.Rdata")) #df: Both

Both$variable <- recode(Both$variable,
                        fishing.mortality="Fishing mortality",
                        ssb="Spawning \n stock biomass",
                        landings="Landings",
                        rec="Recruitment")

cleanup_stocks <- c("Engraulidae - Atlantic, Southeast",
                    "Round sardinella - Atlantic, Southeast",
                    "California anchovy - Pacific, Northeast",
                    "Engraulidae - Pacific, Northwest",
                    "Sardinella spp - Pacific, Northwest",
                    "Slender raindbow sardine - Pacific, Northwest",
                    "Stolephorus spp - Pacific, Northwest",
                    "Sardinella spp - Atlantic, Northeast")

# Need to remove time series that aren't the "dominant" stocks in that ecosystem
remove.x <- vector()
for(i in 1:nrow(Both)){
  if(Both$stock[i] %in% c(Both$domsard[i],Both$domanch[i])){remove.x[i] = T}else{remove.x[i] = F} 
}
Both <- Both[remove.x,]

figS2 <- 
  Both %>% 
  filter(!stock %in% cleanup_stocks) %>%
  filter(datasource=="FAO") %>%
  droplevels() %>%
  ggplot(aes(x=stock,y=max.var,fill=sp)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual("Species",values=sacols) +
  #geom_text(data=sb2,aes(x=stock,y=3,label=round(Percent.diff,1)),vjust=1.1,size=2.8) +
  facet_grid(variable~region,scales='free') +
  theme_classic(base_size=14)%+replace% theme(strip.background  = element_blank()) + 
  xlab("Dominant stock ") +
  ylab("Long-term maximum") +
  theme(axis.text.x = element_text(angle = 60,hjust=1)) +
  ggtitle("FAO landings") +
  theme(plot.margin = unit(c(0,0,0,2), "cm"))

pdf(file.path(figwd,"FAO_maxes.pdf"),width = 11, height=7,useDingbats = F)
figS2
dev.off()

# Add percent differences to figure ---------------------------------------
sb <- subset(Both, datasource=="Barange") %>% 
  select(sp,stock,max.var,variable,region)
sb2 <- dcast(sb,region+variable ~ sp,value.var="max.var")
sb2$Percent.diff <- abs(sb2$Anchovy - sb2$Sardine) / ((abs(sb2$Anchovy) + abs(sb2$Sardine ))/2) * 100
# Originally this was % change, so it could show directionality, but %diff is actually a different value
sb2$Absolute.diff <- sb2$Anchovy - sb2$Sardine



# FIGURE Sx: Plot maximum variables for all the stocks, including non-dominant ones  --------

# head(alldat)
# mf <- melt(alldat,id.vars = c("datasource", "scientificname", "stock", "year", "sp", "region", "subregion"))
# mf$variable <- recode(mf$variable,ssb="Biomass",rec="Recruitment",landings="Landings",totalcatch="Total catch")
# mf %>% group_by(datasource,scientificname,stock,sp,region,variable) %>% summarise(max.var=max(value,na.rm=T)) %>%
#   subset(datasource=="Barange") %>%
#   droplevels() %>%
#   ggplot(aes(x=stock,y=max.var,fill=sp)) + 
#   geom_bar(stat = "identity") + 
#   scale_fill_manual("Species",values=sacols) +
#   #geom_text(data=sb2,aes(x=stock,y=3,label=round(Percent.diff,1)),vjust=1.1,size=2.8) +
#   facet_grid(variable~region,scales='free') +
#   theme_classic(base_size=14)%+replace% theme(strip.background  = element_blank()) + 
#   xlab("Dominant stock ") +
#   ylab("Long-term maximum")  +
#   theme(axis.text.x = element_text(angle = 60,hjust=1))




# Same plot for RAM -------------------------------------------------------
plot.RAM <- RAM.summary
plot.RAM <- filter(plot.RAM,variable != "fishing.mortality")  #Take out fishing mortality
plot.RAM$variable  <- recode(plot.RAM$variable,
                             ssb = "SSB",
                             rec = 'Recruitment',
                             landings = 'Landings')
sb3 <- dcast(plot.RAM,region+variable ~ sp,value.var="max.var")
sb3$Percent.diff <-( (sb3$Anchovy -sb3$Sardine) / sb3$Sardine )*100
#abs(sb2$Anchovy - sb2$Sardine) / ((abs(sb2$Anchovy) + abs(sb2$Sardine ))/2) * 100
# Originally this was % change, so it could show directionality, but %diff is actually a different value
sb3$Absolute.diff <- sb3$Anchovy - sb3$Sardine
toplot <- merge(plot.RAM,sb3)

pdf("peak_stocks_RAM.pdf",width = 4,height=5,useDingbats = FALSE)
ggplot(toplot,aes(x=stock,y=max.var,fill=sp)) + 
  geom_bar(stat="identity",position="dodge",colour="black") + 
  scale_fill_manual("Species",values=sacols) +
  facet_grid(variable~region,scales="free") +
  xlab("Stock") +
  ylab("Long-term maximum") +
  theme_classic(base_size=10) + 
  theme(axis.text.x = element_text(angle = 90,hjust=1)) +
  geom_text(aes(x = 2, y = Inf,label = round(Percent.diff,1)),vjust = 1.1,size = 2.8)
dev.off()
# Do same plot, but scaled to max for each LME ----------------------------


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


# Test for difference in WMR distributions (table S1) ---------------------------------
# Perform tests with full null accross each combination of datasource, region, variable
compares.RBF <- RBF %>% distinct(region, datasource, variable) # get unique comparisons

# make table results from Wilcoxon-Mann-Whitney test comparing observed vs null WMRs
test.table <- vector()
for(s in 1:nrow(compares.RBF)){ 
  giantnull <- get_large_null(dat = RBF,dsource = compares.RBF$datasource[s],reg = compares.RBF$region[s],
                              var = compares.RBF$variable[s],nsims=50)
  test <- test_wmr(obs = get_obs(dat = RBF,dsource = compares.RBF$datasource[s],reg = compares.RBF$region[s],
                                 var = compares.RBF$variable[s]), null.combined = giantnull)
  test$region <- compares.RBF$region[s]
  test$datasource <- compares.RBF$datasource[s]
  test$variable <- compares.RBF$variable[s]
  test.table <- rbind(test.table,test)
}

test.table.out <- test.table %>% select(region, variable, datasource, period, U, N, p = p.value, d = diff, CI.L, CI.U, r = r_p)

test.table.out$period <- recode(test.table.out$period, less.than.5 = "< 5 yr", five.ten = "5-10 yr", ten.plus="10+ yr")
test.table.out$p <- round(test.table.out$p, 3)
test.table.out$p[test.table.out$p < 0.001] <- "< 0.001"
test.table.out$d <- round(test.table.out$d, 3)
test.table.out$CI.L <- round(test.table.out$CI.L, 3)
test.table.out$CI.U <- round(test.table.out$CI.U, 3)
test.table.out$r <- round(test.table.out$d, 2)
#test.table.out$r[test.table.out$r > 0 & test.table.out$r < 0.01] <- "< 0.01"
test.table.out <- arrange(test.table.out, region, datasource, variable)

write.csv(test.table.out, file = here::here("Figures/TableS1_test_WMR.csv"))
