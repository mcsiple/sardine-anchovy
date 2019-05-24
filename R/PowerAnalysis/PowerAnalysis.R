# Power Analysis!

source(here::here("R/PowerAnalysis/AutocorrelatedTimeSeries.R"))
source(here::here("R/PowerAnalysis/Scheuerell_TimeSeries.R"))
source(here::here("R/PowerAnalysis/getSpectral.R"))
source(here::here("R/WMR/NullModel.R"))

library(ggridges) # for "joy" plot
library(forcats)


# same as get_large_null() but just accepts the ts as a list
get_large_null2 <- function(timeseries,nsims){ 
  #' @param timeseries a list of two time series
  #generate as many surrogates as you need to get your full sims:
  yy <- get_surrogates(obs = timeseries,nsurrogates = nsims)
  # Combine multiple null runs to get a good null dist:
  null.combined <- get_wmr(anchovy.ts=yy$Anchovy.surrogates[,1],sardine.ts=yy$Sardine.surrogates[,1])
  for (i in 2:nsims){ #this will take a long time
    newlist <- get_wmr(anchovy.ts=yy$Anchovy.surrogates[,i],sardine.ts=yy$Sardine.surrogates[,i])
    null.combined <- mapply(c, null.combined, newlist, SIMPLIFY=FALSE)
  }
  return(null.combined)
}



longyrs=150

test.y <- t(get.mds.ts(length = longyrs,
                       autocorrs = c(0.7,0.7),
                       rho = 0.1,
                       sds = c(0.8,0.8),
                       driver.period = 50,
                       driver.amp = 0.5,
                       CC = matrix(c(1, -1), ncol = 1)))


ts_test <- list(std_anchovy=as.numeric(scale(test.y[,1])),
                std_sardine=as.numeric(scale(test.y[,2])))

par(mfrow=c(1,1))
plot(1:longyrs,ts_test$std_anchovy, type='l',ylim=c(-4,4),ylab="Biomass")
lines(1:longyrs,ts_test$std_sardine,col='red')



giantnull2 <- get_large_null2(timeseries=ts_test,nsims=50)

 
# Wilcoxon test between distributions
test <- test_wmr(obs = ts_test, null.combined = giantnull2)
true.d <- data.frame(diff=test$diff,CI.L=test$CI.L,CI.U=test$CI.U)
true.d$timescale <- c("less.than.5","five.ten","ten.plus")
true.d$pearson.corr <- cor(ts_test$std_anchovy,ts_test$std_sardine)
true.d

# Now sample from the longer time series for shorter lengths
obslengths <- seq(10,longyrs,by=10)
yearindices <- 1:longyrs
outs <- data.frame(diff=NA,CI.L=NA,CI.U=NA,timescale=NA,obs.length=NA)

for(l in 3:length(obslengths)){ #get errors for wilcox tests with short ts and long time periods, so start with longer ones.
  len = obslengths[l]
  #start.ind <- sample(x = 1:(longyrs-len),size = 1,replace = TRUE ) # if you want to randomize
  start.ind <- 1
  end.ind <- start.ind+(len-1)
  obs <- list(std_anchovy = ts_test$std_anchovy[start.ind:end.ind],
              std_sardine = ts_test$std_sardine[start.ind:end.ind])
  pearson.cor <- cor(obs$std_anchovy,obs$std_sardine,method="pearson")
  giantnull.obs <- get_large_null2(timeseries=obs,nsims=50)
  wmr.obs <- test_wmr(obs = obs, null.combined = giantnull.obs)
  d.obs <- data.frame(diff = wmr.obs$diff,
                       CI.L = wmr.obs$CI.L,
                       CI.U = wmr.obs$CI.U,
                       timescale = wmr.obs$period,
                       obs.length = paste(len),
                       pearson.cor = pearson.cor)
  outs <- bind_rows(outs,d.obs)
  print(outs)
  print(paste("done with obslength",obslengths[l]))
  }

outs <- outs[-1,]
outs$obslength2 <-as.numeric(outs$obs.length)
outs$truediff <- rep(true.d$diff,times=length(obslengths)-2) #length(obslengths)-2
outs$trueloCI <- rep(true.d$CI.L,times=length(obslengths)-2) #length(obslengths)-2
outs$truehiCI <- rep(true.d$CI.U,times=length(obslengths)-2) #length(obslengths)-2

outs$timescale = factor(outs$timescale, levels=c('less.than.5','five.ten','ten.plus'))
levels(outs$timescale) <- c("<5 yr","5-10 yr",">10 yr")


twc <- "darkgrey" # twc stands for "trustworthy color" #"#3182bd"

powerplot <- ggplot(outs,aes(x=obslength2,y=diff)) +
              geom_line(aes(y=truediff),colour=twc)+
              geom_ribbon(data=outs,aes(ymin=trueloCI,ymax=truehiCI),
              alpha=0.5,colour=NA,fill=twc) +
              geom_line(data=outs,aes(x=obslength2,y=diff),lwd=1,colour=twc) +
              geom_ribbon(data=outs,aes(ymin=CI.L,ymax=CI.U),alpha=0.5) +
              facet_wrap(~timescale) +
              xlab("Number of years of observations") +
              ylab("Observed-null distribution \n for WMR (d)") +
              theme_classic(base_size=14) +
              theme(strip.background = element_blank())

# Plot the original time series
# Sardine-anchovy palette
sacols <- c("#3288bd","#d53e4f") #blue = sardine, red = anchovy
ldat <- data.frame(Year=1:150,Sardine=ts_test$std_sardine,Anchovy=ts_test$std_anchovy) %>%
  melt(id.vars="Year")
  
tsplot <- ggplot(ldat,aes(x=Year,y=value,colour=variable)) +
          geom_line(lwd=1) +
          scale_colour_manual("Species",values=rev(sacols)) +
          theme_classic(base_size=14) +
          ylab("\n  Standardized \n biomass ")


tiff(filename = here::here("R/Figures/PowerAnalysis.tiff"),
     width = 8.5,height=4.5,units = 'in',res = 250)

gridExtra::grid.arrange(tsplot, powerplot)

dev.off()  


par(mfrow=c(1,1))
plot(outs$obslength2,outs$pearson.cor,type='b',
     xlab="Length of observations",ylab="Estimated Pearson correlation")
abline(h = true.d$pearson.corr[1],col='blue') 


#save(true.d,file = here::here("TrueD_Asynchronous.RData"))


# Try to do the same but with multiple density curves instead of trying to summarize d --------
# Now sample from the longer time series for shorter lengths
obslengths2 <- seq(25,longyrs,by=25)
#yearindices <- 1:longyrs
bigouts <- data.frame()

for(l in 1:length(obslengths2)){ #get errors for wilcox tests with short ts and long time periods, so start with longer ones.
  len = obslengths2[l]
  start.ind <- 1
  end.ind <- start.ind+(len-1)
  obs <- list(std_anchovy = ts_test$std_anchovy[start.ind:end.ind],
              std_sardine = ts_test$std_sardine[start.ind:end.ind])
  pearson.cor <- cor(obs$std_anchovy,obs$std_sardine,method="pearson")
  #giantnull.obs <- get_large_null2(timeseries=obs,nsims=50)
  wmr.obs.raw <- get_wmr(anchovy.ts = ts_test$std_anchovy[start.ind:end.ind],
                         sardine.ts = ts_test$std_sardine[start.ind:end.ind])
  l5 <- data.frame(period="less.than.5",WMR=as.vector(wmr.obs.raw$less.than.5),N = len)
  l510 <- data.frame(period="five.ten",WMR=as.vector(wmr.obs.raw$five.ten),N = len)
  l10p <- data.frame(period="ten.plus",WMR=as.vector(wmr.obs.raw$ten.plus),N = len)
  
  outs <- bind_rows(l5,l510,l10p)
  bigouts <- bind_rows(bigouts,outs)
} # There will be errors about binding character and factor vectors but it's actually fine

str(bigouts)
bigouts$period <- as.factor(bigouts$period)
bigouts <- bigouts %>% mutate(period=fct_recode(period,
                                                `< 5 yr` = "less.than.5",
                                                `5-10 yr` = "five.ten",
                                                `10+ yr` = "ten.plus"))
bigouts$period <- forcats::fct_relevel(bigouts$period,"< 5 yr" )
  
medians <- bigouts %>% group_by(period,N) %>% summarize(median.wmr =median(WMR,na.rm=T)) %>% as.data.frame()
  
joyplot <- ggplot(bigouts,aes(y=as.factor(N))) +
  geom_density_ridges(aes(x=WMR),fill="navyblue",alpha=0.5,col="navyblue",panel_scaling = T) + 
  theme_classic(base_size=14) %+replace% theme(strip.background  = element_blank()) + 
  ylab("Time series length (number of years)") +
  coord_flip() +
  facet_wrap(~period)
joyplot

tiff(filename = here::here("R/Figures/PowerAnalysis_JoyPlot.tiff"),
     width = 8.5,height=4.5,units = 'in',res = 250)

gridExtra::grid.arrange(tsplot, joyplot)

dev.off()  



# Look for Type I errors --------------------------------------------------

test.typeI <- t(get.mds.ts(length = longyrs,
                       autocorrs = c(0.8,0.8),
                       rho = 0.1,
                       sds = c(0.7,0.7),
                       driver.period = 50,
                       driver.amp = 0.5,
                       CC = matrix(c(1, 1), ncol = 1))) # Notice these are the same sign (not opposite signs as above)


ts_test <- list(std_anchovy=as.numeric(scale(test.typeI[,1])),
                std_sardine=as.numeric(scale(test.typeI[,2])))
# Spearman correlations? Used by Izquierdo-Pena et al...
cor(ts_test$std_anchovy[1:50],ts_test$std_sardine[1:50],method="spearman")


corrs <- vector()
for(i in 2:length(ts_test$std_anchovy)){
  corrs[i] <- cor(ts_test$std_anchovy[1:i],ts_test$std_sardine[1:i],method="spearman")
}

plot(corrs,ylab="Spearman correlation",xlab="Time series length")


#  Try the same thing but w a bunch of simulations ------------------------
nsims <- 1000
corrs.mat <-matrix(NA,nrow = nsims,ncol=longyrs)
  
for(s in 1:nsims){
  test.typeI <- t(get.mds.ts(length = longyrs,
                                       autocorrs = c(0.8,0.8),
                                       rho = 0.1,
                                       sds = c(0.7,0.7),
                                       driver.period = 50,
                                       driver.amp = 0.5,
                                       CC = matrix(c(1, 1), ncol = 1))) 
ts_test <- list(std_anchovy=as.numeric(scale(test.typeI[,1])),
                std_sardine=as.numeric(scale(test.typeI[,2])))
corrs <- vector()
for(i in 2:length(ts_test$std_anchovy)){
  corrs[i] <- cor(ts_test$std_anchovy[1:i],ts_test$std_sardine[1:i],method="spearman")
}
corrs.mat[s,] <- corrs
}

plot(0:150,0:150,xlim=c(0,150),ylim=c(-1,1),type='n',ylab="Spearman correlation",xlab="Number of years with abundance for both species")
for(s in 1:nsims){
  lines(corrs.mat[s,],col='lightgrey')
}
mean.vec <- colMeans(corrs.mat)
median.vec <- apply(X = corrs.mat,FUN = median,MARGIN = 2)
lines(1:150,mean.vec)
#lines(1:150, median.vec,col="blue")
# Add median length of time series
(prob.typeI <- apply(corrs.mat,MARGIN = 2,FUN = function(x) {length(which(x<0))/length(x)} ) )

#length(which(corrs.mat[,2]<0))
library(beyonce)
pal <- beyonce_palette(18)[2:6]


spc <- data.frame(datasource=c("RAM","RAM",rep("Barange",times=5)),
                  region = c("NE Atlantic",
                             "California",
                             "Kuroshio-Oyashio",
                             "NE Atlantic",
                             "Humboldt",
                             "California",
                             "Benguela"),
                  years.of.data = c(21,27,11,18,20,11,24),
                  spearman.cor = c(-0.3,-0.31,-0.08,-0.34,-0.6,-0.65,-0.7),
                  color.to.plot = c(pal[5],pal[2],pal[4],pal[5],pal[3],pal[2],pal[1])
)
points(spc$years.of.data,spc$spearman.cor,pch=19,col=as.character(spc$color.to.plot),cex=2)

# What is the P(type I error) based on the length of ts we have? ----------
load(here::here("R/DataCleanup/allsardineanchovy_3.RData"))
actual.data <- alldat
# How many years of data do we have for both species in one ecosystem?
KO <- subset(actual.data,region=="Kuroshio-Oyashio" & datasource=="Barange")
Ben <- subset(actual.data,region=="Benguela" & datasource=="Barange")
CA <- subset(actual.data,region=="California" & datasource=="Barange")
Hum <- subset(actual.data,region=="Humboldt" & datasource=="Barange")
NEA <- subset(actual.data,region=="NE Atlantic" & datasource=="Barange")

unique(KO$stock)
KO%>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  filter(complete.cases(.)) %>%
  nrow() # This is the number of years of biomass data with both species
xx <-  KO %>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  filter(complete.cases(.))
cor(xx[,2],xx[,3],method="spearman")



unique(Ben$stock)
Ben %>% group_by(stock) %>% summarize(max(ssb,na.rm=T))
Ben%>% select(stock,year,ssb) %>%
  filter(stock %in% c("Northern Benguela sardine","Southern Benguela anchovy"))  %>%
  dcast(year~stock) %>%
  filter(complete.cases(.)) %>%
  nrow() # This is the number of years of biomass data with both species
xx <-  Ben %>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  filter(complete.cases(.))
cor(xx[,2],xx[,3],method="spearman")


unique(CA$stock)
CA %>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  filter(complete.cases(.)) %>%
  nrow() # This is the number of years of biomass data with both species
xx <-  CA %>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  filter(complete.cases(.))
cor(xx[,2],xx[,3],method="spearman")


unique(Hum$stock)
Hum %>% group_by(stock) %>% summarize(max(ssb,na.rm=T))
Hum%>% select(stock,year,ssb) %>%
  filter(stock %in% c("Humboldt anchovy - Central Peru","Humboldt sardine - South Peru N Chile"))  %>%
  dcast(year~stock) %>%
  filter(complete.cases(.)) %>%
  nrow() # This is the number of years of biomass data with both species
xx <-  Hum %>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  filter(complete.cases(.))
cor(xx[,2],xx[,3],method="spearman")



unique(NEA$stock)
NEA %>% group_by(stock) %>% summarize(max(ssb,na.rm=T))
NEA%>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  filter(complete.cases(.)) %>%
  nrow() # This is the number of years of biomass data with both species
xx <-  NEA %>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  filter(complete.cases(.))
cor(xx[,2],xx[,3],method="spearman")


# Number of years with data for both s and a
# Kuroshio-Oyashio: 11 years
# Benguela: 24
# California: 11 
# Humboldt: 20
# NE Atlantic: 18


# Check the same thing for RAM legacy database ----------------------------
# How many years of data do we have for both species in one ecosystem?
KO <- subset(actual.data,region=="Kuroshio-Oyashio" & datasource=="RAM")
Ben <- subset(actual.data,region=="Benguela" & datasource=="RAM")
CA <- subset(actual.data,region=="California" & datasource=="RAM")
Hum <- subset(actual.data,region=="Humboldt" & datasource=="RAM")
NEA <- subset(actual.data,region=="NE Atlantic" & datasource=="RAM")

unique(KO$stock)
KO %>% group_by(stock) %>% summarize(max(ssb,na.rm=T))
# No sardine/anchovy pair available from RAM

unique(Ben$stock)
Ben %>% group_by(stock) %>% summarize(max(ssb,na.rm=T))
# No sardine/anchovy pair available from RAM

unique(CA$stock)
CA %>% group_by(stock) %>% summarize(max(ssb,na.rm=T))
CA %>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  filter(complete.cases(.)) %>%
  nrow() # This is the number of years of biomass data with both species
xx <- CA%>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  filter(complete.cases(.))
cor(xx[,2],xx[,3],method="spearman")



unique(Hum$stock)
Hum %>% group_by(stock) %>% summarize(max(ssb,na.rm=T))
# No sardine/anchovy pair available from RAM

unique(NEA$stock)
NEA %>% group_by(stock) %>% summarize(max(ssb,na.rm=T))
NEA%>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  filter(complete.cases(.)) %>%
  nrow() # This is the number of years of biomass data with both species
xx <- NEA%>% select(stock,year,ssb) %>%
  dcast(year~stock) %>%
  
  filter(complete.cases(.))
cor(xx[,2],xx[,3],method="spearman")

# From RAM legacy database, there are only a couple ts lengths.
# CA: 27 years
# NEA: 21 years
# Among the other stocks, there aren't sardine-anchovy pairs, either one 


# Plot the probability of a Type 1 error vs. how long our time ser --------

tiff(file="R/Figures/TypeI.tiff",width = 5,height = 5,units="in",res = 200)
plot(2:50,prob.typeI[2:50],type='l',ylab="Probability of false detection",col="grey",lwd=2,xlab="Number of years with biomass data for both species")
points(spc$years.of.data,rep(0,times=7),col=as.character(spc$color.to.plot),cex=2,pch=19)
dev.off()
