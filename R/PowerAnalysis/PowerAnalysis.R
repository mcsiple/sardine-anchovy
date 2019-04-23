# Power Analysis!

source(here::here("R/PowerAnalysis/AutocorrelatedTimeSeries.R"))
source(here::here("R/PowerAnalysis/Scheuerell_TimeSeries.R"))
source(here::here("R/PowerAnalysis/getSpectral.R"))
source(here::here("R/WMR/NullModel.R"))


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

par(mfcol=c(3,2))
hist(giantnull2$less.than.5,main="<5",xlab="null WMR")
hist(giantnull2$five.ten,main="5-10",xlab="null WMR")
hist(giantnull2$ten.plus,main="10+",xlab="null WMR")
obs <- get_wmr(anchovy.ts = ts_test$std_anchovy,sardine.ts = ts_test$std_sardine)
hist(obs$less.than.5,main="<5",xlab="WMR")
hist(obs$five.ten,main="5-10",xlab="WMR")
hist(obs$ten.plus,main="10+",xlab="WMR")

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
  #start.ind <- sample(x = 1:(longyrs-len),size = 1,replace = TRUE ) # use this eventually
  start.ind <- 1
  end.ind <- start.ind+(len-1)
  obs <- list(std_anchovy = ts_test$std_anchovy[start.ind:end.ind],
              std_sardine = ts_test$std_sardine[start.ind:end.ind])
  pearson.cor <- cor(obs$std_anchovy,obs$std_sardine,method="pearson")
  giantnull.obs <- get_large_null2(timeseries=obs,nsims=50)
  wmr.obs <- test_wmr(obs = obs, null.combined = giantnull.obs)
  d.obs <- data.frame(diff=wmr.obs$diff,
                       CI.L=wmr.obs$CI.L,
                       CI.U=wmr.obs$CI.U,
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
          scale_colour_manual("Species",values=sacols) +
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
