# Power Analysis!

source(here::here("R/PowerAnalysis/AutocorrelatedTimeSeries.R"))
source(here::here("R/PowerAnalysis/Scheuerell_TimeSeries.R"))
source(here::here("R/PowerAnalysis/getSpectral.R"))
source(here::here("R/WMR/NullModel.R"))


# outputs of get_large_null look like this:
#  list(std_anchovy=std_anchovy,std_sardine=std_sardine)
#outputs of get_obs() look like this:
#  list(std_anchovy=std_anchovy,std_sardine=std_sardine)

#outputs of get_surrogates look like this:
    # nyears = length(std_anchovy) # both vectors should be the same length
    # a.sims <- s.sims <- matrix(NA, nrow = nyears,ncol = nsurrogates)
    # for(s in 1:nsurrogates){
    #   a.sims[,s] <- as.numeric(surrogate(std_anchovy,method = 'phase'))
    #   s.sims[,s] <- as.numeric(surrogate(std_sardine,method = 'phase'))
    # }
    # return(list(Anchovy.surrogates = a.sims,Sardine.surrogates = s.sims))

longyrs=150

#test.y <- generate.sa(true.covar = -.7,nyears=longyrs,diag.sigma=c(.8,.8))
test.y <- t(get.mds.ts(length = longyrs,
                       autocorrs = c(0.6,0.6),
                       rho = 0.4,
                       sds = c(.7,.7),
                       driver.period = 60,
                       driver.amp = 0.5,
                       CC = matrix(c(1, -1), ncol = 1)))

plot(test.y[,1],type='l',ylim=c(-4,4),ylab="Biomass")
lines(test.y[,2],col='red')
ts_test <- list(std_anchovy=as.numeric(scale(test.y[,1])),
                std_sardine=as.numeric(scale(test.y[,2])))

get_large_null2 <- function(timeseries=ts_test,nsims=50){ #just feed it 2 time series
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


giantnull2 <- get_large_null2(timeseries=ts_test,nsims=50)

# ok, this is the null dist for the true dist. so we calculate wmr for the long time series, and this will be the d[] value that the shorter time series should converge to. 

test <- test_wmr(obs = ts_test, null.combined = giantnull2)
true.d <- data.frame(diff=test$diff,CI.L=test$CI.L,CI.U=test$CI.U)
true.d$timescale <- c("less.than.5","five.ten","ten.plus")
true.d$pearson.corr <- cor(ts_test$std_anchovy,ts_test$std_sardine)
true.d

# Now sample from the longer time series for shorter lengths
obslengths <- seq(10,longyrs,by=10)
yearindices <- 1:longyrs
outs <- data.frame(diff=NA,CI.L=NA,CI.U=NA,timescale=NA,obs.length=NA)

for(l in 1:length(obslengths)){ #get errors for wilcox tests with short ts and long time periods, so start with longer ones.
  
  #start.ind <- sample(x = 1:(yearindices-obslengths[l])) # use this eventually
  start.ind <- 1
  obs <- list(std_anchovy = ts_test$std_anchovy[start.ind:obslengths[l]],
              std_sardine = ts_test$std_sardine[start.ind:obslengths[l]])
  pearson.cor <- cor(obs$std_anchovy,obs$std_sardine,method="pearson")
  giantnull.obs <- get_large_null2(timeseries=obs,nsims=50)
  wmr.obs <- test_wmr(obs = obs, null.combined = giantnull.obs)
  d.obs <- data.frame(diff=wmr.obs$diff,
                       CI.L=wmr.obs$CI.L,
                       CI.U=wmr.obs$CI.U,
                       timescale = wmr.obs$period,
                       obs.length = paste(obslengths[l]),
                       pearson.cor = pearson.cor)
  outs <- bind_rows(outs,d.obs)
  print(outs)
  print(paste("done with obslength",obslengths[l]))
  }

outs <- outs[-1,]
outs$obslength2 <-as.numeric(outs$obs.length)
#band <- true.d 
outs$truediff <- rep(true.d$diff,times=length(obslengths)-2)
outs$trueloCI <- rep(true.d$CI.L,times=length(obslengths)-2)
outs$truehiCI <- rep(true.d$CI.U,times=length(obslengths)-2)

outs$timescale = factor(outs$timescale, levels=c('less.than.5','five.ten','ten.plus'))

ggplot(outs,aes(x=obslength2,y=diff,colour=timescale)) +geom_line() +facet_wrap(~timescale) +
  geom_line(aes(y=truediff,colour=timescale) )+
  geom_ribbon(data=outs,aes(ymin=trueloCI,ymax=truehiCI,fill=timescale),colour=NA,alpha=0.5) +
  xlab("Length of observed time series") +
  ylab("Observed-null distribution \n for WMR (- means asynchrony)")
  
plot(outs$obslength2,outs$pearson.cor,xlab="Length of observations",ylab="Estimated Pearson correlation")
abline(h = true.d$pearson.corr[1],col='blue')
