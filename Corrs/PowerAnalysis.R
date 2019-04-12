# Power Analysis!

basedir <- "C:/Users/siplem/Dropbox/Chapter3-SardineAnchovy"
#basedir <- "~/Dropbox/Chapter3-SardineAnchovy" # MAKE SURE BASEDIR MATCHES IN IN THE MAKE AUTOCORRELATED TS FILE
source(file.path(basedir,"Code_SA/sardine-anchovy/Corrs/make autocorrelated MVN time series.R"))
#"Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/Corrs/NullModel.R"
source(file.path(basedir,"Code_SA/sardine-anchovy/Corrs/NullModel.R"))

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

test.y <- generate.sa(true.covar = -.7,nyears=longyrs,diag.sigma=c(.8,.8))
plot(test.y[,1],type='l',ylim=c(-4,4),ylab="Biomass")
lines(test.y[,2],col='red')
ts_test <- list(std_anchovy=scale(test.y[,1]),
                std_sardine=scale(test.y[,2]))

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
true.d$timescale <- c("less5","five10","tenplus")

# Now sample from the longer time series for shorter lengths
obslengths <- seq(10,150,by=10)
yearindices <- 1:longyrs

for(l in 1:length(obslengths)){
  #start.ind <- sample(x = 1:(yearindices-obslengths[l])) # use this eventually
  start.ind <- 1
  obs <- ts_test[start.ind:obslengths[l]]
  giantnull.obs <- get_large_null2(timeseries=obs,nsims=50)
  wmr.obs <- test_wmr(obs = obs, null.combined = giantnull.obs)
  d.obs <- data.frame(diff=wmr.obs$diff,
                       CI.L=wmr.obs$CI.L,
                       CI.U=wmr.obs$CI.U,
                       timescale = wmr.obs$period,
                       obs.length = paste(obslengths[l]))
}
