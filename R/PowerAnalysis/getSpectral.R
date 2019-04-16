# Get spectral scaling parameters
# these are functions that came from Tim originally, they calculate spectral beta and sd for time series.



# functions to calculate tiime series properties --------------------------

# fast functon to calculate beta, does not use the multiple segment method
# spectrum() estimates spectral density
calc.beta.fast<-function(ts.data){
  tmp.spec<-spectrum(ts.data,plot=FALSE)
  freq<-tmp.spec$freq
  spec<-tmp.spec$spec
  tmp.data.frame<-data.frame(freq=log(freq),spec=log(spec))
  lm.out<-lm(spec~freq,data=tmp.data.frame)
  beta<--lm.out$coef[2]
  return(beta)
}

# Coefficient of variation (variability relative to mean)
cv <- function(x){
  val <- (sd(x,na.rm = T)/mean(x,na.rm=T))
  return(val)
}

get.autocorr <- function(x, ar=TRUE){
  #' @description uses arima to estimate the autocorrelation parameter
  #' @param x is a vector time series
  #' @param ar is whether or not you want the AR parameter estimate or sigma squared (kind of janky but it's easier this way with apply()!)
  #' The tricky thing is that I think arima assumes that the time series is stationary, which is def not true for some of these fish time series
  if(all(is.na(x))){return(list(ar1=NA,sig2=NA,acf1=NA))}
  
  if(ar){
  calc.ar <- arima(x = x,order = c(1,1,0)) # this is a first-order autoregressive model
  ar1 <- calc.ar$coef[[1]]
  sig2 <- calc.ar$sigma2 #This is sigma^2 not totally sure how to interpret
  ar.list <- list(ar1=ar1,sig2=sig2)
  return(ar.list)
  }
  
  # UPDATE: ACF might be a better measure of the extent of autocorrelation now that I have thought about it...
  else{
  calc.acf <- acf(x,plot = FALSE)
  acf1 <- calc.acf$acf[2] # lag 1
  acflist <- list(acf1=acf1)
  return(acflist)
  }

}


NA.to.mean <- function(x){
  y = x
  if(all(is.na(x))){ts.filled =rep("NA",times=length(x))}else{
    meanval <- mean(x,na.rm=T)
    y[which(is.na(x))] <- meanval
  }
  return(y)
}

