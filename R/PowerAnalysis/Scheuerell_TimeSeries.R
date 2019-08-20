# This code makes two autocorrelated time series that respond to asynchronous drivers. 
# Code is from Mark, 4/12/19

## function to return a discrete sine wave
sinw <- function(length, period, amp = 1) {
  amp * sin(2 * pi * seq(length) / period)
}

## redefine MVN
mvn <- MASS::mvrnorm

get.mds.ts <- function(length=100,
                       autocorrs,rho,sds,
                       driver.period,driver.amp,
                       CC){
  #' @description makes two autocorrelated time series with correlation rho, their own autocorrelation and sd, and driven by a sine wave driver. Correlation of each ts to the driver is specified by CC.
  #' @param length length of desired time series (numeric vector)
  #' @param autocorrs degree of autocorrelation for each ts (numeric vector)
  #' @param rho correlation between the two time series (numeric)
  #' @param sds SDs for both time series (numeric vector)
  #' @param driver.period the period of the sine wave driver
  #' @param driver.amp the amplitude of the sine wave driver
  #' @param CC effect of the sim wave driver on the true abundance
  #' 
  tt <- length
  bb <- autocorrs
  sd <- sds
  
  ## Var-covar matrix
  Sigma <- diag(sd^2)
  ## off-diag = covariance
  Sigma[2,1] <- Sigma[1,2] <- rho * prod(sd)
  
  ## generate correlated errors
  ee <- t(mvn(tt, rep(0, length(sd)), Sigma))
  
  ## generate driver
  sw <- sinw(tt, driver.period, driver.amp)
  
  ## use the `ee` in a MARSS model
  xx <- ee
  BB <- diag(bb)
  
  cc <- matrix(sw, nrow = 1)
  for(t in 2:tt) {
    xx[,t] <- BB %*% xx[,t-1] + CC %*% cc[,t] + ee[,t]
  }
  return(xx)
}

#get.mds.ts(length = 100,autocorrs = c(0.6,0.6),rho = 0.4,sds = c(0.9,0.9),driver.period = 60,driver.amp = 0.5,CC=matrix(c(1, -1), ncol = 1))

# ## length of ts
# tt <- 100
# 
# ## degree of autocorrelation for each ts
# bb <- c(0.7, 0.7)
# 
# ## desired correlation b/t 2 ts
# rho <- 0.7
# 
# ## SD's for both ts
# sd <- c(0.9, 0.9)
# 
# ## covariance matrix
# ## diagonal = variance
# Sigma <- diag(sd^2)
# ## off-diag = covariance
# Sigma[2,1] <- Sigma[1,2] <- rho * prod(sd)
# 
# ## generate correlated errors
# ee <- t(mvn(tt, rep(0, length(sd)), Sigma))
# 
# ## sine wave as driver
# sw <- sinw(tt, 60, 0.5)
# 
# ## use the `ee` in a MARSS model
# xx <- ee
# BB <- diag(bb)
# CC <- matrix(c(1, -1), ncol = 1)
# cc <- matrix(sw, nrow = 1)
# for(t in 2:tt) {
#   xx[,t] <- BB %*% xx[,t-1] + CC %*% cc[,t] + ee[,t]
# }
# 
# ## how do they look?
# #plot.ts(t(xx))
# #lines(sw)
# 
# ## are they correlated?
# cor(t(xx))
# 
# plot(1:tt,xx[1,],type='l')
# lines(1:tt,xx[2,],col='red')
# lines(1:tt,sw,col='pink')
