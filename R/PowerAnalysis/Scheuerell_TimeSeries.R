# This code makes two autocorrelated time series that respond to asynchronous drivers. 
# Code is from Mark, 4/12/19

## function to return a discrete sine wave
sinw <- function(length, period, amp = 1) {
  amp * sin(2 * pi * seq(length) / period)
}

## redefine MVN
mvn <- MASS::mvrnorm

## length of ts
tt <- 100

## degree of autocorrelation for each ts
bb <- c(0.7, 0.7)

## desired correlation b/t 2 ts
rho <- 0.7

## SD's for both ts
sd <- c(0.9, 0.9)

## covariance matrix
## diagonal = variance
Sigma <- diag(sd^2)
## off-diag = covariance
Sigma[2,1] <- Sigma[1,2] <- rho * prod(sd)

## generate correlated errors
ee <- t(mvn(tt, rep(0, length(sd)), Sigma))

## sine wave as driver
sw <- sinw(tt, 40, 0.5)

## use the `ee` in a MARSS model
xx <- ee
BB <- diag(bb)
CC <- matrix(c(1, -1), ncol = 1)
cc <- matrix(sw, nrow = 1)
for(t in 2:tt) {
  xx[,t] <- BB %*% xx[,t-1] + CC %*% cc[,t] + ee[,t]
}

## how do they look?
plot.ts(t(xx))

## are they correlated?
cor(t(xx))