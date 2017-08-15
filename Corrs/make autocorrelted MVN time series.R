#######################################################################################
# This code is modified from Tim w/ guidance from MDS. The code generates autocorrelated time series that covary with each other. You specify the covariance in the Sigma matrix, which is used to get multivariate process errors. It works nicely to create negatively covarying time series. After the initial test, there is code to figure out whether confidence intervals should be taken from resampling within the same simulated time series ("within"), or among different simulated time series ("between"). After that comparison is clear, we can pick one and use it to determine the amount of data necessary to detect a correlation, or the chance of seeing a spurious one.
#######################################################################################

require(MASS)
# Matrix of autocorrelation terms 
rho <- matrix(c(0.8,0,0,0.8),nrow=2,byrow=T) # used to be 0.8 and 0.8
# Create variance covariance matrix
Sigma <- matrix(0, nrow=2, ncol=2)
diag(Sigma) <- c(0.6,0.6) # These are the sigma values used for both anchovy and sardine in chapter 3

Sigma[1,2]<- -0.75 * Sigma[1,1]  # Add some negative correlation between the two time series
Sigma[2,1]<- -0.75 * Sigma[1,1]

# calculate variance of AR(1)
# beta is the kroneker product of the rho matrix 
beta <- rho  %x% rho
I <- matrix(0,nrow=nrow(beta),ncol=ncol(beta))
diag(I)<-1
Omega <- matrix(solve((I-beta)) %*% as.vector(Sigma) , 
                nrow=2, 
                ncol =2, 
                byrow = F)
# Omega is the variance covariance of Y's

# Set up simulation
par(mfrow=c(1,1))
nyears <- 100
y <- matrix(NA, nrow = 100, ncol=2)
# get first year
for (i in 1:2) y[1,i] <- rnorm(1,0,Omega)
eps <- mvrnorm(n = 100,mu = c(0,0), Sigma = Sigma)

# loop through years
for (i in 2:100) {
  eta.t1 <- y[i-1,]
  eta.t <- diag(rho) * eta.t1 + eps[i,]
  y[i,]<- eta.t
}

test <- y

plot(1:100,y[,1], type = "l", lwd = 2, col = "black")
lines(1:100, y[,2], type = "l", lwd = 2, col = "gray")

# 
#Function to do the same
generate.sa <- function(diag.sigma = c(0.6,0.6),
                        true.covar = 0,
                        nyears = 100){
  rho <- matrix(c(0.8,0,0,0.8),nrow=2,byrow=T) # used to be 0.8 and 0.8
  # Create variance covariance matrix
  Sigma <- matrix(0, nrow=2, ncol=2)
  diag(Sigma) <- diag.sigma # These are the sigma values used for both anchovy and sardine in chapter 3
  Sigma[1,2]<- true.covar * Sigma[1,1]  # Add some negative correlation between the two time series
  Sigma[2,1]<- true.covar * Sigma[1,1]
  # beta is the kroneker product of the rho matrix 
  beta <- rho  %x% rho
  I <- matrix(0,nrow=nrow(beta),ncol=ncol(beta))
  diag(I)<-1
  Omega <- matrix(solve((I-beta)) %*% as.vector(Sigma) , 
                  nrow=2, 
                  ncol =2, 
                  byrow = F)
  # Omega is the variance covariance of Y's
  y <- matrix(NA, nrow = 100, ncol=2)
  # get first year
  for (i in 1:2) y[1,i] <- rnorm(1,0,Omega)
  eps <- mvrnorm(n = 100,mu = c(0,0), Sigma = Sigma)
  # loop through years
  for (i in 2:100) {
    eta.t1 <- y[i-1,]
    eta.t <- diag(rho) * eta.t1 + eps[i,]
    y[i,]<- eta.t
  }
  return(y)
}


# What are chances of spurious correlation? --------
# Need wavelet way to detect negative correlation because covariance can be detected at really short time scales (<10 years). That makes sense.
par(mfrow=c(1,2))




# Wavelet test ------------------------------------------------------------
x <- 1:nrow(test)
y <- test
w = mvcwt(x, y, min.scale = 1, max.scale = 20)
mr = wmr(w)

image(mr, reset.par = FALSE)
contour(mr, bound = NA, add = TRUE)

par(mfrow=c(3,1))
# Scale: <5 yr    
ind <- which(mr$y < 5)
trim.z <- mr$z[nrow(mr$z)-ind,,1]
hist(trim.z,xlim=c(0,1),col=var.color,border=var.color,
     probability = T,main='',xaxt='n',xlab="",ylab="") 
text(tx,1,"<5 yr")

# Scale: 5-10 yr
ind2 <- which(mr$y > 5 & mr$y < 10)
trim.z2 <- mr$z[nrow(mr$z)-ind2,,1]
hist(trim.z2,xlim=c(0,1),col=var.color,border=var.color,
     probability = T,main='',xaxt='n',xlab="")
text(tx,1,"5-10 yr")

#Scale: 10+ yr
ind3 <- which(mr$y > 10)
trim.z3 <- mr$z[nrow(mr$z)-ind3,,1]
if(all(is.na(trim.z3))){plot.new()}else{
  hist(trim.z3,xlim=c(0,1),col=var.color,border=var.color,
       probability = T,main='',
       xlab="Degree of synchrony",xaxt='n',ylab="")
  text(tx,1,"10+ yr")
  axis(1,at=c(0,0.5,1.0), labels=c(0,0.5,1.0))}

# See if KS test works for detecting differences between asynchronous time series and real ones.
ks.test(x = trim.z3, # this is asynchrony at 10+ years, SIM'D from Tim code!
        y = trim.z3b) # this is asynchrony from the real S/A time series


# Figure out which method to use for detecting “power" --------------------
# This first one samples N-year chunks within one time series
set.seed(123)
within <- between <- matrix(NA,nrow = 20,ncol = 100) # ncol = sims, nrow = how many correlations
true.covar <- seq(-1,1,length.out = 20)
years.of.data <- 40     # How many years of data you have

for(c in 1:20){  #loop thru true covariances
  Sigma <- matrix(0, nrow=2, ncol=2)
  diag(Sigma) <- c(0.6,0.6) # These are the sigma values used for anchovy and sardine in chapter 3
  
  Sigma[1,2]<- true.covar[c] * Sigma[1,1]  # Add some negative correlation between the two time series
  Sigma[2,1]<- true.covar[c] * Sigma[1,1]
  
  # calculate variance of AR(1)
  # beta is the kroneker product of the rho matrix 
  beta <- rho  %x% rho
  I <- matrix(0,nrow=nrow(beta),ncol=ncol(beta))
  diag(I)<-1
  Omega <- matrix(solve((I-beta)) %*% as.vector(Sigma) , 
                  nrow=2, 
                  ncol =2, 
                  byrow = F)
  # Omega is the variance covariance of Y's
  
  # Simulation
  nyears <- 100
  y <- matrix(NA, nrow = 100, ncol=2)
  # get first year
  for (i in 1:2) y[1,i] <- rnorm(1,0,Omega)
  eps <- mvrnorm(n = 100,mu = c(0,0), Sigma = Sigma)
  
  # loop through years
  for (i in 2:100) {
    eta.t1 <- y[i-1,]
    eta.t <- diag(rho) * eta.t1 + eps[i,]
    y[i,]<- eta.t
  }
  # Sample within the time series
  # Loop through sims
  for(sim in 1:100){
    start.ind <- sample(x = 1:(100-years.of.data),size = 1)
    ind.to.sample <- start.ind:(start.ind+years.of.data)
    (test <- y[ind.to.sample,]) # For testing, print these a bunch
    MAR.obj <- rbind(test[,1],test[,2])  
    model.sa=list()
    model.sa$Q="unconstrained" # estimate variance and covariance
    kem.sa = MARSS(MAR.obj, model=model.sa, control=list(maxit=1000),silent = TRUE) 
    covariance = kem.sa$par$Q[2]
    correlation = covariance / (sqrt(kem.sa$par$Q[3]) * sqrt(kem.sa$par$Q[1]))
    within[c,sim] <- covariance
  }
}

# It looks like within and between are similar, so we can use either method for simulating.
# MARSS covariance --------------------------------------------------------
# Test whether CIs for variance should be from within the same time series or between different ones.

# Among simulations
for(c in 1:20){  #loop thru true covariances
  Sigma <- matrix(0, nrow=2, ncol=2)
  diag(Sigma) <- c(0.6,0.6) # These are the sigma values used for anchovy and sardine in chapter 3
  
  Sigma[1,2]<- true.covar[c] * Sigma[1,1]  # Add some negative correlation between the two time series
  Sigma[2,1]<- true.covar[c] * Sigma[1,1]
  
  # calculate variance of AR(1)
  # beta is the kroneker product of the rho matrix 
  beta <- rho  %x% rho
  I <- matrix(0,nrow=nrow(beta),ncol=ncol(beta))
  diag(I) <- 1
  Omega <- matrix(solve((I-beta)) %*% as.vector(Sigma) , 
                  nrow=2, 
                  ncol =2, 
                  byrow = F)
  # Loop through sims
  for(sim in 1:100){
    y <- matrix(NA, nrow = 100, ncol=2)
    for (i in 1:2) y[1,i] <- rnorm(1,0,Omega)
    eps <- mvrnorm(n = 100,mu = c(0,0), Sigma = Sigma)
    # loop through years
    for (i in 2:100) {
      eta.t1 <- y[i-1,]
      eta.t <- diag(rho) * eta.t1 + eps[i,]
      y[i,]<- eta.t
    }
    test <- y[1:years.of.data,] 
    MAR.obj <- rbind(test[,1],test[,2])  
    model.sa=list()
    model.sa$Q="unconstrained"
    kem.sa = MARSS(MAR.obj, model=model.sa, control=list(maxit=1000),silent = TRUE) 
    covariance = kem.sa$par$Q[2]
    correlation = covariance / (sqrt(kem.sa$par$Q[3]) * sqrt(kem.sa$par$Q[1]))
    between[c,sim] <- covariance
  }
}

# Figure  -----------------------------------------------------------------
# Clean up this code later if worth it
setwd("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData")
#save(between,file = "Between_Var.RData")
load("Within_Var.RData")
meds <- apply(between,MARGIN = 1,FUN = median, na.rm=T)
loCIs <- apply(between,MARGIN = 1, FUN = quantile, probs = c(0.05,0.95))[1,]
hiCIs <- apply(between,MARGIN = 1, FUN = quantile, probs = c(0.05,0.95))[2,]

par(mar=c(4,4,3,2))
plot(true.covar,meds, ylim=range(c(loCIs, hiCIs)),pch=19,type='b',xlab = "True covariance",ylab="Estimated covariance")
abline(0,1,col="grey")
arrows(true.covar, loCIs, true.covar, hiCIs, length=0.05, angle=90, code=3)

meds <- apply(within,MARGIN = 1,FUN = median, na.rm=T)
loCIs2 <- apply(within,MARGIN = 1, FUN = quantile, probs = c(0.05,0.95), na.rm=T)[1,]
hiCIs2 <- apply(within,MARGIN = 1, FUN = quantile, probs = c(0.05,0.95), na.rm=T)[2,]
arrows(true.covar+0.02, loCIs2, true.covar+0.02, hiCIs2, length=0.05, angle=90, code=3)
points(true.covar+0.02,meds, ylim=range(c(loCIs2, hiCIs2)),type='b',bg='white',pch=21)

legend("topleft",legend = c("Between simulations","Within simulations"),pch=c(19,1))

#kem.sa.CIs = MARSSparamCIs(kem.sa,method="parametric",nboot=200) # if you want confidence intervals for covariance estimates
#print(kem.sa.CIs)

true.covar.vec <- c(-0.9,-0.5,0,0.5,0.9)
tslength.vec <- 10:40
results.mat <- matrix(NA,nrow=length(true.covar.vec),ncol=length(tslength.vec))

for (c in 1:length(true.covar.vec)){
  for(t in 1:length(tslength.vec)){
    test <- generate.sa(true.covar = true.covar.vec[c])
    test <- test[1:tslength.vec[t],] 
    # MAR.obj <- rbind(test[,1],test[,2])  
    # model.sa=list()
    # model.sa$Q="unconstrained"
    # kem.sa = MARSS(MAR.obj, model=model.sa, control=list(maxit=1000),silent = TRUE) 
    # covariance = kem.sa$par$Q[2]
    # correlation = covariance / (sqrt(kem.sa$par$Q[3]) * sqrt(kem.sa$par$Q[1]))
    # results.mat[c,t] <- correlation
  }}

library(beyonce)
print(beyonce_palette(11))
pal <- beyonce_palette(11)
pdf("TSLength_MARSS.pdf",width=11,height=10,useDingbats = FALSE)
plot(tslength.vec,results.mat[1,],ylim=c(-1,1),type='l',col=pal[1],lwd=1.5)
for(i in 2:nrow(results.mat)){
  lines(tslength.vec,results.mat[i,],col=pal[i-1],lwd=1.5)
}

