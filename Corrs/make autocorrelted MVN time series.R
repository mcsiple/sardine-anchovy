require(MASS)
require(MASS)
#matrix of autocorrelation terms 
rho <- matrix(c(0.8,0,0,0.8),nrow=2,byrow=T) # used to be 0.8 and 0.8
# Create variance covariance matrix
Sigma <- matrix(0, nrow=2, ncol=2)
diag(Sigma) <- c(0.6,0.6) # These are the sigma values used for anchovy and sardine in chapter 3

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



# MARSS covariance --------------------------------------------------------

# Test whether CIs for variance should be from within the same time series or between different ones.

within <- between <- vector(length = 100)
true.covar <- seq(-1,1,length.out = 100)

set.seed(123)
years.of.data <- 40 # How many years of data you have

for(i in 1:100){ 
  Sigma <- matrix(0, nrow=2, ncol=2)
  diag(Sigma) <- c(0.6,0.6) # These are the sigma values used for anchovy and sardine in chapter 3
  
  Sigma[1,2]<- true.covar[i] * Sigma[1,1]  # Add some negative correlation between the two time series
  Sigma[2,1]<- true.covar[i] * Sigma[1,1]
  
  # calculate variance of AR(1)
  # beta is the kroneker product of the rho matrix 
  beta <- rho  %x% rho
  I <- matrix(0,nrow=nrow(beta),ncol=ncol(beta))
  diag(I)<-1
  Omega <- matrix(solve((I-beta)) %*% as.vector(Sigma) , 
                  nrow=2, 
                  ncol =2, 
                  byrow = F)
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
  test <- y[1:years.of.data,] 
  MAR.obj <- rbind(test[,1],test[,2])    #Landings are log transformed
  model.sa=list()
  model.sa$Q="unconstrained"
  kem.sa = MARSS(MAR.obj, model=model.sa, control=list(maxit=1000),silent = TRUE) 
  covariance = kem.sa$par$Q[2]
  correlation = covariance / (sqrt(kem.sa$par$Q[3]) * sqrt(kem.sa$par$Q[1]))
  between[i] <- covariance # estimated covariance from one simulation
}

# Within
for (i in 1:100){
  Sigma <- matrix(0, nrow=2, ncol=2)
  diag(Sigma) <- c(0.6,0.6) # These are the sigma values used for anchovy and sardine in chapter 3
  
  Sigma[1,2]<- true.covar[i] * Sigma[1,1]  # Add some negative correlation between the two time series
  Sigma[2,1]<- true.covar[i] * Sigma[1,1]
  
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
}






      #kem.sa.CIs = MARSSparamCIs(kem.sa,method="parametric",nboot=200) # if you want to get confidence intervals for covariance estimates
      #print(kem.sa.CIs)





# Wavelet test ------------------------------------------------------------
# Currently this doesn't work for the simulated time series-- not sure why yet
x <- test[,1]
y <- test[,2]
w = mvcwt(x, y, min.scale = 1, max.scale = 20)
mr = wmr(w)

image(mr, reset.par = FALSE)
contour(mr, bound = NA, add = TRUE)

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

