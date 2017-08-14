require(MASS)
require(MASS)
#matrix of autocorrelation terms 
rho <- matrix(c(0.8,0,0,0.8),nrow=2,byrow=T) # used to be 0.8 and 0.8
# Create variance covariance matrix
Sigma <- matrix(0, nrow=2, ncol =2)
diag(Sigma)<- c(0.8,0.7) #c(0.8,0.7)

Sigma[1,2]<- -0.75 * Sigma[1,1]  # what does this step do??
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

plot(1:100,y[,1], type = "l", lwd = 2, col = "black")
lines(1:100, y[,2], type = "l", lwd = 2, col = "gray")

