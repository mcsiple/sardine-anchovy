# Test - quantiles for distributions of WMR under asynchronous/synchronous dynamics
# Simulate time series as before


require(MASS)
require(mvcwt)

true.cov <- c(-0.25) #,-0.75,-0.5,0
rho <- c(0.9) # 0.8, 0.7
sims = 1:1000
sim.df <- data.frame(expand.grid(sims,true.cov,rho))

set.seed(123)
for(sim in 1:length(sims)){
asyn <- generate.sa(true.covar = true.cov)

x <- 1:nrow(asyn)
y <- asyn
# wave.fun = "Morlet"
scale.exp = 0.5; nscales = get.nscales(x)
min.scale = get.min.scale(x); max.scale = get.max.scale(x)
scales = log2Bins(min.scale, max.scale, nscales)
w = mvcwt(x, y)
mr = wmr(w)
#image(mr, reset.par = FALSE,xlab="Year")
#contour(mr, bound = NA, add = TRUE)
to.trim <- round(scales) # this is the number of cells to trim from the matrix (top and bottom)
mmat <- mr$z[,,1]
# NEW: trim cone of influence from z matrix (these values don't count!)
for(c in 1:ncol(mmat)){
  mmat[1:to.trim[c],c] <- NA
  mmat[nrow(mmat):(nrow(mmat)-to.trim[c]),c] <- NA
}

ind <- which(mr$y < 5)
trim.z.list <- mmat[,ind]
ind2 <- which(mr$y > 5 & mr$y < 10)
trim.z2.list <- mmat[,ind2]
ind3 <- which(mr$y > 10)
trim.z3.list <- mmat[,ind3]

sim.df$median[sim] <- median(trim.z.list,na.rm=T)
sim.df$q75[sim] <- quantile(trim.z.list,probs=0.75,na.rm=T)
sim.df$median2[sim] <- median(trim.z2.list,na.rm=T)
sim.df$q75_2[sim] <- quantile(trim.z2.list,probs=0.75,na.rm=T)
sim.df$median3[sim] <- median(trim.z3.list,na.rm=T)
sim.df$q75_3[sim] <- quantile(trim.z3.list,probs=0.75,na.rm=T)
}


apply(sim.df,MARGIN=2,FUN = median)

par(mfrow=c(3,1))
set.seed(123)
for(i in 1:3){
  asyn <- generate.sa(true.covar = -0.25)
  plot(1:nrow(asyn),asyn[,1],type='l')
  lines(1:nrow(asyn),asyn[,2],col="darkgrey")
  # scales = log2Bins(min.scale, max.scale, nscales)
  # w = mvcwt(x, y)
  # mr = wmr(w)
  # image(mr, reset.par = FALSE,xlab="Year")
  # contour(mr, bound = NA, add = TRUE)
}



