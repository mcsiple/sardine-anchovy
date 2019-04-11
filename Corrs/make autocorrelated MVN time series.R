#######################################################################################
# This code is modified from Tim w/ guidance from MDS. The code generates autocorrelated time series that covary with each other. You specify the covariance in the Sigma matrix, which is used to get multivariate process errors. It works nicely to create negatively covarying time series. After the initial test, there is code to figure out whether confidence intervals should be taken from resampling within the same simulated time series ("within"), or among different simulated time series ("between"). After that comparison is clear, we can pick one and use it to determine the amount of data necessary to detect a correlation, or the chance of seeing a spurious one. See "PowerAnalysis.R" for the submitted file.
#######################################################################################

require(MASS)
require(mvcwt)
require(reshape2)
require(ggplot2)

# basedir <- "~/Dropbox/Chapter3-SardineAnchovy"
# basedir <- "C:/Users/siplem/Dropbox/Chapter3-SardineAnchovy"
# figwd <- file.path(basedir,"Figures")

source("~/Dropbox/SFW/EBFMTesting/Phase 2/getSpectral.R")

    # # Matrix of autocorrelation terms
    # rho <- matrix(c(0.9,0,0,0.9),nrow=2,byrow=T) # you can get these values from estimating ACF[lag 1] from biomass estimates in RAM (or wherever)
    # # Create variance covariance matrix
    # Sigma <- matrix(0, nrow=2, ncol=2)
    # diag(Sigma) <- c(0.6,0.6) # These are the sigma values used for both anchovy and sardine in chapter 3
    # 
    # Sigma[1,2]<- -0.75 * Sigma[1,1]  # Add some negative correlation between the two time series
    # Sigma[2,1]<- -0.75 * Sigma[1,1]
    # 
    # # calculate variance of AR(1)
    # # beta is the kroneker product of the rho matrix 
    # beta <- rho  %x% rho
    # I <- matrix(0,nrow=nrow(beta),ncol=ncol(beta))
    # diag(I)<-1
    # Omega <- matrix(solve((I-beta)) %*% as.vector(Sigma) , 
    #                 nrow=2, 
    #                 ncol =2, 
    #                 byrow = F)
    # # Omega is the variance covariance of Y's
    # 
    # # Set up simulation
    # par(mfrow=c(1,1))
    # nyears <- 100
    # y <- matrix(NA, nrow = 100, ncol=2)
    # # get first year
    # mu1 = mu2 = 1 # means of each time series
    # y[1,] <- mvrnorm(1,mu=c(mu1,mu2),Sigma = Omega)
    # # TEE: Megsie seems like this should be: y<-mvrnorm(1, mu = c(mu1,mu2), Sigma = Omega)
    # 
    # eps <- mvrnorm(n = 100,mu = c(mu1,mu2), Sigma = Sigma)
    # 
    # 
    # # loop through years
    # for (i in 2:100) {
    #   eta.t1 <- y[i-1,]
    #   eta.t <- diag(rho) * eta.t1 + eps[i,]
    #   y[i,]<- eta.t
    # }
    # 
    # #test <- y
    # 
    # plot(1:100,y[,1], type = "l", lwd = 2, col = "red")
    # lines(1:100, y[,2], type = "l", lwd = 2, col = "pink")

# 
#Function to do the same
generate.sa <- function(diag.sigma = c(0.6,0.6), #these are the same result as acf()[1] gives, or very close. They do not match with ARIMA()$ar1
                        true.covar = 0,
                        nyears = 100, rho = matrix(c(0.8,0,0,0.8),nrow=2,byrow=T)){
  # Create variance covariance matrix
  Sigma <- matrix(0, nrow=2, ncol=2)
  diag(Sigma) <- diag.sigma             # These are the sigma values used for both anchovy and sardine in chapter 3
  Sigma[1,2]<- true.covar * Sigma[1,1]  # Add some negative correlation between the two time series
  Sigma[2,1]<- true.covar * Sigma[1,1]
  # beta is the kroneker product of the rho matrix 
  beta <- rho  %x% rho
  I <- matrix(0,nrow=nrow(beta),ncol=ncol(beta))
  diag(I)<-1
  Omega <- matrix(solve((I-beta)) %*% as.vector(Sigma), 
                  nrow=2, 
                  ncol =2, 
                  byrow = F)
  # Omega is the variance covariance of Y's
  y <- matrix(NA, nrow = nyears, ncol=2)
  # get first year
  mu1 = mu2 = 0
  y[1,] <- mvrnorm(1,mu=c(mu1,mu2),Sigma = Omega)
  eps <- mvrnorm(n = nyears,mu = c(0,0), Sigma = Sigma)
  # loop through years
  for (i in 2:nyears) {
    eta.t1 <- y[i-1,]
    eta.t <- diag(rho) * eta.t1 + eps[i,]
    y[i,]<- eta.t
  }
  return(y)
}

#generate.sa(true.covar = -.7,nyears=150)



# What are chances of spurious correlation? --------
# I am making a new version of this for publication, but the code below is still ok for old asynchrony detection method
# See "NullModel.R" for simulations with similar spectral characteristics (null expectation for amount of asynchrony observed)
# See how hard it is to detect asynchrony at different lengths of time series --------
# NOTE: THis one is a bust because the metric isn't sensitive enough. SO the new plan is to generate the null expectations for 0 correlation, then use K-S test to detect asynchrony in a series of simulations.
# result.df <- data.frame(true.covar = rep(c(-0.99,-0.75,-0.5,0,0.5,0.75,0.99),each=8),ts.length = rep(c(10,20,30,40,50,100,150,200)),d5 = NA,d510 = NA,d10=NA) 
# set.seed(123)
# nsims <- 500
# null1 <- null2 <- null3 <- list()
# prop.list1 <- prop.list2 <- prop.list3 <- list()
# 
# for(d in 1:nrow(result.df)){  #
#   # First: make null distribution for what you would expect at corr=0  --------
#   trim.z.list = trim.z2.list = trim.z3.list <- list()
#   tzl1 <- tzl2 <- tzl3 <- list()
#   prop1 = prop2 = prop3 = vector()
#   
#   for(s in 1:nsims){
#    
#     ts <- generate.sa(true.covar=0,
#                       nyears = result.df$ts.length[d]) #Generate null expectation based on length of ts
#     x <- 1:nrow(ts)
#     y <- ts
#     # wave.fun = "Morlet"
#      scale.exp = 0.5; nscales = get.nscales(x)
#      min.scale = get.min.scale(x); max.scale = get.max.scale(x)
#      scales = log2Bins(min.scale, max.scale, nscales)
#      w = mvcwt(x, y)
#      mr = wmr(w)
#      to.trim <- round(scales) # this is the number of cells to trim from the matrix (top and bottom)
#      mmat <- mr$z[,,1]
#      # OMG THIS TOOK ME SO LONG
#      for(c in 1:ncol(mmat)){
#        mmat[1:to.trim[c],c] <- NA
#        mmat[nrow(mmat):(nrow(mmat)-to.trim[c]),c] <- NA
#      }
#      
#     # NEW: trim cone of influence from z matrix (these values don't count!)
#     ind <- which(mr$y < 5)
#     trim.z.list <- mmat[,ind]
#     prop1[s] <- length(which(trim.z.list<0.7)) / length(which(!is.na(trim.z.list)))
#     tzl1 [[s]] <-  trim.z.list
#     
#     ind2 <- which(mr$y > 5 & mr$y < 10)
#     trim.z2.list <- mmat[,ind2]
#     prop2[s] <- length(which(trim.z2.list<0.7)) / length(which(!is.na(trim.z2.list)))
#     tzl2[[s]] <- trim.z2.list
#     
#     ind3 <- which(mr$y > 10)
#     if(length(ind3)==0){
#   		trim.z3 <- NA
#   		sim.vec10[s] <- NA
#   		result.df$d10[d]<- NA	
#   		prop3[s] <- NA
#   		}else{
#     		trim.z3.list <- mmat[,ind3]
#     		prop3[s] <- length(which(trim.z3.list<0.7)) / length(which(!is.na(trim.z3.list)))
#   		}
#     tzl3[[s]] <- trim.z3.list
#   }
#   
#   null1[[d]] <- plyr::ldply(tzl1, data.frame)
#   null2[[d]] <- plyr::ldply(tzl2,data.frame)
#   null3[[d]] <- plyr::ldply(tzl3,data.frame)
#   prop.list1[[d]] <- prop1 # each element of the list = 1 row of the results
#   prop.list2[[d]] <- prop2
#   prop.list3[[d]] <- prop3
#   print(d)
# }
# 
# 
# 
# sim.vec5 <- sim.vec510 <- sim.vec10 <- vector()
# 
# # Get upper limits based on that ts length
# for(d in 1:nrow(result.df)){
# limit1 <- ifelse(all(prop.list1[[d]]==0),NA,quantile(prop.list1[[d]],probs=0.975,na.rm=T))
# limit2 <- ifelse(all(prop.list2[[d]]==0),NA,quantile(prop.list2[[d]],probs=0.975,na.rm=T))
# limit3 <- ifelse(all(prop.list3[[d]]==0),NA,quantile(prop.list3[[d]],probs=0.975,na.rm=T))
# 
#   for(s in 1:nsims){ # simulate a bunch of runs with the true covariance, see if they have higher props under 0.5 than the null.
#   ts <- generate.sa(true.covar=result.df$true.covar[d],
#                     nyears = result.df$ts.length[d])
#   x <- 1:nrow(ts)
#   y <- ts
#   w = mvcwt(x, y)
#   mr = wmr(w)
#   
#   # Keep values inside cone of influence
#   to.trim <- round(scales) # this is the number of cells to trim from the matrix (top and bottom)
#   mmat <- mr$z[,,1]
#   # OMG THIS TOOK ME SO LONG
#   for(c in 1:ncol(mmat)){
#     mmat[1:to.trim[c],c] <- NA
#     mmat[nrow(mmat):(nrow(mmat)-to.trim[c]),c] <- NA
#   }
#   ind <- which(mr$y < 5)
#   trim.z <- mmat[,ind]
#   
#   ind2 <- which(mr$y > 5 & mr$y < 10)
#   trim.z2 <- mmat[,ind2]
#   
#   ind3 <- which(mr$y > 10)
#   
#   if(length(ind3)==0){
#   trim.z3 <- NA
#   sim.vec10[s] <- NA
#   result.df$d10[d]<- NA	
#   }else{
#   trim.z3 <- mmat[,ind3]
#   
#   sim.vec5[s] <- length(which(trim.z<0.7)) / length(which(!is.na(trim.z))) #prop. asynchronous at <5 yr timescale
#   sim.vec510[s] <- length(which(trim.z2<0.7)) / length(which(!is.na(trim.z2)))
#   sim.vec10[s] <- length(which(trim.z3<0.7)) / length(which(!is.na(trim.z3)))
#   }}
#   
#   result.df$d5[d] <- length(which(sim.vec5>limit1))/nsims
#   result.df$d510[d] <- length(which(sim.vec510>limit2))/nsims
#   result.df$d10[d] <- length(which(sim.vec10>limit3))/nsims
#   print(d)
#   }
#   
# #Running right now: rho = 0.8
# 
# head(result.df)
# save(result.df,file="PowerResults_likelyright_rho08.RData")
# 
# # Figure with detections --------------------------------------------------
# mdf <- melt(result.df,id.vars=c("true.covar","ts.length"))
# levels(mdf$variable) <- c("< 5 yr","5-10 yr","10+ yr")
# setwd(figwd)
# pdf("Power_08.pdf",width=8,height=3,useDingbats = FALSE)
# ggplot(mdf,aes(x=ts.length,y=value,colour=as.factor(true.covar),group=true.covar)) +
#   geom_line(lwd=1) + facet_grid(~variable) +
#   scale_color_brewer("True covariance",palette = 7,type = "div") + theme_classic(base_size=14) +
#   xlab("Time series length") + ylab("Probability of \n detecting asynchrony")
# dev.off()
# 
# #ggplot(result.df, aes(x=ts.length,y=d5,colour=true.covar))+ geom_point()
# 
# 
# 
# # Wavelet test ------------------------------------------------------------
# x <- 1:nrow(test)
# y <- test
# 
# plot(x,y[,1],type='l',lwd=1.5)
# lines(x,y[,2],col='darkgrey',lwd=1.5)
# 
# 
# # Remove cone of influence
# scale.exp = 0.5; nscales = get.nscales(x)
# min.scale = get.min.scale(x); max.scale = get.max.scale(x)
# scales = log2Bins(min.scale, max.scale, nscales)
# w = mvcwt(x, y)
# mr = wmr(w)
# image(mr, reset.par = FALSE,xlab="Year")
# contour(mr, bound = NA, add = TRUE)
# to.trim <- round(scales) # this is the number of cells to trim from the matrix (top and bottom)
# mmat <- mr$z[,,1]
# # OMG THIS TOOK ME SO LONG
# for(c in 1:ncol(mmat)){
#   mmat[1:to.trim[c],c] <- NA
#   mmat[nrow(mmat):(nrow(mmat)-to.trim[c]),c] <- NA
# }
# 
# 
# par(mfrow=c(3,1))
# var.color <- beyonce_palette(5)[5]
# 
# # Scale: <5 yr    
# ind <- which(mr$y < 5)
# trim.z <- mmat[,ind]
# hist(trim.z,xlim=c(0,1),col=var.color,border=var.color,
#      probability = T,main='',xaxt='n',xlab="",ylab="") 
# abline(v=0.7,lty=2)
# 
# #axis(1,at=c(0,0.5,1.0), labels=c(0,0.5,1.0))
# #text(tx,1,"<5 yr")
# 
# # Scale: 5-10 yr
# ind2 <- which(mr$y > 5 & mr$y < 10)
# trim.z2 <- mmat[,ind2]
# hist(trim.z2,xlim=c(0,1),col=var.color,border=var.color,
#      probability = T,main='',xaxt='n',xlab="")
# abline(v=0.7,lty=2)
# #axis(1,at=c(0,0.5,1.0), labels=c(0,0.5,1.0))
# #text(tx,1,"5-10 yr")
# 
# #Scale: 10+ yr
# ind3 <- which(mr$y > 10)
# trim.z3 <- mmat[,ind3]
# if(all(is.na(trim.z3))){plot.new()}else{
#   hist(trim.z3,xlim=c(0,1),col=var.color,border=var.color,
#        probability = T,main='',
#        xlab="Degree of synchrony",xaxt='n',ylab="")
#   #text(tx,1,"10+ yr")
#   axis(1,at=c(0,0.5,1.0), labels=c(0,0.5,1.0))}
# abline(v=0.7,lty=2)
# # See if KS test works for detecting differences between asynchronous time series and real ones.
# ks.test(x = trim.z3, # this is asynchrony at 10+ years, SIM'D from Tim code!
#         y = trim.z3b) # this is asynchrony from the real S/A time series
# 
# 
# # Figure out which method to use for detecting “power" --------------------
# # ANSWER: It doesn't matter whether you sample repeatedly within the same time series, or generate a whole new one.
# # This first one samples N-year chunks within one time series
# set.seed(123)
# within <- between <- matrix(NA,nrow = 20,ncol = 100) # ncol = sims, nrow = how many correlations
# true.covar <- seq(-1,1,length.out = 20)
# years.of.data <- 40     # How many years of data you have
# 
# for(c in 1:20){  #loop thru true covariances
#   Sigma <- matrix(0, nrow=2, ncol=2)
#   diag(Sigma) <- c(0.6,0.6) # These are the sigma values used for anchovy and sardine in chapter 3
#   
#   Sigma[1,2]<- true.covar[c] * Sigma[1,1]  # Add some negative correlation between the two time series
#   Sigma[2,1]<- true.covar[c] * Sigma[1,1]
#   
#   # calculate variance of AR(1)
#   # beta is the kroneker product of the rho matrix 
#   beta <- rho  %x% rho
#   I <- matrix(0,nrow=nrow(beta),ncol=ncol(beta))
#   diag(I)<-1
#   Omega <- matrix(solve((I-beta)) %*% as.vector(Sigma) , 
#                   nrow=2, 
#                   ncol =2, 
#                   byrow = F)
#   # Omega is the variance covariance of Y's
#   
#   # Simulation
#   nyears <- 100
#   y <- matrix(NA, nrow = 100, ncol=2)
#   # get first year
#   for (i in 1:2) y[1,i] <- rnorm(1,0,Omega)
#   eps <- mvrnorm(n = 100,mu = c(0,0), Sigma = Sigma)
#   
#   # loop through years
#   for (i in 2:100) {
#     eta.t1 <- y[i-1,]
#     eta.t <- diag(rho) * eta.t1 + eps[i,]
#     y[i,]<- eta.t
#   }
#   # Sample within the time series
#   # Loop through sims
#   for(sim in 1:100){
#     start.ind <- sample(x = 1:(100-years.of.data),size = 1)
#     ind.to.sample <- start.ind:(start.ind+years.of.data)
#     (test <- y[ind.to.sample,]) # For testing, print these a bunch
#     MAR.obj <- rbind(test[,1],test[,2])  
#     model.sa=list()
#     model.sa$Q="unconstrained" # estimate variance and covariance
#     kem.sa = MARSS(MAR.obj, model=model.sa, control=list(maxit=1000),silent = TRUE) 
#     covariance = kem.sa$par$Q[2]
#     correlation = covariance / (sqrt(kem.sa$par$Q[3]) * sqrt(kem.sa$par$Q[1]))
#     within[c,sim] <- covariance
#   }
# }
# 
# # It looks like within and between are similar, so we can use either method for simulating.
# # MARSS covariance --------------------------------------------------------
# # Test whether CIs for variance should be from within the same time series or between different ones.
# 
# # Among simulations
# for(c in 1:20){  #loop thru true covariances
#   Sigma <- matrix(0, nrow=2, ncol=2)
#   diag(Sigma) <- c(0.6,0.6) # These are the sigma values used for anchovy and sardine in chapter 3
#   
#   Sigma[1,2]<- true.covar[c] * Sigma[1,1]  # Add some negative correlation between the two time series
#   Sigma[2,1]<- true.covar[c] * Sigma[1,1]
#   
#   # calculate variance of AR(1)
#   # beta is the kroneker product of the rho matrix 
#   beta <- rho  %x% rho
#   I <- matrix(0,nrow=nrow(beta),ncol=ncol(beta))
#   diag(I) <- 1
#   Omega <- matrix(solve((I-beta)) %*% as.vector(Sigma) , 
#                   nrow=2, 
#                   ncol =2, 
#                   byrow = F)
#   # Loop through sims
#   for(sim in 1:100){
#     y <- matrix(NA, nrow = 100, ncol=2)
#     for (i in 1:2) y[1,i] <- rnorm(1,0,Omega)
#     eps <- mvrnorm(n = 100,mu = c(0,0), Sigma = Sigma)
#     # loop through years
#     for (i in 2:100) {
#       eta.t1 <- y[i-1,]
#       eta.t <- diag(rho) * eta.t1 + eps[i,]
#       y[i,]<- eta.t
#     }
#     test <- y[1:years.of.data,] 
#     MAR.obj <- rbind(test[,1],test[,2])  
#     model.sa=list()
#     model.sa$Q="unconstrained"
#     kem.sa = MARSS(MAR.obj, model=model.sa, control=list(maxit=1000),silent = TRUE) 
#     covariance = kem.sa$par$Q[2]
#     correlation = covariance / (sqrt(kem.sa$par$Q[3]) * sqrt(kem.sa$par$Q[1]))
#     between[c,sim] <- covariance
#   }
# }
# 
# # Figure  -----------------------------------------------------------------
# # Clean up this code later if worth it
# setwd("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData")
# #save(between,file = "Between_Var.RData")
# load("Within_Var.RData")
# meds <- apply(between,MARGIN = 1,FUN = median, na.rm=T)
# loCIs <- apply(between,MARGIN = 1, FUN = quantile, probs = c(0.05,0.95))[1,]
# hiCIs <- apply(between,MARGIN = 1, FUN = quantile, probs = c(0.05,0.95))[2,]
# 
# par(mar=c(4,4,3,2))
# plot(true.covar,meds, ylim=range(c(loCIs, hiCIs)),pch=19,type='b',xlab = "True covariance",ylab="Estimated covariance")
# abline(0,1,col="grey")
# arrows(true.covar, loCIs, true.covar, hiCIs, length=0.05, angle=90, code=3)
# 
# meds <- apply(within,MARGIN = 1,FUN = median, na.rm=T)
# loCIs2 <- apply(within,MARGIN = 1, FUN = quantile, probs = c(0.05,0.95), na.rm=T)[1,]
# hiCIs2 <- apply(within,MARGIN = 1, FUN = quantile, probs = c(0.05,0.95), na.rm=T)[2,]
# arrows(true.covar+0.02, loCIs2, true.covar+0.02, hiCIs2, length=0.05, angle=90, code=3)
# points(true.covar+0.02,meds, ylim=range(c(loCIs2, hiCIs2)),type='b',bg='white',pch=21)
# 
# legend("topleft",legend = c("Between simulations","Within simulations"),pch=c(19,1))
# 
# #kem.sa.CIs = MARSSparamCIs(kem.sa,method="parametric",nboot=200) # if you want confidence intervals for covariance estimates
# #print(kem.sa.CIs)
# 
# true.covar.vec <- c(-0.9,-0.5,0,0.5,0.9)
# tslength.vec <- 10:40
# results.mat <- matrix(NA,nrow=length(true.covar.vec),ncol=length(tslength.vec))
# 
# for (c in 1:length(true.covar.vec)){
#   for(t in 1:length(tslength.vec)){
#     test <- generate.sa(true.covar = true.covar.vec[c])
#     test <- test[1:tslength.vec[t],] 
#     # MAR.obj <- rbind(test[,1],test[,2])  
#     # model.sa=list()
#     # model.sa$Q="unconstrained"
#     # kem.sa = MARSS(MAR.obj, model=model.sa, control=list(maxit=1000),silent = TRUE) 
#     # covariance = kem.sa$par$Q[2]
#     # correlation = covariance / (sqrt(kem.sa$par$Q[3]) * sqrt(kem.sa$par$Q[1]))
#     # results.mat[c,t] <- correlation
#   }}
# 
# library(beyonce)
# print(beyonce_palette(11))
# pal <- beyonce_palette(11)
# pdf("TSLength_MARSS.pdf",width=11,height=10,useDingbats = FALSE)
# plot(tslength.vec,results.mat[1,],ylim=c(-1,1),type='l',col=pal[1],lwd=1.5)
# for(i in 2:nrow(results.mat)){
#   lines(tslength.vec,results.mat[i,],col=pal[i-1],lwd=1.5)
# }

