# Test - quantiles for distributions of WMR under asynchronous/synchronous dynamics
# Simulate time series as before


require(MASS)
require(mvcwt)
figwd <- "/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Figures"
true.covar.vec <- c(-0.9,-0.75,-0.5,-0.25,0)
rho.vec <- c(0.9,0.8,0.7)
sim.all=vector(length=9)

for(r in 1:length(rho.vec)){
  rho <- rho.vec[r]
  rho.sim = matrix(c(rho,0,0,rho),nrow=2,byrow=T)
    for(t in 1:length(true.covar.vec)){
    true.cov <- true.covar.vec[t]
    sims = 1:1000
    sim.df <- data.frame(expand.grid(sims,true.cov,rho.sim[1,1]))
    set.seed(123)
    for(sim in 1:length(sims)){
        asyn <- generate.sa(true.covar = true.cov,rho = rho.sim)
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
    sim.all <- rbind(sim.all,sim.df)
}}

sim.all <- sim.all[-1,]
apply(sim.df,MARGIN=2,FUN = median)

par(mfrow=c(3,1))
set.seed(123)
for(i in 1:3){
  asyn <- generate.sa(true.covar = -0.5,rho=rho.sim)
  plot(1:nrow(asyn),asyn[,1],type='l')
  lines(1:nrow(asyn),asyn[,2],col="darkgrey")
  # scales = log2Bins(min.scale, max.scale, nscales)
  # w = mvcwt(x, y)
  # mr = wmr(w)
  # image(mr, reset.par = FALSE,xlab="Year")
  # contour(mr, bound = NA, add = TRUE)
}

#head(sim.df)


# Plot all the time series with different covariances ---------------------
setwd(figwd)
pdf("ExampleCovariance.pdf",width=9,height=4,useDingbats = FALSE)
par(mfrow=c(1,1))
set.seed(123)
asyn <- generate.sa(true.covar = -0.9,rho=rho.sim)
plot(1:nrow(asyn),asyn[,1],type='l',col="darkblue",xlim=c(0,120))
lines(1:nrow(asyn),asyn[,2],col="red")

true.covar.vec <- c(-0.9,-0.75,-0.5,-0.25)
alpha.vec <- c(0.75,0.5,0.25,0.1)
for(i in 1:4){ 
set.seed(123)
asyn <- generate.sa(true.covar = true.covar.vec[i],rho=rho.sim)
lines(1:nrow(asyn),asyn[,1],col=adjustcolor("darkblue",alpha.f=alpha.vec[i]))
lines(1:nrow(asyn),asyn[,2],col=adjustcolor("red",alpha.f=alpha.vec[i]))
}
legend("right",lty=c(1,1),col=c("darkblue","red"),legend = c("Sardine","Anchovy"))

dev.off()


# Plot a bar plot for the 75% quantiles of all rho and covar values -------
# Requires sim.all dataframe
head(sim.all)
colnames(sim.all)[1:3] <- c("sim","true.covar","rho")
sim.all2 <- sim.all[,-c(4,6,8)]
ms <- melt(sim.all2,id.vars=c("sim",'true.covar','rho'))

pdf("Boxplot_cutoff.pdf",width=7,height = 9,useDingbats = FALSE)
ggplot(ms, aes(x = variable,y=value,fill = variable)) + 
  geom_boxplot() + 
  facet_grid(true.covar~rho) +
  xlab("Time scale of fluctuation") +
  ylab("75th percentile of WMR") +
  theme_classic(base_size=14) +
  scale_fill_grey("Time scale",start = 0.3, end = 0.8) +
  geom_hline(lty=2,colour="red",yintercept = 0.85)
dev.off()



