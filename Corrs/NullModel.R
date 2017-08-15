# Take original time series, compare wavelet modulus ratio distribution between null model (no synchrony)
# If you have two time series of the same spectral characteristics, randomly starting, how often do you observe asynchrony by accident at each time scale?
# This can then be thought of as a "null model"


# load libraries and data -------------------------------------------------
library(fractal) # contains surrogate function
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(DescTools)
library(biwavelet)
load("~/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/RAM_Barange_States.RData") # data frame RB


# Subset data to use ------------------------------------------------------
variables = c("landings","ssb","rec")
regions = c("Benguela", "California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
dsources = c("RAM","Barange")
d = 2; r = 2; v = 1
data.points <- subset(RB,datasource == dsources[d] & 
                        region == regions[r] & 
                        variable == variables[v])
if(nrow(data.points)==0 | length(unique(data.points$Sardine.est))==1 |
   length(unique(data.points$Anchovy.est))==1){print("one time series missing - STOP")}

std_sardine <- data.points$Sardine.est - mean(data.points$Sardine.est)
std_anchovy <- data.points$Anchovy.est - mean(data.points$Anchovy.est)

# Simulate 100 time series of the same spectral characteristics as sardine and anchovy.
sardine_phase <- surrogate(data.points$Sardine.est,method = "phase")
#plot(1:length(sardine_phase),sardine_phase,type='l',ylim=range(c(sardine_phase,std_sardine)))
#lines(1:nrow(data.points),std_sardine,col='red')
anchovy_phase <- surrogate(data.points$Anchovy.est,method = 'phase')
#plot(1:length(anchovy_phase),anchovy_phase,type='l',ylim=range(c(anchovy_phase,std_anchovy)))
#lines(1:nrow(data.points),std_sardine,col='red')


nsims = 1000
nyears = length(std_anchovy)
a.sims <- s.sims <- matrix(NA, nrow = nyears,ncol = nsims)
  for(s in 1:nsims){
    a.sims[,s] <- as.numeric(surrogate(data.points$Anchovy.est,method = 'phase'))
    s.sims[,s] <- as.numeric(surrogate(data.points$Sardine.est,method = 'phase')) 
  }

synch.1 <- synch.5 <- synch.10 <- list()

for(s in 1:nsims){
  x <- 1:nrow(a.sims)
  y <- cbind(a.sims[,s],s.sims[,s])
  w = mvcwt(x, y, min.scale = 1, max.scale = 20)
  mr = wmr(w)
      ind <- which(mr$y < 5)
      synch.1[[s]] <- mr$z[nrow(mr$z)-ind,,1] #
      
      ind2 <- which(mr$y > 5 & mr$y < 10)
      synch.5[[s]] <- mr$z[nrow(mr$z)-ind2,,1] #
      
      ind3 <- which(mr$y > 10)
      synch.10[[s]] <- mr$z[nrow(mr$z)-ind3,,1] #
  }

sa.col <- c("#ef8a62","#67a9cf")

plot(1:nyears,a.sims[,1], type='l',col=sa.col[1],lwd=1.5)
lines(1:nyears,s.sims[,1],type='l',col=sa.col[2],lwd=1.5)

length(which(synch.1[[1]]<0.3))/length(synch.1[[1]])
# On average, at each time scale, what is the density below 0.5 (asynchronous)?
yr1 <- yr5 <- yr10 <- vector()

for(s in 1:nsims){
  yr1[s] <- length(which(synch.1[[s]]<0.5))/length(synch.1[[s]])
  yr5[s] <- length(which(synch.5[[s]]<0.5))/length(synch.5[[s]])
  yr10[s] <- length(which(synch.10[[s]]<0.5))/length(synch.10[[s]])
}

# Estimates for probability of detecting WMR <0.5, for CA landings
sp <- rbind(
quantile(yr1,probs=c(0.05,0.5,0.95)),
quantile(yr5, probs=c(0.05,0.5,0.95)),
quantile(yr10,probs=c(0.05,0.5,0.95))
)

pal <- beyonce_palette(11)
plot(1:3,sp[,2],xaxt='n',ylim=c(0,1),pch=21,bg = pal[c(1,3,5)])
axis(1, at = c(1,2,3), labels = c("<5 yr","5-10 yr","10+ yr"))
arrows(x0 = 1:3, x1 = 1:3,y0 = sp[,1], y1 = sp[,3],col = pal[c(1,3,5)],lwd = 1.5,length = 0.03,angle=90,code = 3)

# What is this density for the true variable? i.e, calculate WMR for the real time series. Is it different from the null model? I.e., is the density below 0.5 similar to the one you would expect from random time series?
# True wavelet form (get data to plot on the graph with the null model!)
x <- 1:nrow(a.sims)
y <- cbind(a.sims[,s],s.sims[,s])
w = mvcwt(x, y, min.scale = 1, max.scale = 20)
mr = wmr(w)
ind <- which(mr$y < 5)
synch.1[[s]] <- mr$z[nrow(mr$z)-ind,,1] #

ind2 <- which(mr$y > 5 & mr$y < 10)
synch.5[[s]] <- mr$z[nrow(mr$z)-ind2,,1] #

ind3 <- which(mr$y > 10)
synch.10[[s]] <- mr$z[nrow(mr$z)-ind3,,1] #



# Compare to one of the runs using a K-S test? 