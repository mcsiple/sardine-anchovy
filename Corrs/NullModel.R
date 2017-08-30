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
library(mvcwt)
library(beyonce)
load("~/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/RAM_Barange_States.RData") # data frame: RB
figwd <- "/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Figures"
setwd(figwd)

# Subset data to use ------------------------------------------------------
variables = c("landings","ssb","rec")
regions = c("Benguela", "California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
dsources = c("RAM","Barange")

sp2 <- data.frame()
true.df <- data.frame()
sa.col <- c("red","darkblue")
ylabel <- c("Landings","SSB","Recruitment")

d = 2; r = 5
for(r in 1:length(regions)){
for(v in 1:3){
data.points <- subset(RB,datasource == dsources[d] & 
                        region == regions[r] & 
                        variable == variables[v])
if(nrow(data.points)==0 | length(unique(data.points$Sardine.est))==1 |
   length(unique(data.points$Anchovy.est))==1){print("one time series missing - STOP")}

#Standardize data
std_sardine <- as.numeric(scale(data.points$Sardine.est)) #data.points$Sardine.est - mean(data.points$Sardine.est)
if(all(is.na(std_sardine))) std_sardine <- rep(NA, times=length(std_sardine)) # In case all values are the same

std_anchovy <- as.numeric(scale(data.points$Anchovy.est)) #data.points$Anchovy.est - mean(data.points$Anchovy.est)
if(all(is.na(std_anchovy))) std_sardine <- rep(NA, times=length(std_anchovy))

# Simulate 100 time series of the same spectral characteristics as sardine and anchovy.
sardine_phase <- surrogate(data.points$Sardine.est,method = "phase")
#plot(1:length(sardine_phase),sardine_phase,type='l',ylim=range(c(sardine_phase,std_sardine)))
#lines(1:nrow(data.points),std_sardine,col='red')
anchovy_phase <- surrogate(data.points$Anchovy.est,method = 'phase')
#plot(1:length(anchovy_phase),anchovy_phase,type='l',ylim=range(c(anchovy_phase,std_anchovy)))
#lines(1:nrow(data.points),std_sardine,col='red')

# Generate surrogate time series from the data
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
  scale.exp = 0.5; nscales = get.nscales(x)
  min.scale = get.min.scale(x); max.scale = get.max.scale(x)
  scales = log2Bins(min.scale, max.scale, nscales)
  w = mvcwt(x, y)
  mr = wmr(w)
  #image(mr)
  # eliminate z values outside cone of influence
          to.trim <- round(scales) # this is the number of cells to trim from the matrix (top and bottom)
          mmat <- mr$z[,,1]
          # OMG THIS TOOK ME SO LONG
          for(c in 1:ncol(mmat)){
            mmat[1:to.trim[c],c] <- NA
            mmat[nrow(mmat):(nrow(mmat)-to.trim[c]),c] <- NA
          }
          
    ind <- which(mr$y < 5)
    synch.1[[s]] <- mmat[,ind] #

    ind2 <- which(mr$y > 5 & mr$y < 10)
    synch.5[[s]] <- mmat[,ind2] #

    ind3 <- which(mr$y > 10)
    synch.10[[s]] <- mmat[,ind3] #
  }

sa.range <- range(c(std_anchovy,std_sardine))
                   
      # plot(1:nyears,std_anchovy, type='l',col=sa.col[1],lwd=1.5,
      #      ylab=paste("Standardized",ylabel[v]),
      #      xlab="Year")
      # lines(1:nyears,std_sardine,type='l',col=sa.col[2],lwd=1.5)
      # legend("topleft",lty=c(1,1),lwd=c(1.5,1.5),legend = c("Anchovy","Sardine"),col=sa.col)

# On average, at each time scale, what is the density below 0.5 (asynchronous)?
# This cutoff can be flexible... it's not always necessarily exactly 0.5 (see appendix)
yr1 <- yr5 <- yr10 <- vector()

for(s in 1:nsims){
  yr1[s] <- length(which(synch.1[[s]]<0.7))/length(which(!is.na(synch.1[[s]])))
  yr5[s] <- length(which(synch.5[[s]]<0.7))/length(which(!is.na(synch.5[[s]])))
  yr10[s] <- length(which(synch.10[[s]]<0.7))/length(which(!is.na(synch.10[[s]])))
}

# Estimates for probability of detecting WMR <0.5, for CA landings
sp <- as.data.frame(rbind(
quantile(yr1,probs=c(0.05,0.25,0.5,0.75,0.95)),
quantile(yr5, probs=c(0.05,0.25,0.5,0.75,0.95)),
quantile(yr10,probs=c(0.05,0.25,0.5,0.75,0.95))
))

# For storing individual sims - these are values of WMR
sp.all <- list("1yr" = do.call(rbind,synch.1), # This stores all the values from each simulation 
               "5yr" = do.call(rbind,synch.5), # TO DO: Organize so that you can save each variable/region combo separately. O_O
               "10" = do.call(rbind,synch.10)) # Then they can be plotted with the distributions from the data (standardized sts)

# Save 
sp$variable <- sp.all <- variables[v]
sp$region <- sp.all <- regions[r]
sp$datasource <- sp.all <- dsources[d]
sp$scale <- c("less.than.5","five.ten","ten.plus")
 sp2 <- rbind(sp2,sp)

 
 # Save simulations
 sp.all2 <- rbind(sp.all2,sp.all)         # Eek! 
 
            # pal <- beyonce_palette(11)
            # plot(1:3,sp[,2],xaxt='n',ylim=c(0,1),pch=21,bg = pal[c(1,3,5)],ylab="Prob(WMR < 0.5)", xlab="")
            # axis(1, at = c(1,2,3), labels = c("<5 yr","5-10 yr","10+ yr"))
            # arrows(x0 = 1:3, x1 = 1:3,y0 = sp[,1], y1 = sp[,3],col = pal[c(1,3,5)],lwd = 1.5,length = 0.03,angle=90,code = 3)


# What is this density for the true variable? i.e, calculate WMR for the real time series. Is it different from the null model? I.e., is the density below 0.5 similar to the one you would expect from random time series?
# True wavelet form (get data to plot on the graph with the null model!)
xx <- 1:length(std_anchovy)
yy <- cbind(std_anchovy,std_sardine)
scale.exp = 0.5; nscales = get.nscales(x)
min.scale = get.min.scale(x); max.scale = get.max.scale(x)
scales = log2Bins(min.scale, max.scale, nscales)
w = mvcwt(xx, yy) # Defaults are the above vars
true <- list()    # z values for histograms of WMR
true.vec <- vector() # to store median densities < 0.5

mr = wmr(w)

to.trim <- round(scales) # this is the number of cells to trim from the matrix to eliminate things outside cone of influence (top and bottom)
mmat <- mr$z[,,1]
# TRIM BABY TRIM
for(c in 1:ncol(mmat)){
  mmat[1:to.trim[c],c] <- NA
  mmat[nrow(mmat):(nrow(mmat)-to.trim[c]),c] <- NA
}

ind <- which(mr$y < 5)
true[[1]] <- mmat[,ind] # subset to values at a period of <5 yrs
true.vec[1] <- length(which(true[[1]]<0.7))/length(which(!is.na(true[[1]]))) # proportion of values at period<5 that are <0.5 (compensatory dynamics side)
ind2 <- which(mr$y > 5 & mr$y < 10)
true[[2]] <- mmat[,ind2] #
true.vec[2] <- length(which(true[[2]]<0.7))/length(which(!is.na(true[[2]])))
ind3 <- which(mr$y > 10)
true[[3]] <- mmat[,ind3] #
true.vec[3] <- length(which(true[[3]]<0.7))/length(which(!is.na(true[[3]])))

tv <- data.frame(scale = c("less.than.5","five.ten","ten.plus"),
                 obs = true.vec,
                 variable = variables[v],
                  region = regions[r],
                  datasource = dsources[d])
true.df <- rbind(true.df,tv)

#points(1:3 + 0.02,true.vec,pch=21, bg = "orange")
#legend("topright",pch = c(21,21),pt.bg = c("grey","orange"),legend=c("null model (no synchrony)","Observation"))
print(sp2)
print(true.df)
}
}




list.results <- list()
#true2 <- true.df
list.results[[1]] <- sp2
list.results[[2]] <- sp2
list.results[[3]] <-sp2
list.results[[4]] <- sp2

df <- ldply (list.results, data.frame)
df <- subset(df, datasource=="Barange")
save(df,file = "NullModelDistributions_07cutoff_std.RData")

pdf("ExpectationPlot_07cutoff.pdf",width=8,height=7,useDingbats = FALSE)
ggplot(df, aes(x=scale,y=X50.)) + 
  #geom_point(size=0.5) + 
  facet_grid(region~variable) + 
  scale_x_discrete(limits=c("ten.plus","five.ten","less.than.5")) +
  coord_flip() +
  geom_linerange(aes(x=scale, ymin = X5.,ymax=X95.),lwd=0.7,colour='darkgrey') +
  geom_linerange(aes(x=scale,ymin=X50.,ymax=X75.),lwd=1.5,colour="darkgrey") +
  theme_classic(base_size=14) +
  geom_point(data = true.df,
             aes(x=scale,y=obs)) +
  theme(strip.background = element_blank()) +
  ylab("Degree of asynchrony") +
  xlab("Time scale")
dev.off()

#dev.off()

# 