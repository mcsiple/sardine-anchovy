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

# Load all the data
load("~/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/RAM_Barange_States.RData") # data frame: RB 
# TO DO: 
#         - Add FAO data to data set
#         - Check whether Barange or RAM have years filled in

# Subset data to use ------------------------------------------------------
variables = c("landings","ssb","rec")
regions = c("Benguela", "California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
dsources = c("RAM","Barange")
d = 2; # set datasource if needed
r=1;v=1; #set region and variable

sp2 <- data.frame()
true.df <- data.frame()


for(r in 1:length(regions)){
    for(v in 1:3){
    data.points <- subset(RB,datasource == dsources[d] & 
                            region == regions[r] & 
                            variable == variables[v])
    
    if(nrow(data.points)==0 | 
       length(unique(data.points$Sardine.est))==1 |
       length(unique(data.points$Anchovy.est))==1){print("one time series missing - STOP")}
    
    # Standardize data
    std_sardine <- as.numeric(scale(data.points$Sardine.est)) # scale so center is 0
    if(all(is.na(std_sardine))) std_sardine <- rep(NA, times=length(std_sardine)) # In case all values are the same
    mean(std_sardine)
    
    std_anchovy <- as.numeric(scale(data.points$Anchovy.est)) 
    if(all(is.na(std_anchovy))) std_sardine <- rep(NA, times=length(std_anchovy))
    mean(std_anchovy)
    
    # Plot an example of standardized time series
    pd <- data.frame(year = c(1:length(std_anchovy),1:length(std_sardine)),
                    sp= c(rep("Anchovy",times=length(std_anchovy)),
                          rep("Sardine",times=length(std_sardine))),
                    std.B = c(std_anchovy,std_sardine))
    
                # ggplot(pd,aes(x=year,y=std.B,colour=sp)) + 
                # geom_line(lwd=1.2) + 
                # scale_colour_manual(values=c("#ef8a62","#67a9cf")) +theme_classic()
                # dev.off()
    
    # Simulate time series of the same spectral characteristics as sardine and anchovy.
    sardine_phase <- surrogate(data.points$Sardine.est,method = 'phase')
    anchovy_phase <- surrogate(data.points$Anchovy.est,method = 'phase')

    # Plot surrogate time series
    pd <- data.frame(year = c(1:length(anchovy_phase),1:length(sardine_phase)),
    sp = c(rep("Anchovy",times=length(anchovy_phase)),
          rep("Sardine",times=length(sardine_phase))),
    std.B = c(as.numeric(scale(anchovy_phase)),as.numeric(scale(sardine_phase))))
    # pdf("EgPlot_Surrogates2.pdf",width=8,height=5,useDingbats = F)
     ggplot(pd,aes(x=year,y=std.B,colour=sp)) + 
       geom_line(lwd=1.2) + 
       scale_colour_manual(values=c("#ef8a62","#67a9cf")) +
       theme_classic()
    # dev.off()
    
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
    # This cutoff can be flexible... it can be something else; we need to figure out what the right cutoff is (see appendix)
    yr1 <- yr5 <- yr10 <- vector()
    wmr.thresh <- 0.5
    
    for(s in 1:nsims){
      yr1[s] <- length(which(synch.1[[s]]<wmr.thresh))/length(which(!is.na(synch.1[[s]])))
      yr5[s] <- length(which(synch.5[[s]]<wmr.thresh))/length(which(!is.na(synch.5[[s]])))
      yr10[s] <- length(which(synch.10[[s]]<wmr.thresh))/length(which(!is.na(synch.10[[s]])))
    }
    
    
    # Plotting
    # Estimates for probability of detecting WMR <threshold, for CA landings
    sp <- as.data.frame(rbind(
    quantile(yr1,probs=c(0.05,0.25,0.5,0.75,0.95)),
    quantile(yr5, probs=c(0.05,0.25,0.5,0.75,0.95)),
    quantile(yr10,probs=c(0.05,0.25,0.5,0.75,0.95))
    ))
    
    # For storing individual sims - these are values of WMR
    # This stores all the values from each simulation
    sp.all <- list("1yr" = do.call(rbind,synch.1),  
                   "5yr" = do.call(rbind,synch.5), 
                   # TO DO: Organize so that you can save each variable/region combo separately. O_O
                   "10" = do.call(rbind,synch.10)) 
    
    # Save 
    sp$variable <- sp.all <- variables[v]
    sp$region <- sp.all <- regions[r]
    sp$datasource <- sp.all <- dsources[d]
    sp$scale <- c("less.than.5","five.ten","ten.plus")
    sp2 <- rbind(sp2,sp)
    
     
     # Save simulations
     #sp.all2 <- rbind(sp.all2,sp.all)         # Eek! 
     
    
    
    # What is this density for the true variable? i.e, calculate WMR for the real time series. Is it different from the null model? I.e., is the density below wmr.thresh similar to the one you would expect from random time series?
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
    # TRIM BABY TRIM - use values from w/in cone of influence
    for(c in 1:ncol(mmat)){
      mmat[1:to.trim[c],c] <- NA
      mmat[nrow(mmat):(nrow(mmat)-to.trim[c]),c] <- NA
    }
    
    ind <- which(mr$y < 5)
    true[[1]] <- mmat[,ind] # subset to values at a period of <5 yrs
    true.vec[1] <- length(which(true[[1]]<wmr.thresh))/length(which(!is.na(true[[1]]))) # proportion of values at period<5 that are <0.7 (compensatory dynamics side)
    ind2 <- which(mr$y > 5 & mr$y < 10)
    true[[2]] <- mmat[,ind2] #
    true.vec[2] <- length(which(true[[2]]<wmr.thresh))/length(which(!is.na(true[[2]])))
    ind3 <- which(mr$y > 10)
    true[[3]] <- mmat[,ind3] #
    true.vec[3] <- length(which(true[[3]]<wmr.thresh))/length(which(!is.na(true[[3]])))
    
    tv <- data.frame(scale = c("less.than.5","five.ten","ten.plus"),
                     obs = true.vec,
                     variable = variables[v],
                      region = regions[r],
                      datasource = dsources[d])
    true.df <- rbind(true.df,tv)
    
    print(sp2)
    print(true.df)
} # end variables loop
} # end regions loop 




list.results <- list()
list.results[[1]] <- sp2
list.results[[2]] <- sp2
list.results[[3]] <-sp2
list.results[[4]] <- sp2

df <- ldply (list.results, data.frame)
df <- subset(df, datasource=="Barange")
save(df,file = "NullModelDistributions_07cutoff_std.RData")
save(true.df,file = "TrueDF.RData")


# Start here if not running sims again! -----------------------------------
load("NullModelDistributions_07cutoff_std.RData")
load("TrueDF.RData")
true.df$ou <- NA
true.df$ou[which(true.df$obs>df$X95.)] <- "Asynchronous"
true.df$ou[which(true.df$obs<df$X5.)] <- "Synchronous"

df <- df %>% subset(variable != "rec")
true.df <- true.df %>% subset(variable !="rec")
true.df$ou[is.na(true.df$ou)] <- "not.significant" 
pali <- c(beyonce_palette(48)[3],"grey",
          beyonce_palette(48)[6])

pdf("ExpectationPlot_07cutoff_black.pdf",width=8,height=7,useDingbats = FALSE)
ggplot(df, aes(x=scale,y=X50.)) + 
  #geom_point(size=0.5) + 
  facet_grid(region~variable) + 
  scale_x_discrete(limits=c("ten.plus","five.ten","less.than.5"),labels=c("Long-term","Medium","Short")) +
  coord_flip() +
  geom_linerange(aes(x=scale, ymin = X5.,ymax=X95.),lwd=1.1,colour='darkgrey') +
  geom_linerange(aes(x=scale,ymin=X50.,ymax=X75.),lwd=2.5,colour="darkgrey") +
  #theme_classic(base_size=14) +
  geom_point(data = true.df,
             aes(x=scale,y=obs,colour=ou),size=4) +
  scale_colour_manual("",labels = c("Asynchronous","not significant","Synchronous"),values=pali) +
  theme(strip.background = element_blank()) +
  ylab("Degree of asynchrony") +
  xlab("Time scale") +
  theme_black(base_size=14)
dev.off()

#dev.off()

# 