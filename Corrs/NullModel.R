<<<<<<< HEAD
# Take original time series, compare wavelet modulus ratio distribution between null model (no synchrony)
# If you have two time series of the same spectral characteristics, randomly starting, how often do you observe asynchrony by accident at each time scale? This can be referred to as the "null model" - how much synchronicity would you see at random.


# load libraries and data -------------------------------------------------
library(fractal) # contains surrogate() fn
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(DescTools)
library(biwavelet)
library(mvcwt)

# Load all the data
#load("~/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/RAM_Barange_States.RData") # data frame: RB 
load("C:/Users/siplem/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/RAM_Barange_States.RData") # data frame: RB 


get_surrogates <- function(dat=RB, dsource, reg, var, nsurrogates){
  #' @description takes a pair of sardine-anchovy time series from the bigger dataset and generates surrogate time series that have the same time series properties but none of the phase information (i.e., no info about the relationship between the two time series).
  #' @param dat the dataset to extract from. NOTE: This time series should be continuous (i.e., no NAs!). Default "RB" which includes RAM and Barange data, already interpolated where necessary.
  #' @param dsource which dataset the data should come from. Eventually, choose between, FAO, RAM, Barange et al.
  #' @param reg region
  #' @param var variable - choose between (ssb, landings, rec)
  #' @param nsurrogates - number of surrogate time series to generate
  #' @return a list of time series info (region, datasource, variable) and a matrix of surrogates for each fish. The surrogates are organized as having rows=years, columns=surrogates
  data.points <- subset(dat,datasource == dsource & 
                          region == reg & 
                          variable == var)
  #print(data.points)
  if(nrow(data.points)==0 | 
     length(unique(data.points$Sardine))==1 |
     length(unique(data.points$Anchovy))==1){stop("Error: one time series missing")}
  
  # Standardize data
  std_sardine <- as.numeric(scale(data.points$Sardine)) # scale so center is 0
  if(all(is.na(std_sardine))) std_sardine <- rep(NA, times=length(std_sardine)) # In case all values are the same
  std_anchovy <- as.numeric(scale(data.points$Anchovy)) 
  if(all(is.na(std_anchovy))) std_sardine <- rep(NA, times=length(std_anchovy))
  
  nyears = length(std_anchovy)
  a.sims <- s.sims <- matrix(NA, nrow = nyears,ncol = nsurrogates)
  for(s in 1:nsurrogates){
    a.sims[,s] <- as.numeric(surrogate(data.points$Anchovy,method = 'phase'))
    s.sims[,s] <- as.numeric(surrogate(data.points$Sardine,method = 'phase'))
  }
  return(list(Region=reg,DataSource=dsource,Variable=var, Anchovy.surrogates = a.sims,Sardine.surrogates = s.sims))
}

( xx <- get_surrogates(dat = RB,dsource = "Barange",reg = "California",var = "ssb",nsurrogates = 10) )



get_wmr <- function(anchovy.ts,sardine.ts){ 
  #' @description takes surrogate or true time series of sardine and anchovy and calculates the wavelet modulus ratio
  #' @param anchovy.ts time series of anchovy abundance, can be observed or surrogate
  #' @param sardine.ts time series of sardine abundance, can be observed or surrogate
  #' @return a list of three matrices, each of which represents one timescale. The contents of each matrix are the wavelet modulus ratios
  x <- 1:length(anchovy.ts)
  y <- cbind(anchovy.ts,sardine.ts)
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
  
  # Synchrony at <5 year scale      
  ind <- which(mr$y < 5)
  synch.1 <- mmat[,ind]
  
  # Synchrony at 5-10 year scale
  ind2 <- which(mr$y > 5 & mr$y < 10)
  synch.5 <- mmat[,ind2]
  
  # Synchrony at 10+ year scale
  ind3 <- which(mr$y > 10)
  synch.10 <- mmat[,ind3]
  
  return(list(less.than.5 = synch.1,
              five.ten = synch.5,
              ten.plus = synch.10))
}

get_wmr(anchovy.ts=xx$Anchovy.surrogates[,1],sardine.ts=xx$Sardine.surrogates[,1])





=======
# Take original time series, compare wavelet modulus ratio distribution between null model (no synchrony)
# If you have two time series of the same spectral characteristics, randomly starting, how often do you observe asynchrony by accident at each time scale? This can be referred to as the "null model" - how much synchronicity would you see at random.


# load libraries and data -------------------------------------------------
library(fractal) # contains surrogate() fn
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(DescTools)
library(biwavelet)
library(mvcwt)

# Load all the data
load(here::here("ProcData/RAM_Barange_States.RData")) # data frame: RB 


get_surrogates <- function(dat=RB, dsource, reg, var, nsurrogates){
  #' @description takes a pair of sardine-anchovy time series from the bigger dataset and generates surrogate time series that have the same time series properties but none of the phase information (i.e., no info about the relationship between the two time series).
  #' @param dat the dataset to extract from. Default "RB" which includes RAM and Barange data, already interpolated where necessary.
  #' @param dsource which dataset the data should come from. Eventually, choose between, FAO, RAM, Barange et al.
  #' @param reg region
  #' @param var variable - choose between (ssb, landings, rec)
  #' @param nsurrogates - number of surrogate time series to generate
  #' @return a list of time series info (region, datasource, variable) and a matrix of surrogates for each fish. The surrogates are organized as having rows=years, columns=surrogates
  data.points <- subset(dat,datasource == dsource & 
                          region == reg & 
                          variable == var)
  #print(data.points)
  if(nrow(data.points)==0 | 
     length(unique(data.points$Sardine.est))==1 |
     length(unique(data.points$Anchovy.est))==1){stop("Error: one time series missing")}
  
  # Standardize data
  std_sardine <- as.numeric(scale(data.points$Sardine.est)) # scale so center is 0
  if(all(is.na(std_sardine))) std_sardine <- rep(NA, times=length(std_sardine)) # In case all values are the same
  std_anchovy <- as.numeric(scale(data.points$Anchovy.est)) 
  if(all(is.na(std_anchovy))) std_sardine <- rep(NA, times=length(std_anchovy))
  
  nyears = length(std_anchovy)
  a.sims <- s.sims <- matrix(NA, nrow = nyears,ncol = nsurrogates)
  for(s in 1:nsurrogates){
    a.sims[,s] <- as.numeric(surrogate(data.points$Anchovy.est,method = 'phase'))
    s.sims[,s] <- as.numeric(surrogate(data.points$Sardine.est,method = 'phase'))
  }
  return(list(Region=reg,DataSource=dsource,Variable=var, Anchovy.surrogates = a.sims,Sardine.surrogates = s.sims))
}

( xx <- get_surrogates(dat = RB,dsource = "Barange",reg = "California",var = "ssb",nsurrogates = 10) )



get_wmr <- function(anchovy.ts,sardine.ts){ 
  #' @description takes surrogate or true time series of sardine and anchovy and calculates the wavelet modulus ratio
  #' @param anchovy.ts time series of anchovy abundance, can be observed or surrogate
  #' @param sardine.ts time series of sardine abundance, can be observed or surrogate
  #' @return a list of three matrices, each of which represents one timescale. The contents of each matrix are the wavelet modulus ratios
  x <- 1:length(anchovy.ts)
  y <- cbind(anchovy.ts,sardine.ts)
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
  
  # Synchrony at <5 year scale      
  ind <- which(mr$y < 5)
  synch.1 <- mmat[,ind]
  
  # Synchrony at 5-10 year scale
  ind2 <- which(mr$y > 5 & mr$y < 10)
  synch.5 <- mmat[,ind2]
  
  # Synchrony at 10+ year scale
  ind3 <- which(mr$y > 10)
  synch.10 <- mmat[,ind3]
  
  return(list(less.than.5 = synch.1,
              five.ten = synch.5,
              ten.plus = synch.10))
}

m.null = get_wmr(anchovy.ts=xx$Anchovy.surrogates[,1],sardine.ts=xx$Sardine.surrogates[,1])


test_wmr <- function(dat=RB, dsource, reg, var, m.null){
  data.points <- subset(dat,datasource == dsource & 
                          region == reg & 
                          variable == var)
  if(nrow(data.points)==0 | 
     length(unique(data.points$Sardine.est))==1 |
     length(unique(data.points$Anchovy.est))==1){stop("Error: one time series missing")}
  
  # Standardize data
  std_sardine <- as.numeric(scale(data.points$Sardine.est)) # scale so center is 0
  if(all(is.na(std_sardine))) std_sardine <- rep(NA, times=length(std_sardine)) # In case all values are the same
  std_anchovy <- as.numeric(scale(data.points$Anchovy.est)) 
  if(all(is.na(std_anchovy))) std_sardine <- rep(NA, times=length(std_anchovy))
  
  # get observed wmr
  m <- get_wmr(std_anchovy, std_sardine)
  
  # compare medians of observed and null wmrs
  test.1 = wilcox.test(m$less.than.5, m.null$less.than.5)
  test.5 = wilcox.test(m$five.ten, m.null$five.ten)
  test.10 = wilcox.test(m$ten.plus, m.null$ten.plus)
  
  test.df = data.frame(Region = reg, DataSource = dsource, Variable = var,
                       period = c("less.than.5","five.ten","ten.plus"), 
                       test.stat = c(test.1$statistic,test.5$statistic,test.10$statistic), 
                       p.value = c(test.1$p.value,test.5$p.value,test.10$p.value))
  
  return(test.df)
}

test_wmr(dat = RB,dsource = "Barange",reg = "California",var = "ssb",m.null = m.null)
>>>>>>> 1a340293de53b92724703c71f7163d02ea244066
