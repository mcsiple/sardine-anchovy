
# Take original time series, compare wavelet modulus ratio distribution between null model (no synchrony)
# If you have two time series of the same spectral characteristics, randomly starting, how often do you observe asynchrony by accident at each time scale? This can be referred to as the "null model" - how much synchronicity would you see at random.


# load libraries and data -------------------------------------------------
library(fractal) # contains surrogate() fn
library(ggplot2)
library(dplyr)
library(reshape2)
library(DescTools)
library(biwavelet)
library(mvcwt)

# Load all the data
#load("~/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/RAM_Barange_States.RData") # data frame: RB 
#load(here::here("ProcData/RAM_Barange_FAO_States.RData")) # data frame: RBF (RBF is MARSS states, RBF2 is replacing NAs with the mean)

# Take original time series, compare wavelet modulus ratio distribution between null model (no synchrony)
# If you have two time series of the same spectral characteristics, randomly starting, how often do you observe asynchrony by accident at each time scale? This can be referred to as the "null model" - how much synchronicity would you see at random.


get_obs <- function(dat=RBF, dsource, reg, var){
  #' @description a function that subsets the data to the variable and region that the user specifies and returns the standardized biomass, landings, etc. time series
  #' #' @param dat the dataset to extract from. Default "RB" which includes RAM and Barange data, already interpolated where necessary.
  #' @param dsource which dataset the data should come from. Eventually, choose between, FAO, RAM, Barange et al.
  #' @param reg region
  #' @param var variable - choose between (ssb, landings, rec)
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
  return(list(std_anchovy=std_anchovy,std_sardine=std_sardine))
}

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


get_surrogates <- function(obs, nsurrogates){
  #' @description takes a pair of sardine-anchovy time series from the bigger dataset and generates surrogate time series that have the same time series properties but none of the phase information (i.e., no info about the relationship between the two time series).
  #' @param obs a list of two vectors, which are standardized sardine (std_sardine) and anchovy (std_anchovy) time series
  #' @param nsurrogates - number of surrogate time series to generate
  #' @return a list of time series info (region, datasource, variable) and a matrix of surrogates for each fish. The surrogates are organized as having rows=years, columns=surrogates
  #' 
  std_anchovy <- obs$std_anchovy
  std_sardine <- obs$std_sardine
  
  nyears = length(std_anchovy) # both vectors should be the same length
  a.sims <- s.sims <- matrix(NA, nrow = nyears,ncol = nsurrogates)
  for(s in 1:nsurrogates){
    a.sims[,s] <- as.numeric(surrogate(std_anchovy,method = 'phase'))
    s.sims[,s] <- as.numeric(surrogate(std_sardine,method = 'phase'))
  }
  return(list(Anchovy.surrogates = a.sims,Sardine.surrogates = s.sims))
}

# ( xx <- get_surrogates(obs = get_obs(dat = RBF,dsource = "Barange",reg = "California",var = "ssb"),
#                        nsurrogates = 10) )
# m.null = get_wmr(anchovy.ts=xx$Anchovy.surrogates[,1],sardine.ts=xx$Sardine.surrogates[,1]) # Sometimes this gives a %dopar% error, but it is ok.
# =======
# 
get_large_null <- function(dat = RB,dsource = "Barange",reg = "California",var = "ssb",nsims){
  #generate as many surrogates as you need to get your full sims:
  yy <- get_surrogates(obs = get_obs(dat = dat,dsource = dsource,reg = reg,var = var),nsurrogates = nsims)
  # Combine multiple null runs to get a good null dist:
  null.combined <- get_wmr(anchovy.ts=yy$Anchovy.surrogates[,1],sardine.ts=yy$Sardine.surrogates[,1])
  for (i in 2:nsims){
    newlist <- get_wmr(anchovy.ts=yy$Anchovy.surrogates[,i],sardine.ts=yy$Sardine.surrogates[,i])
    null.combined <- mapply(c, null.combined, newlist, SIMPLIFY=FALSE)
  }
  return(null.combined)
}


test_wmr <- function(obs, null.combined){
  #' @param obs a list of two vectors, which are standardized sardine (std_sardine) and anchovy (std_anchovy) time series
  std_anchovy <- obs$std_anchovy
  std_sardine <- obs$std_sardine
  # get observed wmr
  m <- get_wmr(std_anchovy, std_sardine)
  
  # use the full large null to compare location of WMRs
  test.1 = wilcox.test(m$less.than.5, m.null$less.than.5, conf.int = TRUE) 
  
  # assign NAs to these values if not enough wmr data to do wilcoxon tests
  if(all(IsZero(m$five.ten))) {
    test.5 <- list(estimate=NA,conf.int=c(NA,NA),p.value=NA)
    print("time series too short to get WMR for 5-10 yr period...returning NA")}else{
      test.5 = wilcox.test(m$five.ten, m.null$five.ten, conf.int = TRUE)
    }
  
  if(all(IsZero(m$ten.plus))) {
    test.10 <- list(estimate=NA,conf.int=c(NA,NA),p.value=NA)
    print("time series too short to get WMR for 10+ yr period...returning NA")}else{
      test.10 = wilcox.test(m$ten.plus, m.null$ten.plus, conf.int = TRUE)
    }
  
  # get relevant test information, 
  # where U is Mann-Whitney test stat (equivalent to Wilcoxson W here): n of all pairs where y =< x
  # diff is the median of diff between all pairs and bounded by CI,
  # Z is a standard normal test statistic derived from U, based on normal approximation
q1 <- qnorm(test.1$p.value)
q5 <- ifelse(is.na(test.5$p.value),NA,qnorm(test.5$p.value))
q10 <- ifelse(is.na(test.10$p.value),NA,qnorm(test.10$p.value))
  
  test.df = data.frame(period = c("less.than.5","five.ten","ten.plus"), 
                       U = c(test.1$statistic,test.5$statistic,test.10$statistic), 
                       diff = c(test.1$estimate,test.5$estimate,test.10$estimate), 
                       CI.L = c(test.1$conf.int[1],test.5$conf.int[1],test.10$conf.int[1]), 
                       CI.U = c(test.1$conf.int[2],test.5$conf.int[2],test.10$conf.int[2]), 
                       Z = c(q1,q5,q10),
                       p.value = c(test.1$p.value,test.5$p.value,test.10$p.value),
                       N = c(length(m$less.than.5)+length(null.combined$less.than.5),
                             length(m$five.ten)+length(null.combined$five.ten),
                             length(m$ten.plus)+length(null.combined$ten.plus)),
                       n1n2 = c(length(m$less.than.5)*length(null.combined$less.than.5),
                             length(m$five.ten)*length(null.combined$five.ten),
                             length(m$ten.plus)*length(null.combined$ten.plus)))
  # effect sizes
  test.df$r_p = test.df$U / test.df$n1n2 # common language effect size: proportion of pairs where y =< x
  test.df$r_bi = 1 - (2*test.df$U / test.df$n1n2) # another method (Wendt (1972)), rank biserial correlation
  test.df$r_z = test.df$Z/sqrt(test.df$N) # rank biserial correlation, relying on normal approximation
  
  return(test.df)
}


test_wmr_sub <- function(obs, null.combined, n.factor = 1){
  #' @param obs a list of two vectors, which are standardized sardine (std_sardine) and anchovy (std_anchovy) time series
  std_anchovy <- obs$std_anchovy
  std_sardine <- obs$std_sardine
  # get observed wmr
  m <- get_wmr(std_anchovy, std_sardine)
  
  # obtain random samples of large null WMR distributions, with sample size = n.factor * length(observed wmr)
  m.null <- list("less.than.5" = NA, "five.ten" = NA, "ten.plus" = NA)
  m.null$less.than.5 <- sample(null.combined$less.than.5, size = n.factor * length(m$less.than.5))
  m.null$five.ten <- sample(null.combined$five.ten, size = n.factor * length(m$five.ten))
  m.null$ten.plus <- sample(null.combined$ten.plus, size = n.factor * length(m$ten.plus))
  # use this sampled null to compare medians of observed and null wmrs
  test.1 = wilcox.test(m$less.than.5, m.null$less.than.5, conf.int = TRUE) 
  if(all(IsZero(m$five.ten))) {
    test.5 <- list(estimate=NA,conf.int=c(NA,NA))
    print("time series too short to get WMR for 5-10 yr period...returning NA")}else{
      test.5 = wilcox.test(m$five.ten, m.null$five.ten, conf.int = TRUE)
      }
  if(all(IsZero(m$ten.plus))) {
    test.10 <- list(estimate=NA,conf.int=c(NA,NA))
    print("time series too short to get WMR for 10+ yr period...returning NA")}else{
      test.10 = wilcox.test(m$ten.plus, m.null$ten.plus, conf.int = TRUE)
    }
  
  # get relevant test information, where diff is the median of diff between samples and bounded by CI
  test.df = data.frame(period = c("less.than.5","five.ten","ten.plus"), 
                       test.stat = c(test.1$statistic,test.5$statistic,test.10$statistic), 
                       diff = c(test.1$estimate,test.5$estimate,test.10$estimate), 
                       CI.L = c(test.1$conf.int[1],test.5$conf.int[1],test.10$conf.int[1]), 
                       CI.U = c(test.1$conf.int[2],test.5$conf.int[2],test.10$conf.int[2]), 
                       p.value = c(test.1$p.value,test.5$p.value,test.10$p.value))
  return(list(test.df, m.null))
}
