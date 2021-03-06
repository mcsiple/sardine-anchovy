
library(dplyr)
library(ggplot2)
library(reshape2)
library(MARSS)

impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

impute.sa <- function(x,thresh=0.3){
  #' @description  This function can be applied to any time series. It imputes the mean for any NA-valued years after the start of the time series
  #' @param x is a vector of values with some missing data
  #' @param thresh is the threshold of how many NAs to accept in the time series (might not use)
  n <- length(x)
  if(all(is.na(x))){new <- x}else{
    firstvalue <- which(!is.na(x))[1] # first non-NA value in the time series
    # only apply the NA fill function after the start of the actual ts values
    newts <- x[firstvalue:n]
    imputed <- impute.mean(newts)
    new <- c(rep(NA,times=(firstvalue-1)),imputed)
  }
  return(new)
}


getMARSSstates <- function(data = alldat, region_or_subregion = "California", scale = "Region", data_source = "FAO", variable = "landings",MARSS.cov = FALSE, plot = FALSE, ccf.calc = FALSE,get.mean.instead = FALSE){
  #' @param data (dataframe) big dataframe of sardine and anchovy time series
  #' @param region_or_subregion (character) one of five regions (subregions are possible but not used for paper)
  #' @param scale (character) "Region" or "Subregion"
  #' @param data_source (character) where the data come from. "Barange" for all the covariance estimates because that is the most complete data source.
  #' @param variable "landings", "ssb", or "rec"
  #' @param MARSS.cov (logical) whether or not to estimate covariance
  #' @param plot (logical) whether to plot the ts
  #' @param ccf.calc (logical) whether to get cross correlation function (haven't used in a while)
  #' @param get.mean.instead (logical) whether to return time series interpolated with long-term mean (deprecated)
  
  if (scale == "Region") {dataset <- filter(data, region == region_or_subregion & datasource==data_source)}
  if (scale == "Subregion") {dataset <- filter(data, subregion == region_or_subregion & datasource==data_source)}   
  
  # special case for FAO, which doesn't have stock names, only land areas: 
  # if(data_source == "FAO") {dataset$stock <- paste(dataset$scientificname, dataset$region)
  #                           dataset <- dataset %>% subset(!is.na(landings))}
  
  if(data_source == "RAM"){ # RAM sometimes has total catch instead of landings-- we use it instead if that's all that's available
    dataset$landings[is.na(dataset$landings)] <- dataset$totalcatch
  }
  
  
  if(variable=="landings"){ 
    lt.maxes <- dataset %>% 
      group_by(sp,stock) %>% 
      summarize(max.var=round(max(landings,na.rm=T),2)) %>% 
      as.data.frame()
  } 
  if(variable=="ssb"){
    lt.maxes <- dataset %>% 
      group_by(sp,stock) %>% 
      summarize(max.var=round(max(ssb,na.rm=T),2)) %>% 
      as.data.frame()
  }
  if(variable=="rec"){
    lt.maxes <- dataset %>% 
      group_by(sp,stock) %>% 
      summarize(max.var=round(max(rec,na.rm=T),2)) %>% 
      as.data.frame()
  }
  if(variable=="fishing.mortality"){
    lt.maxes <- dataset %>% 
      group_by(sp,stock) %>% 
      summarize(max.var=round(max(fishing.mortality,na.rm=T),2)) %>% 
      as.data.frame()
  }
  #print(lt.maxes)
  
  #anchovy stats
  lt.max.a <- max(lt.maxes[which(lt.maxes$sp=="Anchovy"),ncol(lt.maxes)])   
  lt.max.sp <- lt.maxes[lt.maxes$max.var==lt.max.a,2]    #Which anchovy species had largest long term value for this time series (i.e., the "dominant anchovy species")
  
  #sardine stats
  lt.max.s <- max(lt.maxes[which(lt.maxes$sp=="Sardine"),ncol(lt.maxes)])
  lt.max.sp.sar <- lt.maxes[lt.maxes$max.var==lt.max.s,2]     # "Dominant sardine species"
  
  if(is.infinite(lt.max.s) & is.infinite(lt.max.a)){
    stop("both sardine and anchovy are NA")
  }
  
  if(is.infinite(lt.max.s)){
    dom.s.ts <- dataset[which(dataset$stock==lt.max.sp.sar[1]),] # just put in a random stock bc it doesn't matter
  }else{
    dom.s.ts <- dataset[which(dataset$stock==lt.max.sp.sar),]
  }
  
  if(is.infinite(lt.max.a)){
    dom.a.ts <- dataset[which(dataset$stock==lt.max.sp[1]),] # just put in a random stock bc it doesn't matter
  }else{
    dom.a.ts <- dataset[which(dataset$stock==lt.max.sp),]
  }
  
  
  if(variable=="landings"){
    sar = dom.s.ts$landings
    anch = dom.a.ts$landings}
  if(variable=="ssb"){
    sar = dom.s.ts$ssb
    anch = dom.a.ts$ssb}
  if(variable=="rec"){
    sar = dom.s.ts$rec
    anch = dom.a.ts$rec}
  if(variable=="fishing.mortality"){
    sar = dom.s.ts$fishing.mortality
    anch = dom.a.ts$fishing.mortality}
  #Plot the two dominant stocks/ts
  if (plot == TRUE){
    plot(dom.s.ts$year,sar,
         type="l",
         xlab="Year",ylab=variable,main=paste(c(region_or_subregion,variable)),
         xlim=c(min(dataset$year),max(dataset$year)),
         ylim=c(0,max(lt.maxes[,ncol(lt.maxes)])))
  }
  
  # One correlation method: use MARSS to find covariance --------------------------------------------
  if(MARSS.cov == TRUE){
    sard.mars <- FillNAs.ts(cbind(dom.s.ts$year,sar), # fills in missing years with NAs
                            startyear=min(c(dom.s.ts$year,dom.a.ts$year)),
                            endyear=max(c(dom.s.ts$year,dom.a.ts$year)))
    anch.mars <- FillNAs.ts(cbind(dom.a.ts$year,anch),
                            startyear=min(c(dom.s.ts$year,dom.a.ts$year)),
                            endyear=max(c(dom.s.ts$year,dom.a.ts$year)))
    
    # can't log-transform if there are true zeroes in the data, so replace zeroes with a very low value
    # sar and anch are rows, columns are years:
    sard.mars[2,] <- replace(sard.mars[2,], sard.mars[2,]==0, 0.0001)
    anch.mars[2,] <- replace(anch.mars[2,], anch.mars[2,]==0, 0.0001)
    
    # sar and anch are rows, columns are years
    MAR.obj <- log(rbind(sard.mars[2,],anch.mars[2,]))    #Landings are log transformed
    colnames(MAR.obj) <- sard.mars[1,]
    MAR.obj <- zscore(MAR.obj) # z-score the data! then we set u to zero - NEW!!!
    plot(1:ncol(MAR.obj),MAR.obj[1,],type='l',main=paste(region_or_subregion,variable))
    lines(1:ncol(MAR.obj),MAR.obj[2,],col="blue")
    
    model.sa=list()
    model.sa$R="diagonal and equal" # used to be "diagonal and equal"
    model.sa$U="zero"
    model.sa$A="zero" # need to double check this value
    
    # Option 1: ESTIMATE B FULLY
    B1=matrix(list('b1','b12','b21','b2'),2,2) # email
    model.sa$B=B1 # diagonal is the same and off-diags are the same but they are interpreted differently! Subtract 1 from diag to get effect of species on itself (i.e., small diagonal of B means more density dependence). If species are fully density-independent, B_diag = 1) 
    model.sa$Q = "unconstrained"
    
    kem.sa = MARSS(MAR.obj, model=model.sa, control=list(maxit=1e4)) 
    correlation = kem.sa$par$Q[2]/(sqrt(kem.sa$par$Q[3]) * sqrt(kem.sa$par$Q[1]))
    
    Q1 <- as.numeric(kem.sa.CIs$par$Q)[1]
    lo.Q1 <- kem.sa.CIs$par.lowCI$Q[1]
    hi.Q1 <- kem.sa.CIs$par.upCI$Q[1]
    
    Q2 <- as.numeric(kem.sa.CIs$par$Q)[3]
    lo.Q2 <- kem.sa.CIs$par.lowCI$Q[3]
    hi.Q2 <- kem.sa.CIs$par.upCI$Q[3]
    #Covariance
    q12 <- as.numeric(kem.sa.CIs$par$Q[2])
    lo.q12 <- kem.sa.CIs$par.lowCI$Q[2]
    hi.q12 <- kem.sa.CIs$par.upCI$Q[2]
    
    R <- as.numeric(kem.sa.CIs$par$R)
    return(list(
      #Variance
      Q1 = Q1,
      lo.Q1 = lo.Q1,
      hi.Q1 = hi.Q1,
      Q2 = Q2,
      lo.Q2 = lo.Q2,
      hi.Q2 = hi.Q2,
      
      # Option 2 addition (basically just adding covariance estimate):
      q12 =q12,
      lo.q12=lo.q12,
      hi.q12=hi.q12,
      R = R))
  }else correlation = "no MARSS correlation calculated"
  
  if(ccf.calc==TRUE){
    ccf(x=sar, y=anch, na.action = na.pass,type = "correlation")
  }
  
  # This is where this function diverges from the cor() function; return MARSS states OR filled-in w means
  if(get.mean.instead==TRUE){
    sard.mars <- FillNAs.ts(cbind(dom.s.ts$year,sar),
                            startyear=min(c(dom.s.ts$year,dom.a.ts$year)),
                            endyear=max(c(dom.s.ts$year,dom.a.ts$year)))
    anch.mars <- FillNAs.ts(cbind(dom.a.ts$year,anch),
                            startyear=min(c(dom.s.ts$year,dom.a.ts$year)),
                            endyear=max(c(dom.s.ts$year,dom.a.ts$year)))
    sard.mars[2,] <- impute.sa(x = sard.mars[2,])
    anch.mars[2,] <- impute.sa(x = anch.mars[2,])
    
    output <- data.frame(Year = sard.mars[1,],Sardine.est = sard.mars[2,], Anchovy.est = anch.mars[2,],region = region_or_subregion, datasource = data_source, variable = variable)
    return(output)
    
  }else{
    sard.mars <- FillNAs.ts(cbind(dom.s.ts$year,sar),
                            startyear=min(c(dom.s.ts$year,dom.a.ts$year)),
                            endyear=max(c(dom.s.ts$year,dom.a.ts$year)))
    
    anch.mars <- FillNAs.ts(cbind(dom.a.ts$year,anch),
                            startyear=min(c(dom.s.ts$year,dom.a.ts$year)),
                            endyear=max(c(dom.s.ts$year,dom.a.ts$year)))
    
    # replace zeroes with a very low value
    # sar and anch are rows, columns are years:
    sard.mars[2,] <- replace(sard.mars[2,], sard.mars[2,]==0, 0.0001)
    anch.mars[2,] <- replace(anch.mars[2,], anch.mars[2,]==0, 0.0001)
    
    MAR.obj <- log(rbind(sard.mars[2,],anch.mars[2,]))    #Landings are log transformed
    colnames(MAR.obj) <- sard.mars[1,]
    
    model.sa=list()
    model.sa$Q="unconstrained"
    model.sa$R="diagonal and equal" # new from MDS
    model.sa$U="zero" # new from MDS
    B1=matrix(list('b1','b12','b21','b2'),2,2) 
    model.sa$B=B1
    
    kem.sa = MARSS(MAR.obj, model=model.sa, control=list(maxit=1000)) 
    print(kem.sa,what = "states")
    
    #back-transform values
    sard.estimates <- exp(kem.sa$states[1,])
    anch.estimates <- exp(kem.sa$states[2,])
    
    output <- data.frame(Year = sard.mars[1,],Sardine.est = sard.estimates, Anchovy.est = anch.estimates,region = region_or_subregion, datasource = data_source, variable = variable)
    
    return(output)}
}  #End getMARSSstates function


#Test
#getMARSSstates()


