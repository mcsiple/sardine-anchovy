
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(MARSS)

corr.fig <- function(data = alldat, region_or_subregion = "Benguela", scale = "Region", data_source = "Barange", variable = "landings",MARSS.cov = FALSE, plot = FALSE, ccf.calc = FALSE){
  # data = dataset including the time series you're interested in
  # Region = 1 of 5 LMEs (Benguela, California, NE Atlantic, Kuroshio-Oyashio, Humboldt)
  # variable = the variable (rec, biomass, or landings) 
  # Datasource = FAO, RAM, or Barange
  
  if (scale == "Region") {dataset <- filter(data, region == region_or_subregion & datasource==data_source)}
  if (scale == "Subregion") {dataset <- filter(data, subregion == region_or_subregion & datasource==data_source) }   
  
  #dataset <- subset(dataset,datasource == datasource)
  #print(dataset)
  if(variable=="landings"){
    lt.maxes <- ddply(.data=dataset,.(sp,stock),summarize,max.var=round(max(landings,na.rm=TRUE),2))} 
  if(variable=="ssb"){
    lt.maxes <- ddply(.data=dataset,.(sp,stock),summarize,max.var=round(max(ssb,na.rm=TRUE),2))}
  if(variable=="rec"){
    lt.maxes <- ddply(.data=dataset,.(sp,stock),summarize,max.var=round(max(rec,na.rm=TRUE),2))}
  if(variable=="fishing.mortality"){
    lt.maxes <- ddply(.data=dataset,.(sp,stock),summarize,max.var=round(max(fishing.mortality,na.rm=TRUE),2))}
  #print(lt.maxes)
  
  #anchovy stats
  lt.max.a <- max(lt.maxes[which(lt.maxes$sp=="Anchovy"),ncol(lt.maxes)])   
  lt.max.sp <- lt.maxes[lt.maxes$max.var==lt.max.a,2]    #Which anchovy species had biggest long term value for this time series (i.e., the "dominant anchovy species")
  
  #sardine stats
  lt.max.s <- max(lt.maxes[which(lt.maxes$sp=="Sardine"),ncol(lt.maxes)])
  lt.max.sp.sar <- lt.maxes[lt.maxes$max.var==lt.max.s,2]     # "Dominant sardine species"
  
  dom.s.ts <- dataset[which(dataset$stock==lt.max.sp.sar),]
  dom.a.ts <- dataset[which(dataset$stock==lt.max.sp),]
  #   print(head(dom.a.ts))
  #   print(head(dom.s.ts))
  
  if(variable=="landings"){sar = dom.s.ts$landings
  anch = dom.a.ts$landings} 
  if(variable=="ssb"){sar = dom.s.ts$ssb
  anch = dom.a.ts$ssb} 
  if(variable=="rec"){sar = dom.s.ts$rec
  anch = dom.a.ts$rec}
  if(variable=="fishing.mortality"){sar = dom.s.ts$fishing.mortality
  anch = dom.a.ts$fishing.mortality}
  #Plot the two dominant stocks/ts
  if (plot == TRUE){
    plot(dom.s.ts$year,sar,
         type="l",
         xlab="Year",ylab=variable,main=paste(c(region_or_subregion,variable)),
         xlim=c(min(dataset$year),max(dataset$year)),
         ylim=c(0,max(lt.maxes[,ncol(lt.maxes)])))
    #lines(dom.a.ts$year,anch,col='red')
  }
  
  # One correlation method: use MARSS to find covariance --------------------------------------------
  if(MARSS.cov == TRUE){
    source("Fill_NAs_SA.R")
    sard.mars <- FillNAs.ts(cbind(dom.s.ts$year,sar),
                            startyear=min(c(dom.s.ts$year,dom.a.ts$year)),
                            endyear=max(c(dom.s.ts$year,dom.a.ts$year)))
    anch.mars <- FillNAs.ts(cbind(dom.a.ts$year,anch),
                            startyear=min(c(dom.s.ts$year,dom.a.ts$year)),
                            endyear=max(c(dom.s.ts$year,dom.a.ts$year)))
    
    # sar and anch are rows, columns are years
    MAR.obj <- log(rbind(sard.mars[2,],anch.mars[2,]))    #Landings are log transformed
    colnames(MAR.obj) <- sard.mars[1,]
    
    model.sa=list()
    model.sa$Q="unconstrained"
    
    kem.sa = MARSS(MAR.obj, model=model.sa, control=list(maxit=1000)) 
    correlation = kem.sa$par$Q[2]/(sqrt(kem.sa$par$Q[3]) * sqrt(kem.sa$par$Q[1]))
    kem.sa.CIs = MARSSparamCIs(kem.sa,method="parametric",nboot=200)
    print(kem.sa.CIs)
  }
  
  else correlation = "no MARSS correlation calculated"
  
  if(ccf.calc==TRUE){
    ccf(x=sar, y=anch, na.action = na.pass,type = "correlation")
  }
  
  return(list(correlation = correlation,Dom_anch_LTmax=as.character(lt.max.sp),
              Dom_sard_LTmax=as.character(lt.max.sp.sar),dom.anch = dom.a.ts,dom.sard = dom.s.ts,max.table = lt.maxes))
  
  # print(kem.sa.CIs$par.upCI)
  # print(kem.sa.CIs$par.lowCI)
  #legend("topright",legend=c("sardine","anchovy"),lty=c(1,1),col=c('black','red'))
  #table(dataset$Land.Area,dataset$Scientific.name)
  
}  #End corr.fig function

