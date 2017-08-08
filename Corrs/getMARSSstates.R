
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(MARSS)

getMARSSstates <- function(data = alldat, region_or_subregion = "Benguela", scale = "Region", data_source = "Barange", variable = "landings",MARSS.cov = FALSE, plot = FALSE, ccf.calc = FALSE){
  # param  data = dataset including the time series you're interested in
  # Region = 1 of 5 LMEs (Benguela, California, NE Atlantic, Kuroshio-Oyashio, Humboldt)
  # variable = the variable (rec, biomass, or landings) 
  # Datasource = FAO, RAM, or Barange
  
  if (scale == "Region") {dataset <- filter(data, region == region_or_subregion & datasource==data_source)}
  if (scale == "Subregion") {dataset <- filter(data, subregion == region_or_subregion & datasource==data_source)}   
  
  #dataset <- subset(dataset,datasource == datasource)
  #print(dataset)
  # special case for FAO, which doesn't have stock names, only land areas:
  if(data_source == "FAO") {dataset$stock <- paste(dataset$scientificname, dataset$region)
                            dataset <- dataset %>% subset(!is.na(landings))}
  
  if(variable=="landings"){ 
    lt.maxes <- ddply(.data=dataset,.(sp,stock),summarize,max.var=round(max(landings,na.rm=TRUE),2))
    } 
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
    source("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/Fill_NAs_SA.R")
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
  
  # This is where this function diverges from the cor() function; return MARSS states
  source("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/Fill_NAs_SA.R")
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
  print(kem.sa,what = "states")
  sard.estimates <- exp(kem.sa$states[1,])
  anch.estimates <- exp(kem.sa$states[2,])
  
  output <- data.frame(Year = sard.mars[1,],Sardine.est = sard.estimates, Anchovy.est = anch.estimates,region = region_or_subregion, datasource = data_source, variable = variable)
  
  return(output)
}  #End getMARSSstates function


# 
# plot(output$Year,output$Sardine.est,type='l')
# lines(output$Year,output$Anchovy.est,col="red")
# x = output$Year
# y = output[,-1]
# w = mvcwt(x, y, min.scale = 0.25, max.scale = 4)
# mr = wmr(w)
# image(mr, reset.par = FALSE)
# contour(mr, bound = NA, add = TRUE)

