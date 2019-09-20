# Another way to get at whether the time series are asycnhronous at short time scales is to use a DFA and fit a single trend. If sardine and anchovy are asynchronous, they should have opposite loadings at the time scale you're interested in. This starts the same as the "EstimateCovariance.R" code but uses DFA instead of MARSS. 
library(tidyverse)
library(reshape2)
library(MARSS)

source(here::here("R/DataCleanup/Fill_NAs_SA.R"))
load(here::here("R/DataCleanup/allsardineanchovy_3.RData")) # dataframe alldat




# First, simplify data so it's only the "dominant" sardine and anc --------
allfishies <- alldat %>% 
  melt(id.vars=c("datasource","scientificname","stock","year","sp","region","subregion"))

domfishies <- allfishies %>%
  group_by(region,datasource,variable,sp) %>% 
  filter(value==max(value,na.rm=T)) %>%
  mutate(dom_YN = 1) %>%
  select(-value,-year) %>%
  as.data.frame()

domstocks <- left_join(allfishies,domfishies) %>% 
  filter(dom_YN==1)

domstocks %>% 
  filter(datasource=="Barange") %>% 
  ggplot(aes(x=year,y=value,colour=sp)) + 
  geom_line() + 
  facet_grid(variable~region,scales="free_y")


dfa.stocks <- domstocks %>% filter(datasource=="Barange")

v = "ssb" # for testing

regions <- c("Benguela","California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
nregions=length(regions)
maxyr <- max(dfa.stocks$year)
minyr <- min(dfa.stocks$year)


sa.DFA <- function(dataset = dfa.stocks,v = "ssb"){
  #' @description takes data and returns the results of a DFA with a single trend
  #' @param dataset (dataframe) big dataframe of sardine and anchovy time series
  #' @param v "landings", "ssb", or "rec"

    
  # Set up data:
  dfa.obj <- vector()
  for(i in 1:nregions){
    dt <- subset(dataset,variable==v & region==regions[i])
    sard <- dt %>% filter(sp=="Sardine")
    anch <- dt %>% filter(sp=="Anchovy")
    sard.mars <- FillNAs.ts(ObsMat = cbind(sard$year,sard$value), # fills in missing years with NAs
                            startyear=minyr,
                            endyear=maxyr)
    anch.mars <- FillNAs.ts(cbind(anch$year,anch$value),
                            startyear=minyr,
                            endyear=maxyr)
    regioncombo <- rbind(sard.mars[2,],anch.mars[2,])
    print(regions[i])
    dfa.obj <- rbind(dfa.obj,regioncombo)
  }
  
  #dfa.obj <- rbind(minyr:maxyr,dfa.obj)
  #dfa.obj are the observations
  dfa.obj <- zscore(log(dfa.obj))
  
    dfa.model=list()
    Z.vals <- list(
      "z1s",0,0,0,0,
      "z1a",0,0,0,0,
      0,"z2s",0,0,0,
      0,"z2a",0,0,0,
      0,0,"z3s",0,0,
      0,0,"z3a",0,0,
      0,0,0,"z4s",0,
      0,0,0,"z4a",0,
      0,0,0,0,"z5s",
      0,0,0,0,"z5a")
    Z = matrix(Z.vals,nrow=nregions*2,ncol=nregions,byrow=T)
    R = "diagonal and equal" # obs error similar across species and ecosystems
    x0 = U = A = "zero"
    Q = "identity"
    
    
    
    dfa.model <- list(Z=Z,R=R,x0=x0,A=A,Q=Q)

    fit1 <-  MARSS(dfa.obj, model=dfa.model, control=list(maxit=1000,trace=-1))
    par(mfrow=c(1,5))
    for (i in 1:nregions){
      plot(minyr:maxyr,fit1$states[i,],type='l',xlab="Year",ylab="Trend")
    }
    
    return(MARSSparamCIs(fit1,method='hessian',nboot=1000))
}



sa.DFA()


