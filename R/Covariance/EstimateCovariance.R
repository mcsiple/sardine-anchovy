# Get MARSS covariances and CIs for all the stocks and variables ----------

library(tidyverse)
source(here::here("R/DataCleanup/getMARSSstates.R"))
source(here::here("R/DataCleanup/Fill_NAs_SA.R"))
load(here::here("R/DataCleanup/allsardineanchovy_3.RData")) # dataframe alldat

region <- c("Benguela","California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
variable <- c("rec","ssb","landings")
datasource <- c("Barange")

# #This is for testing only! Just for R1
# region <- c("Benguela","California","Humboldt")
# variable <- c("ssb","landings")
# datasource <- c("Barange")

covs <- tidyr::crossing(region,variable, datasource)
covs <- covs %>% add_column(b1.sard=NA,
                    b1.sard.lo=NA,
                    b1.sard.hi=NA,
                    b2.anch=NA,
                    b2.anch.lo=NA,
                    b2.anch.hi=NA,
                    B.12=NA,
                    lo.B12=NA,
                    hi.B12=NA,
                    Q1=NA,
                    lo.Q1=NA,
                    hi.Q1=NA,
                    Q2=NA,
                    lo.Q2=NA,
                    hi.Q2=NA,
                    q12=NA,
                    lo.q12=NA,
                    hi.q12=NA,
                    R=NA)

for (i in 1:nrow(covs)){ 
  output <- getMARSSstates(data = alldat,region_or_subregion = covs$region[i],scale = "Region",data_source = covs$datasource[i],variable = covs$variable[i],ccf.calc=FALSE,get.mean.instead = TRUE,MARSS.cov = T)
  # ONLY FOR OPTION 1
  # covs$b1.sard[i] = output$b1.sard
  # covs$b1.sard.lo[i] = output$b1.sard.lo
  # covs$b1.sard.hi[i] = output$b1.sard.hi
  # covs$b2.anch[i] = output$b2.anch
  # covs$b2.anch.lo[i] = output$b2.anch.lo
  # covs$b2.anch.hi[i] = output$b2.anch.hi
  # # Interaction
  # covs$B.12[i] = output$B.12
  # covs$lo.B12[i] = output$lo.B12
  # covs$hi.B12[i] = output$hi.B12
  #Covariance
  covs$Q1[i] = output$Q1
  covs$lo.Q1[i] = output$lo.Q1
  covs$ hi.Q1[i] = output$hi.Q1
  covs$Q2[i] = output$Q2
  covs$lo.Q2[i] = output$lo.Q2
  covs$hi.Q2[i] = output$hi.Q2
  
  # NEW just for option 2
  covs$q12[i] = output$q12
  covs$lo.q12[i] = output$lo.q12
  covs$hi.q12[i] = output$hi.q12
  
  covs$R[i] = output$R
  print(covs)
}

#save(covs,file="MARSScovsBarange_BCovTest.RData") # this one estimates B diagonals and one off-diag
#save(covs,file="MARSScovsBarange_QCovTest.RData")

# For Reviewer 1: should we estimate the B matrix instead of Q? -----------
# First, simplify data so it's only the "dominant" S and A 
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

regions <- c("Benguela","California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
nregions=length(regions)
maxyr <- max(dfa.stocks$year)
minyr <- min(dfa.stocks$year)

dataset = dfa.stocks
variables <- c("landings","ssb","rec")

for(j in 1:3){
    v = variables[j]
    # Set up data:
    marss.obj <- vector()
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
      marss.obj <- rbind(marss.obj,regioncombo)
    }
    
    
    marss.obj <- zscore(log(marss.obj))
    
    marss.model=list()
    
    Z = "identity"
    R = "diagonal and equal" # obs error similar across species and ecosystems
    x0 = U = A = "zero"
    #B = "unconstrained"
    # B.vals1 <- list( # B.vals2 assumes that the interspecific effects are the same in all locations, and has a lower AICc. The simplified one (fit3, below) has a lower AICc so we are going to use that one.
    #   "b1s","b1as",0,0,0,0,0,0,0,0,
    #   "b1sa","b1a",0,0,0,0,0,0,0,0,
    #   0,0,"b2s","b2as",0,0,0,0,0,0,
    #   0,0,"b2sa","b2a",0,0,0,0,0,0,
    #   0,0,0,0,"b3s","b3as",0,0,0,0,
    #   0,0,0,0,"b3sa","b3a",0,0,0,0,
    #   0,0,0,0,0,0,"b4s","b4as",0,0,
    #   0,0,0,0,0,0,"b4sa","b4a",0,0,
    #   0,0,0,0,0,0,0,0,"b5s","b5as",
    #   0,0,0,0,0,0,0,0,"b5sa","b5a")
    # B = matrix(B.vals1,nrow=nregions*2,ncol=nregions*2,byrow=T)
    
    # B.vals2 <- list( # could also estimate that the interspecific effects are the same in all locations
    #   "b1s","bas",0,0,0,0,0,0,0,0,
    #   "bsa","b1a",0,0,0,0,0,0,0,0,
    #   0,0,"b2s","bas",0,0,0,0,0,0,
    #   0,0,"bsa","b2a",0,0,0,0,0,0,
    #   0,0,0,0,"b3s","bas",0,0,0,0,
    #   0,0,0,0,"bsa","b3a",0,0,0,0,
    #   0,0,0,0,0,0,"b4s","bas",0,0,
    #   0,0,0,0,0,0,"bsa","b4a",0,0,
    #   0,0,0,0,0,0,0,0,"b5s","bas",
    #   0,0,0,0,0,0,0,0,"bsa","b5a")
    # B = matrix(B.vals2,nrow=nregions*2,ncol=nregions*2,byrow=T)
    # 
    B.vals3 <- list( # could also estimate that the interspecific effects are the same in all locations
      "bs","bas",0,0,0,0,0,0,0,0,
      "bsa","ba",0,0,0,0,0,0,0,0,
      0,0,"bs","bas",0,0,0,0,0,0,
      0,0,"bsa","ba",0,0,0,0,0,0,
      0,0,0,0,"bs","bas",0,0,0,0,
      0,0,0,0,"bsa","ba",0,0,0,0,
      0,0,0,0,0,0,"bs","bas",0,0,
      0,0,0,0,0,0,"bsa","ba",0,0,
      0,0,0,0,0,0,0,0,"bs","bas",
      0,0,0,0,0,0,0,0,"bsa","ba")
    B = matrix(B.vals3,nrow=nregions*2,ncol=nregions*2,byrow=T) 
    
    Q.vals <- list(
      "q1s","q1",0,0,0,0,0,0,0,0,
      "q1","q1a",0,0,0,0,0,0,0,0,
      0,0,"q2s","q2",0,0,0,0,0,0,
      0,0,"q2","q2a",0,0,0,0,0,0,
      0,0,0,0,"q3s","q3",0,0,0,0,
      0,0,0,0,"q3","q3a",0,0,0,0,
      0,0,0,0,0,0,"q4s","q4",0,0,
      0,0,0,0,0,0,"q4","q4a",0,0,
      0,0,0,0,0,0,0,0,"q5s","q5",
      0,0,0,0,0,0,0,0,"q5","q5a")
    Q = matrix(Q.vals,nrow=nregions*2,ncol=nregions*2,byrow=T)
    
    
    
    marss.model <- list(Z=Z,R=R,x0=x0,A=A,Q=Q,B=B,U=U)
    
    fit <-  MARSS(marss.obj, model=marss.model, control=list(maxit=1000,trace=-1))
    
    CIs <- MARSSparamCIs(fit,nboot=1000)
    
    Bestimates <- data.frame(param = c(rownames(CIs$par$B)),
                             med = CIs$par$B,
                             loCI = CIs$par.lowCI$B,
                             hiCI=CIs$par.upCI$B,
                             AICc = fit$AICc,
                             Obs.error = CIs$par$R,
                             lo.Obs.error = CIs$par.lowCI$R,
                             hi.Obs.error = CIs$par.upCI$R)
    Qestimates <- data.frame(param = rownames(CIs$par$Q),
                             med = CIs$par$Q,
                             loCI = CIs$par.lowCI$Q,
                             hiCI=CIs$par.upCI$Q,
                             AICc = fit$AICc,
                             Obs.error = CIs$par$R,
                             lo.Obs.error = CIs$par.lowCI$R,
                             hi.Obs.error = CIs$par.upCI$R)
    
    qp <- Qestimates %>% mutate(Var_Cov= ifelse(grepl(x = param, pattern = "s") | grepl(x = param, pattern = "a")  ,"var","cov"),
                                Region = ifelse(grepl(x=param,pattern="1"),"Benguela",
                                                ifelse(grepl(x=param,pattern="2"),"California",
                                                       ifelse(grepl(x=param,pattern="3"),"Humboldt",
                                                              ifelse(grepl(x=param,pattern="4"),"Kuroshio-Oyashio",
                                                                     "NE Atlantic")))))
    bp <- Bestimates %>% mutate(Var_Cov= ifelse(grepl(x=param,"as") | grepl(x=param,"sa")  ,"offdiag","diag"),
                                Region = "All")
    qp$Region = forcats::fct_rev(qp$Region) # reverse factor order
    bp$Region = forcats::fct_rev(bp$Region)
    
    if(v=="landings"){
      qp$variable = bp$variable = "Landings"
      q.landings = qp
      b.landings = bp
      
    }
    if(v=="ssb"){
      qp$variable = bp$variable = "Biomass"
      q.ssb = qp
      b.ssb = bp
    }
    if(v=="rec"){
      qp$variable = bp$variable = "Recruitment"
      q.rec = qp
      b.rec = bp
    }
} #end variables loop

q.all <- bind_rows(q.landings,q.ssb,q.rec)
q.all$variable <- forcats::fct_relevel(q.all$variable,"Landings", "Biomass", "Recruitment")

b.all <- bind_rows(b.landings,b.ssb,b.rec)
b.all$variable <- forcats::fct_relevel(b.all$variable,"Landings", "Biomass", "Recruitment")

covar.plot <- q.all %>% filter(Var_Cov=="cov") %>%
              ggplot(aes(x=Region,y=med)) + 
              geom_point(size=3) +
              geom_linerange(aes(ymin=loCI,ymax=hiCI),size=1) + 
              geom_hline(yintercept=0,lty=2) +
              xlab("Ecosystem") +
              ylab("Covariance") +
              coord_flip() +
              facet_wrap(~variable) +
              theme_classic(base_size=14) +
              theme(strip.background = element_blank())

pdf("Fig_4a_covariance_new.pdf",width = 9,height =3,useDingbats = F)
covar.plot
dev.off()


# Look at det(B)^2 to confirm the importance of species interaction --------
# One B per region - janky but functional
b.ssb <- b.all %>% filter(variable=="Biomass")
B1 <- matrix(b.all$med[1:4],nrow=2,ncol=2)
B2 <- matrix(b.all$med[c(5,2,3,6)],nrow=2,ncol=2)
B3 <- matrix(b.all$med[c(7,2,3,8)],nrow=2,ncol=2)
B4 <- matrix(b.all$med[c(9,2,3,10)],nrow=2,ncol=2)
B5 <- matrix(b.all$med[c(11,2,3,12)],nrow=2,ncol=2)
det(B1)^2
det(B2)^2
det(B3)^2
det(B4)^2
det(B5)^2