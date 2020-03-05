# Get MARSS covariances and CIs for all the stocks and variables ----------

library(tidyverse)
library(here)
library(reshape2)
library(MARSS)

source(here("R/DataCleanup/Fill_NAs_SA.R"))
load(here("R/DataCleanup/allsardineanchovy_3.RData")) # dataframe alldat


# First, simplify data so it's only the "dominant" S and A (highest median biomass)
allfishies <- alldat %>% 
  melt(id.vars=c("datasource",
                 "scientificname",
                 "stock",
                 "year",
                 "sp",
                 "region",
                 "subregion"))

dom.by.ssb <- alldat %>%
  group_by(datasource,region,sp,stock) %>%
  summarize(median.ssb = median(ssb,na.rm=T),
            median.landings = median(landings,na.rm=T),
            median.tc = median(totalcatch,na.rm=T)) %>%
  mutate(is.dominant = ifelse(median.ssb == max(median.ssb,na.rm = T),1,0)) %>%
  filter(is.dominant == 1) %>%
  as.data.frame()

# Dominant stocks by median biomass are the same as the dominant ones by max biomass
domfishies <- allfishies %>%
  group_by(region,datasource,variable,sp) %>%
  filter(value==max(value,na.rm=T)) %>%
  mutate(is.dominant = 1) %>%
  select(-value,-year) %>%
  as.data.frame()

domstocks <- allfishies %>%
  group_by(region,datasource,variable,sp,stock) %>%
  summarize(med.value = median(value,na.rm=T)) %>%
  group_by(region,datasource,variable,sp) %>%
  filter(med.value == max(med.value,na.rm=T)) %>%
  mutate(is.dominant = 1) %>%
  select(-med.value) %>%
  right_join(allfishies) %>%
  filter(is.dominant == 1)

domstocks %>% 
  filter(datasource=="Barange") %>% 
  ggplot(aes(x=year,y=value,colour=sp)) + 
  geom_line() + 
  facet_grid(variable~region,scales="free_y")

regions <- c("Benguela","California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
nregions <- length(regions)
variables <- c("landings","ssb","rec")
nvariables <- length(variables)

dataset <- domstocks %>% filter(datasource=="Barange") # Fig 3A based on Barange et al. 2009 ts
maxyr <- max(dataset$year)
minyr <- min(dataset$year)

for(j in 1:length(variables)){
  # Set up data:
  marss.obj <- vector()
  v <- variables[j] 
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
  
  B.vals <- list( # interspecific effects and autocorrelation are the same in all locations - this one has the lowest AICc among all variations for structure of B.
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
  B = matrix(B.vals,nrow=nregions*2,ncol=nregions*2,byrow=T)
  
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
  
  fit <-  MARSS(marss.obj, model=marss.model, control=list(maxit=1e4,trace=-1))
  
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

# This is table s4
all <- bind_rows(q.all,b.all) %>%
  mutate(Var_Cov = as.factor(Var_Cov),
         Region = as.factor(Region))

save(all, file = "BQestimates.RData")
load("BQestimates.RData")

table.s4 <- all %>%
  select(variable,Region,param,med,loCI,hiCI,Var_Cov) %>%
  mutate_at(vars(med, loCI,hiCI), list(~ round(., 2)))

write.csv(table.s4,"TableS4.csv")

# Figure 
covar.plot <- q.all %>% 
  filter(Var_Cov=="cov") %>%
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
# One B per region
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

