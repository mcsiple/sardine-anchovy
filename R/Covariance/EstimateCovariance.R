# Get MARSS covariances and CIs for all the stocks and variables ----------

library(tidyverse)
source(here::here("R/DataCleanup/getMARSSstates.R"))
source(here::here("R/DataCleanup/Fill_NAs_SA.R"))
load(here::here("R/DataCleanup/allsardineanchovy_3.RData")) # dataframe all

region <- c("Benguela","California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
variable <- c("rec","ssb","landings")
datasource <- c("Barange")

covs <- tidyr::crossing(region,variable, datasource)
covs$loCI <- covs$med <-  covs$hiCI <- NA

for (i in 1:nrow(covs)){ 
  output <- getMARSSstates(data = alldat,region_or_subregion = covs$region[i],scale = "Region",data_source = covs$datasource[i],variable = covs$variable[i],ccf.calc=FALSE,get.mean.instead = TRUE,MARSS.cov = T)
  # covs$med[i] <- output$Q12
  # covs$loCI[i] <- output$loQ12
  # covs$hiCI[i] <- output$hiQ12
  covs$b1.sard = output$b1.sard
  covs$b1.sard.lo = output$b1.sard.lo
  covs$b1.sard.hi = output$b1.sard.hi
  covs$b2.anch = output$b2.anch
  covs$b2.anch.lo = output$b2.anch.lo
  covs$b2.anch.hi = output$b2.anch.hi
  # Interaction
  covs$B.12 = output$B.12
  covs$lo.B12 = output$lo.B12
  covs$hi.B12 = output$hi.B12
  #Covariance
  covs$Q1 = output$Q1
  covs$lo.Q1 = output$lo.Q1
  covs$ hi.Q1 = output$hi.Q1
  covs$Q2 = output$Q2
  covs$lo.Q2 = output$lo.Q2
  covs$hi.Q2 = output$hi.Q2
  print(covs)
}

#save(covs,file="MARSScovsBarange_May13_2019.RData")
load(file="MARSScovsBarange_May13_2019.RData")

covs2 <- covs %>% mutate(variable = fct_recode(variable, 
                                               "Spawning stock biomass" = "ssb",
                                               "Landings"="landings",
                                               "Recruitment"="rec"))

covs2$variable_f = factor(covs2$variable, levels=c('Recruitment',"Spawning stock biomass","Landings"))

(marsscovplot <- ggplot(covs2,aes(x=region,y=med)) + 
    facet_wrap(~variable_f) + 
    geom_point(size=2.5) +
    geom_linerange(aes(ymin=loCI, ymax=hiCI),lwd=1) +
    theme_classic(base_size = 14) +
    theme(strip.background = element_blank()) + #take out boxes around strips
    xlab("Region") +
    ylab("Covariance") +
    geom_hline(yintercept = 0,lty=2) + coord_flip() )

tiff("MARSSCovs_v2.tiff",width = 8,height = 2.5,units = 'in',res = 300)
marsscovplot
dev.off()

pdf("MARSSCovs_v2.pdf",width = 8,height = 2.5,useDingbats = FALSE)
marsscovplot
dev.off()
