# Prep data for measuring corelations and replacement
# This code allows the user to fill in missing values in sardine and anchovy time series with estimated states from MARSS (default) and long-term means. Having a value for every year (including years where no data was collected) is required for estimating WMR. 
# This is a simple way to fill in a few spaces, most datasets do not have a lot of NAs.
# Five LMEs: NE Atlantic, Benguela, California, Kuroshio-Oyashio, Humboldt

load(here::here("R/DataCleanup/allsardineanchovy_3.RData")) # data frame: alldat - this is the NEW (2019) data pile- including new time series from RAM and FAO

# Load functions
source(here::here("R/DataCleanup/getMARSSstates.R"))
source(here::here("R/DataCleanup/Fill_NAs_SA.R")) 

data.list <- list()
region.list <- list()
regions <- unique(alldat$region)
dsources <- unique(alldat$datasource) #Barange, RAM, FAO
variables <- c("rec","ssb","landings")

d=1 #janky but easier for troubleshooting if you do one data source at a time
for(r in 1:length(regions)){
  var.list <- list()
  region = regions[r]
  data_source = dsources[d]
  if(data_source=="FAO"){
      var.list[[1]] <- NA # for FAO only
      var.list[[2]] <- NA  # for FAO only
  }
  for(v in 1:3){
    var.list[[v]] <- getMARSSstates(data = alldat,
                                    region_or_subregion = region,
                                    scale = "Region", 
                                    data_source = data_source,
                                    variable = variables[v],
                                    ccf.calc = FALSE,
                                    MARSS.cov = FALSE,
                                    plot = FALSE,
                                    get.mean.instead = FALSE)
  }
  region.list[[r]] <- do.call(rbind,var.list)
}



if(d==1) Barange <- do.call(rbind,region.list)
if(d==2) RAM <- do.call(rbind,region.list)
if(d==3) {FAO <- do.call(rbind,region.list)
          FAO <- FAO[complete.cases(FAO),] }# remove the rows that are all NA}


RBF <- rbind(Barange,RAM,FAO)
save(RBF,file="RAM_Barange_FAO_States.RData")



# Fill in with lt means instead of states ---------------------------------------------------

data.list <- list()
region.list <- list()
regions <- unique(alldat$region)
dsources <- unique(alldat$datasource) #Barange, RAM, FAO
variables <- c("rec","ssb","landings")

d=3 # if d = 3 (data source == FAO), set v loop to 3:3 since there are only landings data for FAO
for(r in 1:length(regions)){
  var.list <- list()
  region = regions[r]
  data_source = dsources[d]
  if(data_source=="FAO"){
    var.list[[1]] <- NA # for FAO only
    var.list[[2]] <- NA  # for FAO only
  }
  for(v in 3:3){
    var.list[[v]] <- getMARSSstates(data = alldat,
                                    region_or_subregion = region,
                                    scale = "Region", 
                                    data_source = data_source,
                                    variable = variables[v],
                                    ccf.calc = FALSE,get.mean.instead = TRUE)
  }
  region.list[[r]] <- do.call(rbind,var.list)
}


if(d==1) Barange2<- do.call(rbind,region.list)
if(d==2) RAM2 <- do.call(rbind,region.list)
if(d==3) {FAO2 <- do.call(rbind,region.list)
          FAO2 <- FAO2[complete.cases(FAO2),] }

Barange2 <- do.call(rbind,region.list)
RAM2 <- do.call(rbind,region.list)
FAO2 <- do.call(rbind,region.list)
FAO2 <- FAO[complete.cases(FAO2),] # remove the rows that are all NA
RBF2 <- rbind(Barange,RAM,FAO)
save(RBF2,file="RAM_Barange_FAO_MeanFillIn.RData")
