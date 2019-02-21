# Prep data - long term mean vs. MARSS states

# Try with MARSS states from sardine and anchovy time series --------------
# In order to get numbers for all years for patchy data sets, use MARSS to get states. 
# This is a simple way to fill in a few spaces, most datasets do not have a lot of NAs.
# Five LMEs: NE Atlantic, Benguela, California, Kuroshio-Oyashio, Humboldt
load("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Datasets/allsardineanchovy_3.RData") # data frame: alldat - this is the NEW (2019) data pile
source("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/Corrs/getMARSSstates.R")
data.list <- list()
region.list <- list()
regions <- unique(alldat$region)
dsources <- unique(alldat$datasource) #Barange, RAM, FAO
variables <- c("rec","ssb","landings")

d=1
for(r in 1:length(regions)){
  var.list <- list()
  region = regions[r]
  data_source = dsources[d]
  for(v in 1:3){
    var.list[[v]] <- getMARSSstates(data = alldat,
                                    region_or_subregion = region,
                                    scale = "Region", 
                                    data_source = data_source,
                                    variable = variables[v],
                                    ccf.calc = FALSE)
  }
  region.list[[r]] <- do.call(rbind,var.list)
}



Barange <- do.call(rbind,region.list)
RAM <- do.call(rbind,region.list[1:3])
FAO <- do.call(rbind,region.list)
RB <- rbind(Barange,RAM)
save(RB,file="RAM_Barange_FAO_States.RData")
# save(Barange,file = "BarangeStates.RData")
# save(RAM, file = "RAMStates.RData")

