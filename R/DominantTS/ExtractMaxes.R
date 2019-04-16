#library(plyr)
library(dplyr)

#   -----------------------------------------------------------------------

# This function takes variable, region, and returns the dominant sardine and anchovy in that ecosystem
# "Dominant" = has the max value for biomass, catches, etc.

extract.maxes <- function(data = alldat, region_or_subregion, scale = "Region", data_source, variable = "landings"){
  #' @param data = dataset including the time series you're interested in
  #' @param region_or_subregion name of the region or subregion that you want to extract max values for;
  #' either a subregion or 1 of 5 LMEs (Benguela, California, NE Atlantic, Kuroshio-Oyashio, Humboldt)
  #' @param scale = whether the name in "region_or_subregion" is the name of a full LME or a subregion
  #' @param data_source = FAO, RAM, or Barange
  #' @param variable = the variable (rec, biomass, or landings)
  
  if (scale == "Region") {dataset <- filter(data, region == region_or_subregion & datasource==data_source)}
  if (scale == "Subregion") {dataset <- filter(data, subregion == region_or_subregion & datasource==data_source) }   
  
  if(variable=="landings"){
    lt.maxes <- dataset %>% 
      group_by(sp,stock) %>% 
      summarize(max.var=round(max(landings,na.rm=TRUE),2))}
   
  if(variable=="ssb"){
    lt.maxes <- dataset %>% 
      group_by(sp,stock) %>% 
      summarize(max.var=round(max(ssb,na.rm=TRUE),2))}
      
  if(variable=="rec"){
    lt.maxes <- dataset %>% 
      group_by(sp,stock) %>% 
      summarize(max.var=round(max(rec,na.rm=TRUE),2))}
  if(variable=="fishing.mortality"){ # Max (F) doesn't make sense, but including it anyway to be thorough
    lt.maxes <- dataset %>% 
      group_by(sp,stock) %>% 
      summarize(max.var=round(max(fishing.mortality,na.rm=TRUE),2))}
  #print(lt.maxes)
  
  #anchovy stats
  lt.max.a <- max(lt.maxes[which(lt.maxes$sp=="Anchovy"),ncol(lt.maxes)],na.rm=TRUE)   
  lt.max.sp <- as.character(lt.maxes[lt.maxes$max.var==lt.max.a,2])    # Which anchovy species had biggest long term value for this time series (i.e., the "dominant anchovy species")
  
  #sardine stats
  lt.max.s <- max(lt.maxes[which(lt.maxes$sp=="Sardine"),ncol(lt.maxes)],na.rm=TRUE)
  lt.max.sp.sar <- as.character(lt.maxes[lt.maxes$max.var==lt.max.s,2])     # "Dominant sardine species"
  
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
 
  
  return(list(Dom_anch_LTmax=as.character(lt.max.sp),
              Dom_sard_LTmax=as.character(lt.max.sp.sar),
              dom.anch = dom.a.ts,dom.sard = dom.s.ts,
              max.table = as.data.frame(lt.maxes)))
 
}  #End of extract.maxes fxn

#extract.maxes(data = alldat,region_or_subregion = "Benguela",data_source = "Barange" )

