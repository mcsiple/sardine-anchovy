#Fills in NAs for missing years in a time series

FillNAs <- function(ObsMat,startyear,endyear){
  FillObs <- matrix(nrow=nrow(ObsMat),ncol=length(startyear:endyear))
  colnames(FillObs) <- startyear:endyear
  for (i in 1:ncol(FillObs)){
    if (colnames(FillObs)[i] %in% colnames(ObsMat) ){FillObs[,i] <- ObsMat[,which(colnames(ObsMat)==colnames(FillObs)[i])]}
  }
  return(FillObs)
}

#Fills NAs in a time series (for use w MARSS!)
FillNAs.ts <- function(ObsMat,startyear,endyear){ 
  #In ObsMat, the first column is years and the second is ts data 
  #startyear=min(ObsMat[,1])
  #endyear=max(ObsMat[,1])
  expanded.length <- length(startyear:endyear)
  Fill.ts <- cbind(startyear:endyear,rep(NA,times=expanded.length))
  
  for (i in 1:nrow(Fill.ts)){
    if(Fill.ts[i,1] %in% ObsMat[,1]){
      Fill.ts[i,2] <- ObsMat[which(ObsMat[,1]==Fill.ts[i,1]),2]
      }
  }
  
  return(t(Fill.ts))
  }

