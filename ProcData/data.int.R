#Function data.int()

#Takes a dataset (from GraphClick) and interpolates fishing mortality or recruitment for all the years during which Fm or rec was estimated, then prints/saves interpolated data


data.proc <- function(data,fish.or.rec){   # fish.rec determines whether function is for landings/mortality or SSB/recruitment
  years <- round(data[,3])
  vars <- data[,4]
  print(cbind(years, vars))
  # Need to make sure missing years are really missing - expand ts but save NAs
  allyrs <- min(years,na.rm=T):max(years,na.rm=T)
  var.new <- vector(length=length(allyrs))
  for(y in 1:length(allyrs)){
    if(allyrs[y] %in% years){
          yr <- which(years==allyrs[y])
          var.new[y] <- vars[yr]
    }else(var.new[y] <- NA)
    }
  new.data <- cbind(allyrs,var.new)
  if (fish.or.rec == "fish")  colnames(new.data) <- c("FYear","Mf")
  else (colnames(new.data) <- c("RYear","Recruitment"))
  return(new.data)
}

data.int <- function(data,fish.or.rec){ # May need to fix this one
  years <- round(data[,3])
  minyear <- floor(min(years,na.rm=T))
  maxyear <- floor(max(years,na.rm=T))
  if(FALSE %in% is.na(data[,4])) {ap <- approx(x=years,y=data[,4],xout=minyear:maxyear)}    #If there is no recruitment data, fill recruitment column with NAs.
  else { ap <- cbind(years,rep(NA,times=length(years))) }
  if (fish.or.rec == "fish")  names(ap) <- c("FYear","Mf")
  else (names(ap) <- c("RYear","Recruitment"))
  data <- c(data[,1:2],ap)
}