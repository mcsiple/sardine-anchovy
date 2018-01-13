# SKIP TO LINE 60 IF NO ADD'L DATA PROCESSING NEEDED
# install.packages("mvcwt")
library(mvcwt)


# Example from vignette ---------------------------------------------------
# With sample data from mvcwt
data(lrlake)
x = subset(lrlake, Basin == "Treatment", LRL.Day) / 365.25 # Time series sampling points (day of year in this case)
y = subset(lrlake, Basin == "Treatment", -(1:8)) # Time series only (each column = one species)
#y[20:30,9] <- NA #Test what happens where there are holes. Update: it is bad.
w = mvcwt(x, y, min.scale = 0.25, max.scale = 4)
mr = wmr(w)
image(mr, reset.par = FALSE)
contour(mr, bound = NA, add = TRUE)

plot(1:nrow(y),y[,1],type='l')
lines(1:nrow(y),y[,2],col='red')


# Try with MARSS states from sardine and anchovy time series --------------
# In order to get numbers for all years for patchy data sets, use MARSS to get states. 
# This is a simple way to fill in a few spaces, most datasets do not have a lot of NAs.
# Five LMEs: NE Atlantic, Benguela, California, Kuroshio-Oyashio, Humboldt
load("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/R files/allsardineanchovy.RData") # data frame: alldat
source("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/Corrs/getMARSSstates.R")
data.list <- list()
region.list <- list()
regions <- unique(alldat$region)
dsources <- unique(alldat$datasource)
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
RB <- rbind(Barange,RAM)
save(RB,file="RAM_Barange_States.RData")
# save(Barange,file = "BarangeStates.RData")
# save(RAM, file = "RAMStates.RData")


###########################################################################
###########################################################################
###########################################################################
# Start here if you don’t need to re-process anything! --------------------
###########################################################################
###########################################################################


# Test to double check that MARSS states didn’t lose any info (i.e --------
# If the model didn't converge, there might be bad estimates in the dataset RB,
# but we should still be able to use the data.
summary(alldat)
head(alldat)
m <- melt(alldat,id.vars=c("datasource","scientificname","year"))
ggplot(alldat,aes(x=year,y=ssb,colour=subregion)) + geom_point() + facet_wrap(~region)
#Need to be able to subset to just the dominant stock for each type/region **** TO DO


# Plot prefs
#devtools::install_github("dill/beyonce") # David Lawrence Miller's Beyoncé palettes
library(beyonce)
print(beyonce_palette(40))

var.colors <- data.frame(var = c("landings","ssb","rec"),
                         col = beyonce_palette(40)[c(1,4,7)])



# Load data etc.
load("~/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/ProcData/RAM_Barange_States.RData") # data frame RB

      variables = c("landings","ssb","rec")
      regions = c("Benguela", "California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
      dsources = c("RAM","Barange")
      

      #Cycle thru variables
      d <- 2 
      r <- 1
      v <- 2
      for(r in 1:length(regions)){
        par(mfcol=c(4,3),mar=c(2, 4, 1, 2) + 0.1)
      for(v in 1:3){
      var.color <- as.character(var.colors$col[v])
      tx <- 0.2 # Where text should be placed on histogram
      data.points <- subset(RB,datasource == dsources[d] & region == regions[r] & variable == variables[v])
      if(nrow(data.points)==0 | length(unique(data.points$Sardine.est))==1 |
         length(unique(data.points$Anchovy.est))==1) {plot.new()} else{
      # RAM only has both sardine and anchovy for any of the variables (e.g., ssb, rec) for Humboldt. 
      std_sardine <- as.numeric(scale(data.points$Sardine.est))  #Standardize: subtract mean, divide by stdev
      std_anchovy <- as.numeric(scale(data.points$Anchovy.est)) # Standardize
      x = data.points$Year       # Time series sampling points (day of year in this case)
      y = cbind(std_sardine,std_anchovy) # Time series only (each column = one species)
      
      # If using surrogate time series:
      #x = 1:length(anchovy_phase)
      #y = cbind(as.numeric(anchovy_phase),as.numeric(sardine_phase))
      
      # If using simulated data (NOTE: need starting values)
      scale.exp = 0.5; nscales = get.nscales(x)
      min.scale = get.min.scale(x); max.scale = get.max.scale(x)
      scales = log2Bins(min.scale, max.scale, nscales)
      w = mvcwt(x, y)
      mr = wmr(w)
      # eliminate z values outside cone of influence
      to.trim <- round(scales) # this is the number of cells to trim from the matrix (top and bottom)
      mmat <- mr$z[,,1]
      # OMG THIS TOOK ME SO LONG
      for(c in 1:ncol(mmat)){
        mmat[1:to.trim[c],c] <- NA
        mmat[nrow(mmat):(nrow(mmat)-to.trim[c]),c] <- NA
      }
    
      #image(mr, reset.par = FALSE)
      #contour(mr, bound = NA, add = TRUE)
      
      
      # Plots
      plot(1:nrow(y),y[,1],type='l',ylim=c(min(c(y[,1],y[,2])),max(c(y[,1],y[,2]))),
           xlab="Year",ylab="Standardized var",
           main=paste(regions[r]," - ",
                      variables[v]),col=sa.col[2],lwd=1.5)
      lines(1:nrow(y),y[,2],col=sa.col[1],lwd=1.5)
      
      mr$z[mr$z>0.95] <- NA
      
      # Scale: <5 yr    
      ind <- which(mr$y < 5)
      trim.z <- mmat[ind,]
      hist(trim.z,xlim=c(0,1),col=var.color,border=var.color,
           probability = T,main='',
           xaxt='n',xlab="",ylab="") 
      text(tx,1,"<5 yr")
      
      # Scale: 5-10 yr
      ind2 <- which(mr$y > 5 & mr$y < 10)
      trim.z2 <- mmat[ind2,]
      hist(trim.z2,xlim=c(0,1),col=var.color,border=var.color,
           probability = T,main='',
           xaxt='n',xlab="")
      text(tx,1,"5-10 yr")
      
      #Scale: 10+ yr
      ind3 <- which(mr$y > 10)
      trim.z3 <- mmat[ind3,]
      if(all(is.na(trim.z3))){plot.new()}else{
      hist(trim.z3,xlim=c(0,1),col=var.color,border=var.color,
           probability = T,main='',
           xlab="Degree of synchrony",
           xaxt='n',ylab="")
      text(tx,1,"10+ yr")
      axis(1,at=c(0,0.5,1.0), labels=c(0,0.5,1.0))}
      
      }
        } # end variables loop
        } #end regions loop
      # plot(response.ts[,1],type='l')
      # lines(response.ts[,2],type='l',col='red')
       
