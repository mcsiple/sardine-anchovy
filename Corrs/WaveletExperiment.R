# install.packages("mvcwt")
library(mvcwt)

# Example from vignette ---------------------------------------------------
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
# Five LMEs: NE Atlantic, Benguela, California, Kuroshio-Oyashio, Humboldt
load("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/R files/allsardineanchovy.RData")
source("/Users/mcsiple/Dropbox/Chapter3-SardineAnchovy/Code_SA/sardine-anchovy/Corrs/getMARSSstates.R")
data.list <- list()
region.list <- list()
regions <- unique(alldat$region)
dsources <- unique(alldat$datasource)
d=2 # Start with 
#for (d in 1:length(dsources)){
 for(r in 1:length(regions)){
# Need to fix FAO data before using this function on them.
      var.list <- list()
      region = regions[r]
      data_source = dsources[d]
      variables = c("rec","ssb","landings")
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
save(Barange,file = "BarangeStates.RData")
save(RAM, file = "RAMStates.RData")



# Start here if you don’t need to re-process anything! --------------------
load("BarangeStates.RData")
load("RAMStates.RData")


      data.points <- subset(RAM, region == "California" & variable == "ssb")
      # RAM only has both sardine and anchovy for any of the variables (e.g.,ss, rec) for Humboldt. 
      
      x = data.points$Year       # Time series sampling points (day of year in this case)
      y = data.points[,c("Sardine.est","Anchovy.est")] # Time series only (each column = one species)
    
      w = mvcwt(x, y, min.scale = 1, max.scale = 20)
      mr = wmr(w)
      #mr.boot = wmr.boot(w, smoothing = 1, reps = 1000, mr.func = "wmr")
      image(mr, reset.par = FALSE)
      contour(mr, bound = NA, add = TRUE)
      

      # If using simulated data (but need starting values)
      subyrs <- 1:100 #1:nrow(response.ts)
      x <- subyrs 
      y <- response.ts[subyrs,1:2]
      w = mvcwt(x, y, min.scale = 1, max.scale = 20)
      mr = wmr(w)
      
      # Plots
      par(mfrow=c(4,1),mar = c(2,2,1,3))
      plot(1:nrow(y),y[,1],type='l',ylim=c(-2,max(c(y[,1],y[,2]))))
      lines(1:nrow(y),y[,2],col='red')
      
      tx <- 0.2
      # Scale: <5 yr    
      mr$z[mr$z>0.95] <- NA
      ind <- which(mr$y < 5)
      trim.z <- mr$z[nrow(mr$z)-ind,,1]
      hist(trim.z,xlim=c(0,1),col="lightblue",probability = T,main='',xaxt='n') #
      text(tx,1,"<5 yr")
      
      # Scale: 5-10 yr
      ind2 <- which(mr$y > 5 & mr$y < 10)
      trim.z2 <- mr$z[nrow(mr$z)-ind2,,1]
      hist(trim.z2,xlim=c(0,1),col="lightblue",probability = T,main='',xaxt='n')
      text(tx,1,"5-10 yr")
      
      #Scale: 10+ yr
      ind3 <- which(mr$y > 10)
      trim.z3 <- mr$z[nrow(mr$z)-ind3,,1]
      hist(trim.z3,xlim=c(0,1),col="lightblue",probability = T,main='',xlab="Degree of synchrony",xaxt='n')
      text(tx,1,"10+ yr")
      axis(1,at=c(0,0.5,1.0), labels=c(0,0.5,1.0))
      
      # plot(response.ts[,1],type='l')
      # lines(response.ts[,2],type='l',col='red')
       
      