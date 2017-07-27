# install.packages("mvcwt")
library(mvcwt)

# Function takes sequences as columns and compute continuous wavelet transform of each 
    # mvcwt(x, y, scale.exp = 0.5, nscales = get.nscales(x),
    #       min.scale = get.min.scale(x), max.scale = get.max.scale(x),
    #       scales = log2Bins(min.scale, max.scale, nscales), loc = regularize(x),
    #       wave.fun = "Morlet")


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
unique(alldat$region)
# Need to fix FAO data before using this function on them.
      var.list <- list()
      region = "Benguela" # For Humboldt: warning for ssb and landings, not recruitment
      data_source = "RAM"
      variables = c("rec","ssb","landings")
      for(v in 1:3){
      var.list[[v]] <- getMARSSstates(data = alldat,
                                 region_or_subregion = region,
                                 scale = "Region", 
                                 data_source = "RAM",
                                 variable = variables[v],
                                 ccf.calc = FALSE)
      }
      names(var.list) = variables
      data.points <- subset(alldat, region == "Benguela", datasource = data_source)
      # RAM only has both sardine and anchovy for any of the variables (e.g.,ss, rec) for Humboldt. 
      
      x = var.list$rec$Year # Time series sampling points (day of year in this case)
      y = var.list$rec[,c("Sardine.est","Anchovy.est")] # Time series only (each column = one species)
    
      w = mvcwt(x, y, min.scale = 0.25, max.scale = 4)
      mr = wmr(w)
      image(mr, reset.par = FALSE)
      contour(mr, bound = NA, add = TRUE)
      
      
      plot(1:nrow(y),y[,1],type='l',ylim=c(0,max(c(y[,1],y[,2]))))
      lines(1:nrow(y),y[,2],col='red')
      