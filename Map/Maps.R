# Try to get global map to put LME shapefile on
library(rgdal)
library(leaflet)

mapwd <- "Map/LME66"

lmes <- readOGR("LMEs66.shp")

bins <- c(0, 200, 500, 1000, Inf)
pal <- colorBin("YlOrRd", domain = lmes$Shape_Area, bins = bins)

leaflet(lmes) %>%
  addPolygons(
    fillColor = ~pal(Shape_Area),
    color = "darkgrey", 
              weight = 1, smoothFactor = 0.5, 
              opacity = 1.0, fillOpacity = 0.5)

# Look at the data
head(lmes@data)

class(lmes) # "SpatialPolygonsDataFrame"
names(lmes)

# How to color polygons by some variable? in lmes@data or otherwise