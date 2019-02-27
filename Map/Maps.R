# Try to get global map to put LME shapefile on
library(rgdal)
library(leaflet)

mapwd <- "Map/LME66"

lmes <- readOGR("LMEs66.shp")
leaflet(lmes) %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5)
