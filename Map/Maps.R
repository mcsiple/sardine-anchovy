# Try to get global map to put LME shapefile on
library(rgdal)
library(leaflet)
library(data.table)
library(mapview) #For saving leaflet maps as pdf or html 
mapwd <- "Map/LME66"

lmes <- readOGR(file.path(mapwd,"LMEs66.shp"))
# Look at the data
head(lmes@data)
class(lmes) # "SpatialPolygonsDataFrame"
names(lmes)
# Data frame with the degree of asynchrony - eventually this will be a csv
asyn.test <- data.frame(lme=c(3,13,21,22),
                        region.name=c("California","Humboldt","NE Atlantic","NE Atlantic"),
                        asynchrony = c(0.5,0.6,0.7,0.2))
gde_15<- lmes
data <- asyn.test
map_data_fortified <- fortify(gde_15, region = "LME_NUMBER") %>% 
  mutate(id = as.numeric(id)) #fortify is depracated
map_data_fortified <- broom::tidy(gde_15, region = "LME_NUMBER")
map_data <- map_data_fortified %>% left_join(data, by = c("id" = "bfs_id"))


# Using leaflet
bins <- c(0, 200, 500, 1000, Inf)
pal <- colorBin("YlOrRd", domain = lmes$Shape_Area, bins = bins)
(lmap <- leaflet(lmes) %>%
  addPolygons(
    fillColor = ~pal(Shape_Area),
    color = "darkgrey", 
              weight = 1, smoothFactor = 0.5, 
              opacity = 1.0, fillOpacity = 0.5))

mapshot(lmap, file = paste0(getwd(), "/map.pdf"),
        remove_controls = c("homeButton", "layersControl"))


# How to color polygons by some variable? in lmes@data or otherwise


# I got the following code from: https://github.com/rstudio/leaflet/issues/169
lmes.dt <- data.table(lmes@data)
asyn.test <- data.table(lme=c(3,13,21,22),
                        region.name=c("California","Humboldt","NE Atlantic","NE Atlantic"),
                        asynchrony = c(0.5,0.6,0.7,0.2))
# Create an explicit attribute to keep polygons IDs (useful to "re-attach" the table to the polygons later)
lmes.dt[, rn := row.names(lmes)]
setkey(lmes.dt, LME_NUMBER)
setkey(asyn.test, lme)
lmes.dt2 <- asyn.test[lmes.dt]

setkey(lmes.dt2, rn)
lmes@data <- data.frame(lmes.dt2[row.names(lmes)])
# Make sure polygons IDs and data.frame row.names match
bins <- c(0, 0.5,1)
pal <- colorBin("YlOrRd", domain = lmes$asynchrony, bins = bins)
(lmap <- leaflet(lmes) %>%
    addPolygons(
      fillColor = ~pal(asynchrony),
      color = "darkgrey", 
      weight = 1, smoothFactor = 0.5, 
      opacity = 1.0, fillOpacity = 0.5))


#states <- spChFIDs(states, states$rn)

# Trying to use something that isn’t leaflet ------------------------------
# Theme and other code can eb found here, but it's complicated and I don't think I need this much detail
#https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
# lmes2 <- sf::st_as_sf(lmes)
# plot(lmes2)
asyn.dat <- data.frame(LME_NUMBER=c(3,13,21,22),
           region.name=c("California","Humboldt","NE Atlantic","NE Atlantic"),
           asynchrony = c(0.5,0.6,0.7,0.2))
data <- asyn.dat
#map_data_fortified <- broom::tidy(lmes)
mapa <- readOGR(file.path(mapwd,"LMEs66.shp"))
mapa@data$id <- rownames(mapa@data)
mapa@data   <- dplyr::left_join(mapa@data, asyn.dat, by="LME_NUMBER") #
mapa.df     <- broom::tidy(mapa)
mapa.df     <- dplyr::left_join(mapa.df,mapa@data, by="id")

library(ggplot2)
theplot <- ggplot(mapa.df, aes(x=long, y=lat, group=group))+
  geom_polygon(aes(fill=asynchrony))+
  coord_fixed()

theplot
