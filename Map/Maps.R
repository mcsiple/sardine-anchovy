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

# Color polygons by asynchrony variable in lmes@data or otherwise
# I got the following code from: https://github.com/rstudio/leaflet/issues/169
lmes.dt <- data.table(lmes@data)

# ***** Replace this below table with the full asynchrony table, as soon as we have a good metric
asyn.test <- data.table(lme=c(3,13,21,22),
                        region.name=c("California",
                                      "Humboldt",
                                      "NE Atlantic",
                                      "NE Atlantic"),
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

# Using MAPS
library(maps)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Assign LMEs to study regions --------------------------------------------
mylmes <- data.frame(lme=c(3,13,29,21,22,24,59,60,20,49,50,51))
mylmes$region <- c("California","Humboldt","Benguela",
                   rep("NE Atlantic",times=6),
                   rep("Kuroshio-Oyashio",times=3))

lmes3 <- sf::st_as_sf(lmes) %>% #select(region.name)
  filter(lme %in% unique(mylmes$lme)) #hoping that this speeds things up, but who knows

pal<- brewer.pal(n = 10,name = "RdBu")

# Now add countries in the background -------------------------------------
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
ggplot(data=world) + 
  geom_sf(color="lightgrey",fill="lightgrey") +
  geom_sf(data = lmes3,aes(fill=asynchrony),size=0.1) + 
  scale_fill_gradientn(colours = pal) + 
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "white"))


# Add asynchrony information to table with lmes ---------------------------
load(file.path(datwd,"LogDiffsMax.RData")) #df is st, summary of log diffs.
# ml column is log ratio of anchovy/sardine. It's the median of the log ratios from all data soures

# Asynchrony table
head(st)
replaceability <- left_join(st,mylmes,by="region") 
ssbreplace.dt <- filter(replaceability,newvar=="Spawning stock biomass") %>% as.data.table()



# Trying to use something that isn’t leaflet ------------------------------
# Theme and other code can eb found here, but it's complicated and I don't think I need this much detail
#https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
# lmes2 <- sf::st_as_sf(lmes)
# plot(lmes2)
# asyn.dat <- data.frame(LME_NUMBER=c(3,13,21,22),
#                        region.name=c("California","Humboldt","NE Atlantic","NE Atlantic"),
#                        asynchrony = c(0.5,0.6,0.7,0.2))
# data <- asyn.dat
# #map_data_fortified <- broom::tidy(lmes)
# mapa <- readOGR(file.path(mapwd,"LMEs66.shp"))
# mapa@data$id <- rownames(mapa@data)
# mapa@data   <- dplyr::left_join(mapa@data, asyn.dat, by="LME_NUMBER") #
# mapa.df     <- broom::tidy(mapa)
# mapa.df     <- dplyr::left_join(mapa.df,mapa@data, by="id")
# 
# library(ggplot2)
# theplot <- ggplot(mapa.df, aes(x=long, y=lat, group=group))+
#   geom_polygon(aes(fill=asynchrony))+
#   coord_fixed()
# 
# theplot

