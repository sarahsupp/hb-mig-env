# Code for eBird migration project, to prepare the full eBird dataset into an "effort" dataset
# developed by S. Supp and L. Graham for new analysis combining eBird data with remotely sensed environment

# DISCLAIMER: Code is under development
# 17 July 2019

# Inputs: Full eBird dataset for all species 2006-2018

# Outputs: Effort dataset with columns for DAY, YEAR, CELL, COUNT, assuming a dggridR resolution 
#          of 6 (12452 km2) for the Fuller4H projection (not the default in ddgridR)


# Import libraries
library(tidyverse)
library(ggmap)
library(geosphere)
library(maps)
library(dggridR)
library(rgdal)

# read in full eBird data to be converted into effort data (downloaded on 22 August 2019)
effort = read_csv("") #FIXME, need full eBird to begin with

# clean up effort data by grouping aggregating the records from group record, and using the mean lat-longs in case of recorder error
effort <- effort %>%
  select(common_name, scientific_name, latitude, longitude, observation_date,
         time_observations_started, observer_id, protocol_type, project_code,
         duration_minutes, effort_distance_km, effort_area_ha, number_observers, 
         sampling_event_identifier, group_identifier) %>%
  mutate(sampling_event_identifier = 
           ifelse(is.na(group_identifier), sampling_event_identifier, group_identifier)) %>% 
  group_by(common_name, observation_date, sampling_event_identifier) %>%
  mutate(latitude2 = mean(latitude), longitude2 = mean(longitude), 
         latdiff = latitude-latitude2, londiff = longitude-longitude2) %>%
  distinct(sampling_event_identifier, .keep_all = TRUE)  %>%
  mutate(year = year(observation_date),
         month = month(observation_date),
         day = yday(observation_date))

# FAL suggests using Fuller grid with resolution=6 (12452 km2)
hex_6 <- dgconstruct(projection = "FULLER", aperture = 4, topology = "HEXAGON", res=6)

#Get the corresponding grid cells and center coordinates for all observations in the dataset (lat-long pair)
effort$cell <- dgGEO_to_SEQNUM(hex_6, effort$LONGITUDE, effort$LATITUDE)$seqnum
effort$cell_lat <- dgSEQNUM_to_GEO(hex_6, effort$cell)$lat_deg #FIXME: pretty sure this is accurate, check
effort$cell_lon <-dgSEQNUM_to_GEO(hex_6, effort$cell)$lon_deg #FIXME: pretty sure this is accurate, check

#Converting SEQNUM to GEO gives the center coordinates of the cells
cellcenters <- dgSEQNUM_to_GEO(hex_6, effort$cell)

#Get and plot the number of observations in each cell and in each year
ebirdcounts <- effort %>%
  group_by(cell, year) %>% 
  summarise(count=n())

ggplot(ebirdcounts, aes(x=count)) + 
  geom_histogram(binwidth=10) + 
  facet_wrap(~year)

#Get the grid cell boundaries for cells which had bird observations #FIXME: Will it be a problem that cell is repeated here, across multiple years? If so, how fix?
grid <- dgcellstogrid(hex_6, ebirdcounts$cell, frame=TRUE, wrapcells=TRUE)

#Update the grid cells' properties to include the number of observations in each cell #FIXME: and in each year
grid <- merge(grid, ebirdcounts, by.x="cell", by.y="cell")

#Get polygons for the spatial range and make a map of bird observations
countries <- map_data("usa") #add Canada and Mexico

ggplot() + 
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=count), alpha=0.4)    +
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  #  geom_point  (aes(x=cellcenters$lon_deg, y=cellcenters$lat_deg), size=0.5) +
  scale_fill_gradient(low="blue", high="red") + 
  facet_wrap(~year)




#dgearthgrid(hex_6, savegrid="cell_out-12.shp") # FIXME: Do I actually need this step or can I just pass the grid object directly onto the dataset, as with ISEA3H?

# TODO: FAL other code, to modify for our needs here:

#hex <- readOGR("cell_out-12.shp")
# hex_6@data$POLYFID <- 1:40962
# 
# # eBirdData$SUB_ID <- ifelse(!is.na(eBirdData$GROUP_ID), eBirdData$GROUP_ID, eBirdData$SUB_ID)
# # eBirdData <- aggregate(cbind(LATITUDE, LONGITUDE) ~ PRIMARY_COM_NAME + YEAR + DAY + SUB_ID, data=eBirdData, FUN=mean, na.rm=TRUE)
# 
# locs <- SpatialPoints(cbind(eBirdData$LONGITUDE, eBirdData$LATITUDE))
# proj4string(locs) <- proj4string(hex_6)
# ovr <- over(locs, hex_6)
# 
# eBirdData$POLYFID <- ovr$POLYFID

eft <- aggregate(SUB_ID ~ POLYFID + YEAR + DAY, data=eBirdData, FUN=length)