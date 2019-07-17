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

# read in full eBird data to be converted into effort data
effort = #FIXME, need full eBird to begin with

# FAL suggests using Fuller grid with resolution=6 (12452 km2)
hex_6 <- dgconstruct(projection = "FULLER", aperture = 4, topology = "HEXAGON", res=6)

dgearthgrid(hex_6, savegrid="cell_out-12.shp") # FIXME: Do I actually need this step or can I just pass the grid object directly onto the dataset, as with ISEA3H?


# TODO: FAL other code, to modify for our needs here:
library(rgdal)

hex <- readOGR("cell_out-12.shp")
hex@data$POLYFID <- 1:40962

eBirdData$SUB_ID <- ifelse(!is.na(eBirdData$GROUP_ID), eBirdData$GROUP_ID, eBirdData$SUB_ID)
eBirdData <- aggregate(cbind(LATITUDE, LONGITUDE) ~ PRIMARY_COM_NAME + YEAR + DAY + SUB_ID, data=eBirdData, FUN=mean, na.rm=TRUE)

locs <- SpatialPoints(cbind(eBirdData$LONGITUDE, eBirdData$LATITUDE))
proj4string(locs) <- proj4string(hex)
ovr <- over(locs, hex)

eBirdData$POLYFID <- ovr$POLYFID

eft <- aggregate(SUB_ID ~ POLYFID + YEAR + DAY, data=eBirdData, FUN=length)