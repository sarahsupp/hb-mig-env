# Code for eBird migration project, to prepare the full eBird dataset into an "effort" dataset
# developed by S. Supp and L. Graham for new analysis combining eBird data with remotely sensed environment

# DISCLAIMER: Code is under development
# 17 July 2019

# Inputs: Full eBird dataset for all species 2006-2018

# Outputs: Effort dataset with columns for DAY, YEAR, CELL, COUNT, assuming a dggridR resolution 
#          of 8 (7774.20 km2) or 7 (23322.62 km2). 

# Import libraries
library(tidyverse)
library(ggmap)
library(geosphere)
library(maps)
library(dggridR)

# read in full eBird data to be converted into effort data
effort = #FIXME, need effort file