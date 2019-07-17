# Code for eBird migration project, to prepare eBird species data for analysis
# developed by S. Supp and L. Graham for new analysis combining eBird data with remotely sensed environment

# DISCLAIMER: Code is under development
# 17 July 2019

# Inputs: Downloaded eBird dataset for each species 2006-2018

# Outputs: Species dataset with columns for SCIENTIFIC NAME, COMMON NAME, LATITUDE, LONGITUDE, 
#      OBSERVATION DATE (change to Year and Julian Day of Year columns), TIME OBSERVATIONS STARTED,
#      OBSERVER ID, PROTOCOL TYPE, PROJECT CODE, DURATION MINUTES (previously was hours), 
#      EFFORT DISTANCE KM, EFFORT AREA HA, NUMBER OBSERVERS


# Import libraries

library(auk)
library(tidyverse)

# CALLIOPE HUMMINGBIRD
dat <- auk_ebd("data/ebd_calhum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/calhum.txt", overwrite=TRUE) %>%
  read_ebd() %>%
  select(common_name, scientific_name, latitude, longitude, observation_date,
         time_observations_started, observer_id, protocol_type, project_code,
         duration_minutes, effort_distance_km, effort_area_ha, number_observers, 
         sampling_event_identifier, group_identifier) %>%
  mutate(sampling_event_identifier = 
           ifelse(!is.na(group_identifier), group_identifier, sampling_event_identifier)) %>% #TODO: Check, but seems to accomplish what I want? For all records which are part of a groupID, the groupID is assigned to the sampling event ID
  group_by(common_name, observation_date, sampling_event_identifier) %>%
  mutate(latitude2 = mean(latitude), longitude2 = mean(longitude), 
         latdiff = latitude-latitude2, londiff = longitude-longitude2) #TODO: no differences observed... is this because groups had total agreement (realistic)? or because I've done something wrong. 
# FIXME: Also, I should have some rows with the same groupID but that doesn't seem to be the case? Because I used mutate, some duplicates should be present, will need to remove duplicate group_identifiers, now that the means have been calculated

#FAL code for dealing with group_identifier
# eBirdData$SUB_ID <- ifelse(!is.na(eBirdData$GROUP_ID), eBirdData$GROUP_ID, eBirdData$SUB_ID)
# eBirdData <- aggregate(cbind(LATITUDE, LONGITUDE) ~ PRIMARY_COM_NAME + YEAR + DAY + SUB_ID, data=eBirdData, FUN=mean, na.rm=TRUE)


day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Calliope Hummingbird")

# BROAD TAILED HUMMINGBIRD
dat <- auk_ebd("data/ebd_brthum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/brthum.txt") %>%
  read_ebd() %>%
  select(common_name, scientific_name, latitude, longitude, observation_date,
         time_observations_started, observer_id, protocol_type, project_code,
         duration_minutes, effort_distance_km, effort_area_ha, number_observers, 
         sampling_event_identifier, group_identifier)

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Broad-tailed Hummingbird")


# RUFOUS HUMMINGBIRD
dat <- auk_ebd("data/ebd_rufhum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/rufhum.txt", overwrite=TRUE) %>%
  read_ebd() %>%
  select(common_name, scientific_name, latitude, longitude, observation_date,
         time_observations_started, observer_id, protocol_type, project_code,
         duration_minutes, effort_distance_km, effort_area_ha, number_observers, 
         sampling_event_identifier, group_identifier) %>%
mutate(sampling_event_identifier = 
           ifelse(is.na(group_identifier), sampling_event_identifier, group_identifier)) %>% #TODO: Check that this accomplishes what I want? For all records which are part of a groupID, the groupID is assigned to the sampling event ID #FIXME: I don't think this is actually working, because it says all items are unique, which shouldn't be true
  group_by(common_name, observation_date, sampling_event_identifier) %>%
  mutate(latitude2 = mean(latitude), longitude2 = mean(longitude), 
         latdiff = latitude-latitude2, londiff = longitude-longitude2) %>%
  distinct(sampling_event_identifier, .keep_all = TRUE)
#TODO: no differences observed... is this because groups had total agreement (realistic)? or because I've done something wrong. Also, I should have some rows with the same groupID but that doesn't seem to be the case?


day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Rufous Hummingbird")

# BLACK CHINNED HUMMINGBIRD
dat <- auk_ebd("data/ebd_bkchum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/bkchum.txt") %>%
  read_ebd() %>%
  select(common_name, scientific_name, latitude, longitude, observation_date,
         time_observations_started, observer_id, protocol_type, project_code,
         duration_minutes, effort_distance_km, effort_area_ha, number_observers, 
         sampling_event_identifier, group_identifier)

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Black-chinned Hummingbird")


# RUBY THROATED HUMMINGBIRD
dat <- auk_ebd("data/ebd_rthhum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/rthhum.txt") %>%
  read_ebd() %>%
  select(common_name, scientific_name, latitude, longitude, observation_date,
         time_observations_started, observer_id, protocol_type, project_code,
         duration_minutes, effort_distance_km, effort_area_ha, number_observers, 
         sampling_event_identifier, group_identifier)

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Ruby-throated Hummingbird")


#------------------------

hbird = read_tsv("data/ebd_calhum_200601_201812_relFeb-2019/ebd_calhum_200601_201812_relFeb-2019.txt")
#TODO: Keep only 1 record from each group identifier
#TODO: Keep only columns for SCIENTIFIC NAME, COMMON NAME, LATITUDE, LONGITUDE, 
#      OBSERVATION DATE (change to Year and Julian Day of Year columns), TIME OBSERVATIONS STARTED,
#      OBSERVER ID, PROTOCOL TYPE, PROJECT CODE, DURATION MINUTES (previously was hours), 
#      EFFORT DISTANCE KM, EFFORT AREA HA, NUMBER OBSERVERS
# TODO: Would this kind of clean up be better to do in a separate script, 
#       then import the cleaned files directly here? Seems better...?