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
library(lubridate)

# CALLIOPE HUMMINGBIRD
dat <- auk_ebd("data/ebd_calhum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/calhum.txt", overwrite=TRUE) %>%
  # unique=FALSE ensures that all records are brought in, not just one from each group identifier
  read_ebd(unique=FALSE) %>%
  # select the columns that will be used
  select(common_name, scientific_name, latitude, longitude, observation_date,
         time_observations_started, observer_id, protocol_type, project_code,
         duration_minutes, effort_distance_km, effort_area_ha, number_observers, 
         sampling_event_identifier, group_identifier) %>%
  # Reassigns Group records to Sampling Event records, so we may next account for discrepancies in recording location among group members.
  mutate(sampling_event_identifier = 
           ifelse(is.na(group_identifier), sampling_event_identifier, group_identifier)) %>% 
  group_by(common_name, observation_date, sampling_event_identifier) %>%
  # check for discrepancies in lat and lon recording, and calculate the mean value
  mutate(latitude2 = mean(latitude), longitude2 = mean(longitude), 
         latdiff = latitude-latitude2, londiff = longitude-longitude2) %>%
  distinct(sampling_event_identifier, .keep_all = TRUE) %>%
  mutate(year = year(observation_date),
         month = month(observation_date),
           day = yday(observation_date)) #TODO: Are times in eBird relative to the location, or converted already into UTC?
write_csv(dat, path = "data/calhum.csv")

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Calliope Hummingbird")

# BROAD TAILED HUMMINGBIRD
dat <- auk_ebd("data/ebd_brthum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/brthum.txt", overwrite=TRUE) %>%
  read_ebd(unique=FALSE) %>%
  select(common_name, scientific_name, latitude, longitude, observation_date,
         time_observations_started, observer_id, protocol_type, project_code,
         duration_minutes, effort_distance_km, effort_area_ha, number_observers, 
         sampling_event_identifier, group_identifier) %>%
  mutate(sampling_event_identifier = 
           ifelse(is.na(group_identifier), sampling_event_identifier, group_identifier)) %>% 
  group_by(common_name, observation_date, sampling_event_identifier) %>%
  mutate(latitude2 = mean(latitude), longitude2 = mean(longitude), 
         latdiff = latitude-latitude2, londiff = longitude-longitude2) %>%
  distinct(sampling_event_identifier, .keep_all = TRUE) %>%
  mutate(year = year(observation_date),
         month = month(observation_date),
         day = yday(observation_date))
write_csv(dat, path = "data/brthum.csv")

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Broad-tailed Hummingbird")


# RUFOUS HUMMINGBIRD
dat <- auk_ebd("data/ebd_rufhum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/rufhum.txt", overwrite=TRUE) %>%
  read_ebd(unique=FALSE) %>%
  select(common_name, scientific_name, latitude, longitude, observation_date,
         time_observations_started, observer_id, protocol_type, project_code,
         duration_minutes, effort_distance_km, effort_area_ha, number_observers, 
         sampling_event_identifier, group_identifier) %>%
mutate(sampling_event_identifier = 
           ifelse(is.na(group_identifier), sampling_event_identifier, group_identifier)) %>% #TODO: Check that this accomplishes what I want? For all records which are part of a groupID, the groupID is assigned to the sampling event ID #FIXME: I tested this on a smaller version of the dataframe and it seems to work fine. Is there a data issue?
  group_by(common_name, observation_date, sampling_event_identifier) %>%
  mutate(latitude2 = mean(latitude), longitude2 = mean(longitude), 
         latdiff = latitude-latitude2, londiff = longitude-longitude2) %>%
  distinct(sampling_event_identifier, .keep_all = TRUE) %>%
  mutate(year = year(observation_date),
         month = month(observation_date),
         day = yday(observation_date))
write_csv(dat, path = "data/rufhum.csv")

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Rufous Hummingbird")


# BLACK CHINNED HUMMINGBIRD
dat <- auk_ebd("data/ebd_bkchum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/bkchum.txt") %>%
  read_ebd(unique=FALSE) %>%
  select(common_name, scientific_name, latitude, longitude, observation_date,
         time_observations_started, observer_id, protocol_type, project_code,
         duration_minutes, effort_distance_km, effort_area_ha, number_observers, 
         sampling_event_identifier, group_identifier) %>%
  mutate(sampling_event_identifier = 
           ifelse(is.na(group_identifier), sampling_event_identifier, group_identifier)) %>% 
  group_by(common_name, observation_date, sampling_event_identifier) %>%
  mutate(latitude2 = mean(latitude), longitude2 = mean(longitude), 
         latdiff = latitude-latitude2, londiff = longitude-longitude2) %>%
  distinct(sampling_event_identifier, .keep_all = TRUE) %>%
  mutate(year = year(observation_date),
         month = month(observation_date),
         day = yday(observation_date))
write_csv(dat, path = "data/bkchum.csv")

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Black-chinned Hummingbird")


# RUBY THROATED HUMMINGBIRD
dat <- auk_ebd("data/ebd_rthhum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/rthhum.txt") %>%
  read_ebd(unique=FALSE) %>%
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
write_csv(dat, path = "data/rthhum.csv")

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Ruby-throated Hummingbird")


#------------------------

hbird = read_tsv("data/ebd_calhum_200601_201812_relFeb-2019/ebd_calhum_200601_201812_relFeb-2019.txt")
