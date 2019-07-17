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
  auk_filter(file = "data/calhum.txt") %>%
  read_ebd()

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Calliope Hummingbird")

# BROAD TAILED HUMMINGBIRD
dat <- auk_ebd("data/ebd_brthum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/brthum.txt") %>%
  read_ebd()

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Broad-tailed Hummingbird")


# RUFOUS HUMMINGBIRD
dat <- auk_ebd("data/ebd_rufhum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/rufhum.txt") %>%
  read_ebd()

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Rufous Hummingbird")

# BLACK CHINNED HUMMINGBIRD
dat <- auk_ebd("data/ebd_bkchum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/bkchum.txt") %>%
  read_ebd()

day_summary <- dat %>% group_by(observation_date) %>%
  summarise(count = n())

ggplot(day_summary, aes(x = observation_date, y = count)) +
  geom_line() + ggtitle("Black-chinned Hummingbird")


# RUBY THROATED HUMMINGBIRD
dat <- auk_ebd("data/ebd_rthhum_200601_201812_relFeb-2019.txt") %>%
  auk_protocol(c("Stationary", "Traveling", "Incidental")) %>%
  auk_date(c("2008-01-01", "2018-12-31")) %>%
  auk_filter(file = "data/rthhum.txt") %>%
  read_ebd()

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