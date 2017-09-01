library(lubridate); library(dplyr)

rec.data <- readRDS("rec_events.rds")

# Noise Summary
noise.data <- rec.data %>%
  filter(Description == 'Average noise',
         Date.Time > ymd('20161113')) %>%
  mutate(array = ifelse(grepl('I', Site), 'Inner',
                        ifelse(grepl('O', Site), 'Outer',
                               'MD WEA')),
         Data = as.numeric(Data),
         day = lubridate::floor_date(Date.Time, 'day'))

noise.data %>%
  group_by(day, Site, array) %>%
  summarize(max = max(Data)) %>%
  group_by(day, array) %>%
  summarize(avg = mean(max)) %>%
  group_by(array) %>%
  summarize(t = sum(avg > 300))


noise.data %>%
  group_by(day, array) %>%
  summarize(avg = mean(Data)) %>%
  group_by(array) %>%
  summarize(sum(avg > 300))



rec.data <- readRDS("rec_events.rds")
rec.data <- rec.data %>%
  filter(Description == 'Tilt angle',
         Date.Time > ymd('20161113')) %>%
  mutate(array = ifelse(grepl('I', Site), 'Inner',
                        ifelse(grepl('O', Site), 'Outer',
                               'MD WEA')),
         Data = as.numeric(Data),
         day = lubridate::floor_date(Date.Time, 'day')) %>%
  group_by(day, array) %>%
  summarize(avg = mean(Data))
