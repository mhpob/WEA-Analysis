library(lubridate); library(dplyr)

rec.data <- readRDS("data and imports/rec_events.rds")

# Temperature summary
wt.data <- rec.data %>%
  filter(Description == 'Average temperature',
         Date.Time > ymd('20161113')) %>%
  mutate(array = ifelse(grepl('I', Site), 'Inner',
                        ifelse(grepl('O', Site), 'Outer',
                               'MD WEA')),
         Data = as.numeric(Data),
         day = lubridate::floor_date(Date.Time, 'day'))

wt.data %>%
  filter(Date.Time > ymd('20170330')) %>%
  group_by(day, array) %>%
  summarize(daily = mean(Data))%>%
  ungroup() %>%

  # group_by(array) %>%
  summarize(min = min(daily),
            avg = mean(daily),
            max = max(daily))




# Noise Summary
noise.data <- rec.data %>%
  filter(Description == 'Average noise',
         Date.Time > ymd('20161113')) %>%
  mutate(array = ifelse(grepl('I', Site), 'Inner',
                        ifelse(grepl('O', Site), 'Outer',
                               'MD WEA')),
         Data = as.numeric(Data),
         day = lubridate::floor_date(Date.Time, 'day'))

# Number of days with loud noise pulses
noise.data %>%
  group_by(day, Site, array) %>%
  summarize(max = max(Data)) %>%
  group_by(day, array) %>%
  summarize(avg = mean(max)) %>%
  group_by(array) %>%
  summarize(sum(avg > 300))

# Number of days with loud average noise levels
noise.data %>%
  group_by(day, array) %>%
  summarize(avg = mean(Data)) %>%
  group_by(array) %>%
  summarize(sum(avg > 300))


# Tilt summary
tilt.data <- rec.data %>%
  filter(Description == 'Tilt angle',
         Date.Time > ymd('20161113')) %>%
  mutate(array = ifelse(grepl('I', Site), 'Inner',
                        ifelse(grepl('O', Site), 'Outer',
                               'MD WEA')),
         Data = as.numeric(Data),
         day = lubridate::floor_date(Date.Time, 'day'))

# Mean daily tilt angle
tilt.data %>%
  group_by(day, array) %>%
  summarize(daily = mean(Data))%>%
  ungroup() %>%
  summarize(min = min(daily),
            avg = mean(daily),
            max = max(daily))

tilt.data %>%
  group_by(day, array) %>%
  summarize(daily = mean(Data))%>%
  group_by(array) %>%
  summarize(min = min(daily),
            avg = mean(daily),
            var = var(daily),
            max = max(daily))
