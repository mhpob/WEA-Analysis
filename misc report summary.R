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
         day = lubridate::floor_date(Date.Time, 'day'),
         year = ifelse(month(Date.Time) %in% 10:12 |
                         (month(Date.Time) == 9 & day(Date.Time) %in% 16:30),
                       year(Date.Time) - 1999, year(Date.Time) - 2000),
         season = ifelse(month(Date.Time) %in% 4:8 |
                           (month(Date.Time)==3) |
                           (month(Date.Time)==8),
                         'SprSum', 'AutWin'),
         season = paste0(season, year))


wt.data %>%
  filter(season== "SprSum18") %>%
  group_by(day, array) %>%
  summarize(daily = mean(Data))%>%
  ungroup() %>%

  # group_by(array) %>%
  summarize(min = min(daily),
            avg = mean(daily),
            max = max(daily))

wt.data %>%
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
         day = lubridate::floor_date(Date.Time, 'day'),
         year = ifelse(month(Date.Time) %in% 10:12 |
                         (month(Date.Time) == 9 & day(Date.Time) %in% 16:30),
                       year(Date.Time) - 1999, year(Date.Time) - 2000),
         season = ifelse(month(Date.Time) %in% 4:8 |
                           (month(Date.Time) == 3) |
                           (month(Date.Time) == 8),
                         'SprSum', 'AutWin'),
         season = paste0(season, year))

noise.data %>%
  group_by(day, array) %>%
  summarize(daily = mean(Data))%>%
  group_by(array) %>%
  summarize(min = min(daily),
            avg = mean(daily),
            max = max(daily))

noise.data %>%
  group_by(day, array) %>%
  summarize(daily = mean(Data))%>%
  ungroup() %>%
  summarize(min = min(daily),
            avg = mean(daily),
            max = max(daily))

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
         day = lubridate::floor_date(Date.Time, 'day'),
         year = ifelse(month(Date.Time) %in% 10:12 |
                         (month(Date.Time) == 9 & day(Date.Time) %in% 16:30),
                       year(Date.Time) - 1999, year(Date.Time) - 2000),
         season = ifelse(month(Date.Time) %in% 4:8 |
                           (month(Date.Time) == 3 & day(Date.Time) %in% 16:31) |
                           (month(Date.Time) == 9 & day(Date.Time) %in% 1:15),
                         'SprSum', 'AutWin'),
         season = paste0(season, year))

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
