library(ggplot2); library(lubridate); library(dplyr)

# Noise & tilt munging ----
rec <- readRDS("data and imports/rec_events.rds")

noise_tilt <- rec %>%
  filter(Description %in% c('Average noise', 'Tilt angle',
                            'Average temperature'),
         Date.Time > '2016-11-15') %>%
  mutate(Data = as.numeric(Data),
         hr12 = floor_date(Date.Time, unit = '12hour')) %>%
  group_by(hr12, Description, Site, Lat, Long) %>%
  summarize(avg = mean(Data)) %>%
  as.data.frame()


#Aggregate by array
arr <- function(data, part){grepl(part, data[, 'Site'])}
noise_tilt$Array <- ifelse(arr(noise_tilt, 'A'), 'MD WEA',
                           ifelse(arr(noise_tilt, 'O'), 'Outer',
                                  'Inner'))

agg_noise_tilt <- aggregate(avg ~ hr12 + Description + Array,
                            data = noise_tilt, FUN = mean) %>%
  filter(Array == 'MD WEA') %>%
  select(-Array)

# Wind munging ----
ocmd.wind <- readRDS("data and imports/ocmd_wind.rds")

hr12.wind <- ocmd.wind %>%
  filter(date.time >= '2016-11-13',
         date.time <= '2017-08-26') %>%
  mutate(hr12 = floor_date(date.time, unit = '12hour')) %>%
  group_by(hr12) %>%
  summarize(avg = mean(wspd),
            Description = 'Wind speed')

# Combine data ----
phys_data <- rbind(agg_noise_tilt, hr12.wind)
phys_data$Description <- factor(phys_data$Description,
                                levels = c('Wind speed', 'Average noise',
                                           'Tilt angle', 'Average temperature'),
                                ordered = T)

ggplot() + geom_line(data = phys_data, aes(x = hr12, y = avg)) +
  facet_wrap(~ Description, ncol = 1, scales = 'free_y') +
  scale_x_datetime(date_breaks = '1 month', date_labels = '%b',
                   minor_breaks = NULL) +
  labs(x = 'Date', y = '12-hour mean value') +
  theme_bw()
