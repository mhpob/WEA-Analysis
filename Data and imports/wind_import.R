library(ggplot2); library(lubridate); library(dplyr)

# Trying to get wind data from the NOAA NDBC. This will inevitably break. Sorry.
# http://www.ndbc.noaa.gov/station_history.php?station=ocim2

ocmd.2016 <- read.table('http://www.ndbc.noaa.gov/view_text_file.php?filename=ocim2h2016.txt.gz&dir=data/historical/stdmet/')


URLs.2017 <- paste0('http://www.ndbc.noaa.gov/view_text_file.php?filename=ocim2',
                    seq(1, 7, 1),
                    '2017.txt.gz&dir=data/stdmet/',
                    c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul'),
                    '/')

ocmd.2017.1_7 <- lapply(URLs.2017, read.table)
ocmd.2017.1_7 <- do.call(rbind, ocmd.2017.1_7)


ocmd.2017.8 <- read.table('http://www.ndbc.noaa.gov/data/stdmet/Aug/ocim2.txt')


ocmd.wind <- rbind(ocmd.2016, ocmd.2017.1_7, ocmd.2017.8)

names(ocmd.wind) <- c('year', 'month', 'day', 'hr', 'min', 'wdir',
                      'wspd', 'gst', 'wvht', 'dpd', 'apd', 'mwd',
                      'pres', 'atmp', 'wtmp', 'dewp', 'vis', 'tide')


rm(ocmd.2016, URLs.2017, ocmd.2017.1_7, ocmd.2017.8)

ocmd.wind <- ocmd.wind %>%
  mutate(date.time = paste(year, month, day,
                           hr, min),
         date.time = ymd_hm(date.time)) %>%
  select(date.time, wdir, wspd, gst, pres, atmp, wtmp) %>%
  filter(wspd != 99,
         date.time >= '2016-11-01')

# saveRDS(ocmd.wind, file = 'data and imports/ocmd_wind.rds')

# hrly.wind <- ocmd.wind %>%
#   mutate(hr.agg = floor_date(date.time, unit = 'hour')) %>%
#   group_by(hr.agg) %>%
#   summarize(mean.wind = mean(wspd))
#
# ggplot(data = hrly.wind, aes(x = hr.agg, y = mean.wind)) +
#   geom_line()
#
# hrly6.wind <- ocmd.wind %>%
#   mutate(hr.agg = floor_date(date.time, unit = '6hour')) %>%
#   group_by(hr.agg) %>%
#   summarize(mean.wspd = mean(wspd))
#
# ggplot(data = hrly6.wind, aes(x = hr.agg, y = mean.wspd)) +
#   geom_line() +
#   scale_x_datetime(date_breaks = '1 month', date_labels = '%b',
#                    minor_breaks = NULL)

hrly12.wind <- ocmd.wind %>%
  mutate(hr.agg = floor_date(date.time, unit = '12hour')) %>%
  group_by(hr.agg) %>%
  summarize(mean.wspd = mean(wspd))

ggplot(data = hrly12.wind, aes(x = hr.agg, y = mean.wspd)) +
  geom_line() +
  scale_x_datetime(date_breaks = '1 month', date_labels = '%b',
                   minor_breaks = NULL)

