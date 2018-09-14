library(lubridate); library(dplyr)

# Trying to get meteorological data from the NOAA NDBC.
# OCMD
# http://www.ndbc.noaa.gov/station_history.php?station=ocim2
# Offshore buoy
# http://www.ndbc.noaa.gov/station_page.php?station=44009

hist.read <- function(station, year){
  if(grepl(year, Sys.Date())){
    monthnum <- as.numeric(strftime(Sys.time(), '%m'))
    station <- tolower(station)
    urls <- c(paste0('https://www.ndbc.noaa.gov/view_text_file.php?filename=',
                     station,
                     seq(1, (monthnum - 3), 1),
                     year,
                     '.txt.gz&dir=data/stdmet/',
                     month.abb[1:(monthnum - 3)],
                     '/'),
              paste0('https://www.ndbc.noaa.gov/data/stdmet/',
                     month.abb[monthnum - 2],
                     '/', station, '.txt'))
    data <- lapply(urls, function(x){
      tryCatch(
        data.table::fread(x, skip = 2, na.strings = c('NA', 'MM'),
                          col.names = c('year', 'month', 'day', 'hr', 'min', 'wdir',
                                        'wspd', 'gst', 'wvht', 'dpd', 'apd', 'mwd',
                                        'pres', 'atmp', 'wtmp', 'dewp', 'vis', 'tide')),
        error = function(e) NULL)})
    data2 <- data.table::fread(paste0('https://www.ndbc.noaa.gov/data/realtime2/',
                                      toupper(station), '.txt'),
                    skip = 2, drop = 18, na.strings = c('NA', 'MM'),
                    col.names = c('year', 'month', 'day', 'hr', 'min', 'wdir',
                                  'wspd', 'gst', 'wvht', 'dpd', 'apd', 'mwd',
                                  'pres', 'atmp', 'wtmp', 'dewp', 'vis', 'tide'))
    data <- do.call(rbind, data)
    data <- rbind(data, data2)

  } else{
    url <- paste0('https://www.ndbc.noaa.gov/view_text_file.php?filename=',
                  station, 'h', year,
                  '.txt.gz&dir=data/historical/stdmet/')
    data <- data.table::fread(url, skip = 2, na.strings = c('NA', 'MM'),
                     col.names = c('year', 'month', 'day', 'hr', 'min', 'wdir',
                                   'wspd', 'gst', 'wvht', 'dpd', 'apd', 'mwd',
                                   'pres', 'atmp', 'wtmp', 'dewp', 'vis', 'tide'))
  }

  data$date.time <- paste(data$year, data$month, data$day, data$hr, data$min, sep = '-')
  data$date.time <- lubridate::ymd_hm(data$date.time, tz = 'America/New_York')

  data$station <- station
  data$pres[data$pres == 9999] <- NA

  temp.data99 <- data[, c('wspd', 'gst', 'wvht', 'dpd', 'apd')]
  temp.data99[temp.data99 == 99] <- NA

  temp.data999 <- data[, c('wdir', 'mwd', 'atmp', 'wtmp', 'dewp')]
  temp.data999[temp.data999 == 999] <- NA

  data <- cbind(data[, c('date.time', 'station', 'pres')],
                temp.data99, temp.data999)
  data <- data[, c('date.time', 'station', 'wdir', 'wspd', 'gst', 'wvht', 'dpd',
                   'apd', 'mwd', 'pres', 'atmp', 'wtmp', 'dewp')]

  data
}

ocmd.2016 <- hist.read('ocim2', 2016)
ocmd.2017 <- hist.read('ocim2', 2017)
ocmd.2018 <- hist.read('ocim2', 2018)
buoy.2016 <- hist.read('44009', 2016)
buoy.2017 <- hist.read('44009', 2017)
buoy.2018 <- hist.read('44009', 2018)

met.data <- rbind(buoy.2016, buoy.2017, buoy.2018, ocmd.2016, ocmd.2017, ocmd.2018)

saveRDS(met.data, file = 'data and imports/NDBC_data.rds')
rm(buoy.2016, buoy.2017, buoy.2018, ocmd.2016, ocmd.2017, ocmd.2018)


# hrly.wind <- met.data %>%
#   mutate(hr.agg = floor_date(date.time, unit = 'hour')) %>%
#   group_by(hr.agg) %>%
#   summarize(mean.wind = mean(wspd))
#
# ggplot(data = hrly.wind, aes(x = hr.agg, y = mean.wind)) +
#   geom_line()
#
# hrly6.wind <- met.data %>%
#   mutate(hr.agg = floor_date(date.time, unit = '6hour')) %>%
#   group_by(hr.agg) %>%
#   summarize(mean.wspd = mean(wspd))
#
# ggplot(data = hrly6.wind, aes(x = hr.agg, y = mean.wspd)) +
#   geom_line() +
#   scale_x_datetime(date_breaks = '1 month', date_labels = '%b',
#                    minor_breaks = NULL)

# hrly12.wind <- met.data %>%
#   mutate(hr.agg = floor_date(date.time, unit = '12hour')) %>%
#   group_by(hr.agg) %>%
#   summarize(mean.wspd = mean(wspd))
#
# ggplot(data = hrly12.wind, aes(x = hr.agg, y = mean.wspd)) +
#   geom_line() +
#   scale_x_datetime(date_breaks = '1 month', date_labels = '%b',
#                    minor_breaks = NULL)

