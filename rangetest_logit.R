
# Import range test ping data ----
rangetestimport <- function(path, pattern){
  filenames <- list.files(path = path, pattern = pattern, full.names = T)
  filenames <- grep('csv', filenames, value = T)
  detections <- lapply(filenames, data.table::fread,
                       colClasses = c('Sensor Value' = 'numeric',
                                      'Sensor Unit' = 'character'))
  detections <- as.data.frame(do.call(rbind, detections))
  names(detections) <- c('date.utc', 'receiver', 'transmitter', 't.name',
                         't.serial', 'sensor.value', 'sensor.unit', 'station',
                         'lat', 'long')
  detections$date.utc <- as.POSIXct(detections$date.utc,
                                    format = '%Y-%m-%d %H:%M:%S',
                                    tz = 'UTC')
  detections$date.east <- detections$date.utc
  attr(detections$date.east, 'tzone') <- 'America/New_York'
  detections
}

rangetest <- rbind(
  rangetestimport('p:/obrien/biotelemetry/detections/offshore md/Fish Migration/201804',
                  '30[378]'), #inner (IS2) range test
  rangetestimport('p:/obrien/biotelemetry/detections/offshore md/Fish Migration/201804',
                  '30[569]')) #wea (AN3) range test

# Data munging ----
library(dplyr)
data <- rangetest %>%
  filter(date.east > '2017-12-21',
         date.east < '2018-04-11',
         grepl('-60', transmitter)) %>%
  mutate(array = ifelse(grepl('A', station), 'MD WEA', 'Inner'),
         internal = case_when(grepl('60767', transmitter) ~ 'IS2',
                              grepl('60771', transmitter) ~ 'IS2_250',
                              grepl('60772', transmitter) ~ 'IS2_800',
                              grepl('60769', transmitter) ~ 'AN3_250',
                              grepl('60770', transmitter) ~ 'AN3_800',
                              grepl('60773', transmitter) ~ 'AN3'),
         day = lubridate::date(date.east),
         hour6 = lubridate::floor_date(date.east, unit = '6hour'))

array.spl <- split(data, data$array)
# lapply(lapply(array.spl, function(x) table(x$internal, x$station)),
#        function(x) round(x/diag(x), 2))

calc_det_freq <- function(data, time_unit){
  spl_time <- split(data, data[time_unit])

  # Count transmission/detection combinations, divide by total transmissions
  spl_time <- lapply(spl_time, function(x){
    temp_x <- xtabs(~ internal + station, data = x)
    temp_x <- round(temp_x / diag(temp_x), 2)
    as.data.frame(temp_x, stringsAsFactors = F)
  })

  for(i in seq_along(spl_time)){
    spl_time[[i]]$date <- names(spl_time)[i]
  }

  df.data <- do.call(rbind, spl_time)
  row.names(df.data) <- NULL

  df.data$d1 <- as.numeric(sapply(
      strsplit(df.data$internal, '[_]'), '[', 2))
  df.data$d2 <- as.numeric(sapply(
      strsplit(df.data$station, '[_]'), '[', 2))
  df.data[is.na(df.data)] <- 0
  df.data$Freq[df.data$Freq > 1] <- 1
  df.data$distance <- abs(df.data$d1 - df.data$d2)


  if(time_unit == 'day'){
    df.data$date <- as.Date(df.data$date)
  } else {
    df.data$date <- lubridate::ymd_hms(df.data$date, tz = 'America/New_York')
  }

  names(df.data) <- tolower(names(df.data))
  df.data[, c('date', 'distance', 'freq')]
}

array.spl <- lapply(array.spl, calc_det_freq, time_unit = 'day')
for(i in seq_along(array.spl)){
  array.spl[[i]]$array <- names(array.spl)[i]
}

# Modeling ----
lapply(array.spl, function(x){
  glm(freq ~ distance, family = 'binomial', data = x)
})

d50_glm <- function(data){
  spl.data <- split(data, data$date)
  d50 <- lapply(spl.data, function(x){
    model_fit <- glm(freq ~ distance, data = x, family = 'binomial')
    as.numeric(-coef(model_fit)[1]/coef(model_fit)[2])
  })
  d50 <- data.frame(d50 = do.call(rbind, d50))
  d50$date <- as.Date(row.names(d50))
  row.names(d50) <- NULL
  d50[, c('date', 'd50')]
}

d50 <- lapply(array.spl, d50_glm)
for(i in seq_along(d50)) d50[[i]]$array <- names(d50)[i]
d50 <- do.call(rbind, d50)

library(ggplot2)
ggplot() + geom_line(data = d50, aes(x = date, y = d50, color = array), lwd = 2) +
  labs(x = NULL, y = 'D50 (m)', color = 'Array') +
  theme_bw()

d50plot <- function(array, day){
  ggplot(data = filter(array.spl[[array]], date == day),
         aes(x = distance, y = freq)) +
    geom_point() +
    geom_smooth(method = 'glm', method.args = list(family = 'binomial'), se = F) +
    geom_point(data = filter(d50, array == array, date == day),
               aes(x = d50, y = 0.5), col = 'red', size = 3) +
    geom_segment(data = filter(d50, array == array, date == day),
                 aes(x = 0, y = 0.5, xend = d50, yend = 0.5)) +
    geom_segment(data = filter(d50, array == array, date == day),
                 aes(x = d50, y = 0, xend = d50, yend = 0.5)) +
    lims(x = c(0, 1000), y = c(0, 1)) +
    labs(x = 'Distance', y = 'Frequency of Detection') +
    coord_flip()
}

d50plot('Inner', '2018-03-20')
