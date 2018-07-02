
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
       # function(x) round(x/diag(x), 2))

calc_det_freq <- function(data, time_unit){
  spl_time <- split(data, data[time_unit])
  spl_t.array <- lapply(spl_time, function(x) split(x, x$array))

  # Count transmission/detection combinations, divide by total transmissions
  spl_time <- lapply(spl_t.array, function(x){

    temp_x <- lapply(x, function(y){
      temp_y <- xtabs(~ internal + station, data = y)
      temp_y <- round(temp_y / diag(temp_y), 2)
      as.data.frame(temp_y, stringsAsFactors = F)
    })

    for(i in 1:2){
      temp_x[[i]]$array <- names(temp_x)[i]
    }

    temp_x <- do.call(rbind, temp_x)
    row.names(temp_x) <- NULL

    temp_x$d1 <- as.numeric(sapply(
      strsplit(temp_x$internal, '[_]'), '[', 2))
    temp_x$d2 <- as.numeric(sapply(
      strsplit(temp_x$station, '[_]'), '[', 2))
    temp_x[is.na(temp_x)] <- 0
    temp_x$distance <- abs(temp_x$d1 - temp_x$d2)

    temp_x[, names(temp_x) %in% c('array', 'distance', 'Freq')]
  })

  for(i in seq_along(spl_time)){
    spl_time[[i]]$date <- names(spl_time)[i]
  }

  df.data <- do.call(rbind, spl_time)
  row.names(df.data) <- NULL

  if(time_unit == 'day'){
    df.data$date <- as.Date(df.data$date)
  } else {
    df.data$date <- lubridate::ymd_hms(df.data$date, tz = 'America/New_York')
  }

  names(df.data) <- tolower(names(df.data))
  df.data[, c('date', 'array', 'distance', 'freq')]
}

array.spl <- lapply(array.spl, calc_det_freq, time_unit = 'day')
for(i in seq_along(array.spl)){
  array.spl[[i]]$array <- names(array.spl)[i]
  array.spl[[i]]$Freq[array.spl[[i]]$Freq > 1] <- 1
}

# Modeling ----
lapply(array.spl, function(x){
  glm(Freq ~ distance, family = 'binomial', data = x)
})

d50_glm <- function(data){
  spl.data <- split(data, data$date)
  spl.data <- lapply(spl.data, function(x){
    glm(Freq ~ distance, data = x, family = 'binomial')
  })
  d50 <- lapply(spl.data, function(x) as.numeric(-coef(x)[1]/coef(x)[2]))
  d50 <- data.frame(d50 = do.call(rbind, d50))
  d50$date <- as.Date(row.names(d50))
  row.names(d50) <- NULL
  d50
}

j <- d50_glm(array.spl[[1]])
