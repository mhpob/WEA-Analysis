library(mcp)

library(dplyr)
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
  detections <- detections[grepl('-60', detections$transmitter),]
  detections$date.utc <- as.POSIXct(detections$date.utc,
                                    format = '%Y-%m-%d %H:%M:%S',
                                    tz = 'UTC')
  detections$date.east <- detections$date.utc
  attr(detections$date.east, 'tzone') <- 'America/New_York'
  detections
}

rangetest_201804 <- rbind(
  #inner (IS2) range test
  rangetestimport('p:/obrien/biotelemetry/detections/offshore md/Fish Migration/201804',
                  '30[378]'),
  #wea (AN3) range test
  rangetestimport('p:/obrien/biotelemetry/detections/offshore md/Fish Migration/201804',
                  '30[569]')) %>%
  mutate(internal = case_when(grepl('60767', transmitter) ~ 'IS2',
                              grepl('60771', transmitter) ~ 'IS2_250',
                              grepl('60772', transmitter) ~ 'IS2_800',
                              grepl('60773', transmitter) ~ 'AN3',
                              grepl('60769', transmitter) ~ 'AN3_250',
                              grepl('60770', transmitter) ~ 'AN3_800'))

# Need to do Aug18 & Dec18-downloaded data separately since receivers moved. Not
# assigning sites to transmitters here; going to use the NAs above as a signal to
# reassign after movement.
rangetest_201812 <- rbind(
  #inner (IS2) range test
  rangetestimport('p:/obrien/biotelemetry/detections/offshore md/Fish Migration/201808',
                  '30[367]'),
  rangetestimport('p:/obrien/biotelemetry/detections/offshore md/Fish Migration/201812',
                  '30[367]'),
  #wea (AN3) range test
  rangetestimport('p:/obrien/biotelemetry/detections/offshore md/Fish Migration/201808',
                  '30[59]|461'),
  rangetestimport('p:/obrien/biotelemetry/detections/offshore md/Fish Migration/201812',
                  '30[59]|461'))


# Data munging ----
data <- bind_rows(rangetest_201804, rangetest_201812) %>%
  distinct(date.utc, receiver, transmitter, .keep_all = T) %>%
  filter(date.east > '2017-12-21',
         date.east < '2018-12-05') %>%
  mutate(array = ifelse(grepl('A', station), 'MD WEA', 'Inner'),
         internal = case_when(is.na(internal) &
                                grepl('60767', transmitter) ~ 'IS2',
                              is.na(internal) &
                                grepl('60771', transmitter) ~ 'IS2_250',
                              is.na(internal) &
                                grepl('60770', transmitter) ~ 'IS2_800',
                              is.na(internal) &
                                grepl('60769', transmitter) ~ 'AN3',
                              is.na(internal) &
                                grepl('60773', transmitter) ~ 'AN3_250',
                              is.na(internal) &
                                grepl('60925', transmitter) ~ 'AN3_800',
                              # There was an error in redeployment from 201804-201808:
                              # the AN3_800 transmitter was not activated.
                              T ~ internal),
         day = lubridate::date(date.east),
         hour6 = lubridate::floor_date(date.east, unit = '6hour'))

array.spl <- split(data, data$array)
lapply(lapply(array.spl, function(x) table(x$internal, x$station)),
       function(x) round(x/diag(x), 2))

calc_det_freq <- function(data, time_unit){
  spl_time <- split(data, data[time_unit])

  # Count transmission/detection combinations, divide by total transmissions
  spl_time <- lapply(spl_time, function(x){
    temp <- xtabs(~ internal + station, data = x)
    temp_df <- as.data.frame(temp, stringsAsFactors = F)
    temp_df$trials <- rep(diag(temp), each = 3)
    temp_df$fail <- temp_df$trials - temp_df$Freq
    temp_df
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
  df.data$distance <- abs(df.data$d1 - df.data$d2)

  df.data$freq <- df.data$Freq / df.data$trials
  df.data$freq[df.data$freq > 1] <- 1

  df.data$fail[df.data$fail < 0] <- 0


  if(time_unit == 'day'){
    df.data$date <- as.Date(df.data$date)
  } else {
    df.data$date <- lubridate::ymd_hms(df.data$date, tz = 'America/New_York')
  }

  names(df.data)[names(df.data) == 'Freq'] <- 'success'
  df.data[, c('date', 'distance', 'trials', 'success', 'fail', 'freq')]
}

array.spl <- lapply(array.spl, calc_det_freq, time_unit = 'day')
for(i in seq_along(array.spl)){
  array.spl[[i]]$array <- names(array.spl)[i]
}

det.freq <- bind_rows(array.spl)
data_subs <- det.freq
data_subs$index <- as.numeric(data_subs$date) - 17520
data_subs <- data_subs %>%
  group_by(index distance) %>%
  summarize(success = sum(success),
            trials = sum(trials))

model = list(
  success|trials(trials) ~ 1,   # constant rate
  ~ 1 + distance,
  ~ 1 + distance
)
empty <- mcp(model, data = data_subs, family = binomial(), sample = F,)

prior = list(
  cp_1 = "dirichlet(2)",
  cp_2 = "dirichlet(2)"
)

options(mc.cores = 4)
fit = mcp(model, data = data_subs, prior = prior, family = binomial(), adapt = 700)
plot(fit)

# https://lindeloev.github.io/mcp/articles/priors.html


model = list(
  recalled | trials(items) ~ 1,
  ~ 1,
  ~1
)
prior = list(
  int_1 = "dnorm(2, 1)",  # high recall for easy trials
  cp_1 = "dnorm(4, 1)",  # performance dicontinuity around 4
  items_2 = "dnorm(-0.4, 1) T( , -0.2)"  # performance deteriorates considerably
)

empty = mcp(model, family = binomial(), sample = FALSE)

set.seed(42)
df = data.frame(items = rep(1:9, each = 40))
df$recalled = empty$simulate(
  df$items,
  cp_1 = 4,
  int_1 = 2.5,
  items_2 = -0.4,
  sigma = 0.5)
