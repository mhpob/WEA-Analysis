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
# lapply(lapply(array.spl, function(x) table(x$internal, x$station)),
#        function(x) round(x/diag(x), 2))

calc_det_freq <- function(data, time_unit){
  spl.data <- split(data, data[time_unit])

  # Count transmission/detection combinations, divide by total transmissions
  spl.data <- lapply(spl.data, function(x) xtabs(~ internal + station, data = x))
  spl.data <- lapply(spl.data, function(x) round(x / diag(x), 2))
  spl.data <- lapply(spl.data, as.data.frame, stringsAsFactors = F)
  for(i in seq_along(spl.data)) spl.data[[i]]$date <- names(spl.data)[i]

  df.data <- do.call(rbind, spl.data)

  df.data <- df.data[df.data$Freq != 1,]
  df.data$d1 <- as.numeric(sapply(strsplit(df.data$internal, '[_]'), '[', 2))
  df.data$d2 <- as.numeric(sapply(strsplit(df.data$station, '[_]'), '[', 2))
  df.data[is.na(df.data)] <- 0
  df.data$distance <- abs(df.data$d1 - df.data$d2)

  row.names(df.data) <- NULL
  df.data <- df.data[, c('date', 'distance', 'Freq')]
  if(time_unit == 'day'){
    df.data$date <- as.Date(df.data$date)
  } else {
    df.data$date <- lubridate::ymd_hms(df.data$date, tz = 'America/New_York')
  }

  df.data <- rbind(df.data,
                   data.frame(date = rep(unique(df.data$date), times = 2),
                              distance = 0,
                              Freq = 1))

  df.data
}

array.spl <- lapply(array.spl, calc_det_freq, time_unit = 'day')
for(i in seq_along(array.spl)) array.spl[[i]]$array <- names(array.spl)[i]
det.freq <- do.call(rbind, array.spl)
row.names(det.freq) <- NULL
det.freq$Freq[det.freq$Freq > 1] <- 1

# Modeling ----
# Logistic regression, fix the asymptote at 1. D50 is distance at 50% detection
wea <- nls(formula = Freq ~ 1 / (1 + exp(-k * (distance - d50))),
           start = list(k = -0.01, d50 = 500),
           data = det.freq, subset = (array == 'MD WEA'))
inner <- nls(formula = Freq ~ 1 / (1 + exp(-k * (distance - d50))),
           start = list(k = -0.01, d50 = 500),
           data = det.freq, subset = (array == 'Inner'))
pred.wea <- predict(wea, data.frame(distance = seq(0, 1000, 50)))
pred.inn <- predict(inner, data.frame(distance = seq(0, 1000, 50)))
predicted <- data.frame(array = rep(c('MD WEA', 'Inner'), each = 21),
                        predicted = c(pred.wea, pred.inn),
                        distance = rep(seq(0, 1000, 50), 2))

# D50 over time
d50_calc <- function(data){
  spl.data <- split(data, data$date)
  spl.data <- lapply(spl.data, function(x){
    tryCatch(
    coef(nls(formula = Freq ~ 1 / (1 + exp(-k * (distance - d50))),
             start = list(k = -0.01, d50 = 500),
             data = x))[2],
    error = function(e) NA)
  })
  if(T %in% is.na(spl.data)) warning('Some D50 unable to be calculated. NA inserted.')
  df.data <- as.data.frame(do.call(rbind, spl.data))
  df.data$date <- as.Date(row.names(df.data))
  row.names(df.data) <- NULL
  df.data
}

array.spl <- split(det.freq, det.freq$array)
d50 <- lapply(array.spl, d50_calc)
for(i in seq_along(d50)) d50[[i]]$array <- names(d50)[i]
d50 <- do.call(rbind, d50)
row.names(d50) <- NULL


# Plotting ----
library(ggplot2)
ggplot(data = det.freq, aes(x = as.Date(date), y = Freq,
                     color = as.factor(distance))) +
  geom_point() +
  stat_smooth() +
  coord_cartesian(ylim = c(0, 1))+
  facet_wrap(~ array) +
  labs(x = NULL, y = 'Frequency of detection', color = 'Distance')

ggplot() +
  geom_hline(aes(yintercept = 0.5)) +
  geom_boxplot(data = det.freq,
               aes(x = distance, y = Freq, group = interaction(distance, array),
                   fill = array)) +
  scale_fill_manual(values = c('#f8766d', '#00ba38')) +
  geom_smooth(data = det.freq, aes(x = distance, y = Freq, lty = array),
              method = 'glm', method.args = list(family = 'binomial')) +
  geom_smooth(data = det.freq, aes(x = distance, y = Freq, lty = array), col = 'green',
              method = 'glm',
              method.args = list(family = binomial(link = 'probit'))) +
  geom_smooth(data = det.freq, aes(x = distance, y = Freq, lty = array), col = 'red',
              method = 'nls',
              method.args = list(formula = y ~ 1 / (1 + exp(-k * (x - d50))),
                                 start = list(k = -0.01, d50 = 500)),
              se = F) +
  # facet_wrap(~ array, nrow = 2) +
  labs(x = 'Distance', y = 'Frequency of detection',
       lty = 'Array', fill = 'Array') +
  theme_bw()

ggplot()+geom_line(data = d50, aes(x = date, y = d50), lwd = 1) +
  facet_wrap(~array, ncol = 1) +
  labs(x = NULL, y = 'D50', color = 'Array') +
  theme_bw()

# Other data ----
# Pairs and TS plot freq by noise, tilt, temperature, sst, deltat (daily)
#     by noise, tilt, temperature (hourly)

rec.data <- readRDS("data and imports/rec_events.rds")
rec.data <- rec.data %>%
  mutate(date.local = .POSIXct(Date.Time, tz = 'America/New_York')) %>%
  filter(date.local > '2017-12-21',
         date.local < '2018-04-11',
         grepl('30[3-9]', Receiver),
         grepl('Average [nt]|Tilt', Description)) %>%
  mutate(array = ifelse(grepl('A', Site), 'MD WEA', 'Inner'),
         Data = as.numeric(Data),
         date = lubridate::date(date.local)) %>%
  group_by(date, array, Description) %>%
  summarize(min = min(Data),
            mean = mean(Data),
            max = max(Data))

freq.cast <- reshape2::dcast(data = det.freq, date + array ~ distance,
                             fun.aggregate = mean, value.var = 'Freq')
rec.cast <- reshape2::dcast(rec.data, date + array ~ Description,
                            fun.aggregate = mean, value.var = 'mean')

env.vars <- full_join(freq.cast, d50) %>% full_join(rec.cast)
pairs(env.vars[env.vars$array == 'MD WEA',4:10])
pairs(env.vars[env.vars$array == 'Inner',4:10])

ggplot() + geom_point(data = env.vars,
                      aes(x = `Average noise`, y = d50, color = array)) +
  labs(x = 'Noise (mV)', y = 'Distance at 50% Detection Probability',
       color = 'Array') +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9))

# SST and DeltaT
sst <- readRDS('data and imports/sst.rds')
sst <- sst %>%
  filter(grepl('IS2|AN3', site)) %>%
  mutate(array = ifelse(grepl('I', site), 'Inner', 'MD WEA')) %>%
  group_by(array, date) %>%
  summarize(sst = mean(sst))

env.vars <- env.vars %>%
  left_join(sst) %>%
  mutate(dt = sst - `Average temperature`)

# Wind/wave direction and magnitude
met.data <- readRDS('data and imports/ndbc_data.rds')
met.data <- met.data %>%
  filter(date.time > '2017-12-21',
         date.time < '2018-04-11',
         station == '44009') %>%
  mutate(date = lubridate::date(date.time)) %>%
  group_by(date, station) %>%
  summarize_all(mean, na.rm = T) %>%
  select(-date.time)

env.vars <- left_join(env.vars, met.data)

ggplot() + geom_point(data = env.vars,
                      aes(x = wvht, y = d50, color = array))
pairs(env.vars[env.vars$array == 'Inner', c(7:10, 12, 14:15, 17:22)])

in.vars <- env.vars[env.vars$array == 'Inner' & !is.na(env.vars$d50),]
wea.vars <- env.vars[env.vars$array == 'MD WEA',]

cor(in.vars$d50, in.vars$wvht)



