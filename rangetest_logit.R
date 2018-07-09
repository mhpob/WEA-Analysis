
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

# D50 Modeling ----
lapply(array.spl, function(x){
  glm(cbind(success, fail) ~ distance, family = binomial, data = x)
})

dpct_glm <- function(data, pct = 50){
  spl.data <- split(data, data$date)
  p <- pct/100

  dpct <- lapply(spl.data, function(x){
    # D50 from binomial with logit link
    model_fit <- glm(cbind(success, fail) ~ distance, data = x, family = 'binomial')
    as.numeric(
      (log10(p / (1 - p)) - coef(model_fit)[1]) / coef(model_fit)[2]
    )
  })

  dpct <- data.frame(dpct = do.call(rbind, dpct))
  dpct$date <- as.Date(row.names(dpct))
  row.names(dpct) <- NULL
  dpct[, c('date', 'dpct')]
}

d_probs <- lapply(list(`5` = 5, `25` = 25, `50` = 50, `75` = 75, `95` = 95),
                  function(x){
                    lapply(array.spl, dpct_glm, x) %>%
                      bind_rows(.id = 'array')
                  }) %>%
  bind_rows(.id = 'pct') %>%
  tidyr::spread(pct, dpct)

# Summary ----
# Overall
filter(det.freq, distance != 0) %>%
  group_by(array) %>%
  summarize(mean = mean(freq),
            std = sd(freq),
            min = min(freq),
            max = max(freq))

# Per distance
filter(det.freq, distance != 0) %>%
  group_by(array, distance) %>%
  summarize(mean = mean(freq),
            std = sd(freq),
            min = min(freq),
            max = max(freq))

# D50 summary
d_probs %>%
  group_by(array) %>%
  summarize(mean = mean(`50`),
            std = sd(`50`),
            min = min(`50`),
            max = max(`50`))

# Test array differences, block by distance
model_array <- glm(cbind(success, fail) ~ array + distance, family = binomial, data = det.freq)

# Model curve comparison
model <- glm(cbind(success, fail) ~ distance * array, family = binomial,
             data = det.freq)
drop1(model, ~., test = 'Chisq')



# Plotting ----
library(ggplot2)
ggplot(data = d_probs) +
  geom_ribbon(aes(x = date, ymax = `25`, ymin = `75`), fill = 'grey70') +
  geom_line(aes(x = date, y = `50`), lwd = 1) +
  geom_line(aes(x = date, y = `95`), color = 'palegreen3', lwd = 1) +
  geom_line(aes(x = date, y = `5`), color = 'palevioletred', lwd = 1) +
  labs(x = NULL, y = 'Detection distance (m)') +
  facet_wrap(~ array, nrow = 2) +
  theme_bw()

ggplot(data = d_probs) +
  geom_histogram(aes(`50`, fill = array), binwidth = 25, position = 'identity',
                 alpha = 0.7) +
  theme_bw()

d50plot <- function(array, day){
  array <- enquo(array)
  day <- enquo(day)
  d50dat <- filter(d50, array == !! array, date == !! day)

  ggplot(data = filter(det.freq, array == !! array, date == !! day),
         aes(x = distance, y = freq)) +
    geom_point() +
    geom_smooth(method = 'glm', method.args = list(family = 'binomial'), se = F) +
    geom_point(data = d50dat,
               aes(x = d50, y = 0.5), col = 'red', size = 3) +
    geom_segment(data = d50dat,
                 aes(x = 0, y = 0.5, xend = d50, yend = 0.5)) +
    geom_segment(data = d50dat,
                 aes(x = d50, y = 0, xend = d50, yend = 0.5)) +
    lims(x = c(0, 1000), y = c(0, 1)) +
    labs(x = 'Distance', y = 'Frequency of Detection') +
    coord_flip() +
    theme_bw()
}

# TS modeling ----
library(forecast)
fit <- auto.arima(ts(d_probs[d_probs$array == 'MD WEA', '95']), trace = T)
fit2 <- auto.arima(ts(d_probs[d_probs$array == 'Inner', '95']), trace = T)
tsdiag(fit) # inspect model fit
sqrt(fit$sigma2) #stdev


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
                             fun.aggregate = mean, value.var = 'freq')
rec.cast <- reshape2::dcast(rec.data, date + array ~ Description,
                            fun.aggregate = mean, value.var = 'mean')

env.vars <- full_join(freq.cast, d50) %>% full_join(rec.cast)

GGally::ggpairs(data = env.vars[, c(2, 7, 9:11)], aes(color = array))


ggplot(data = env.vars, aes(x = `Average noise`, y = dpct, color = array)) +
  geom_point(size = 3) +
  geom_smooth(method = 'glm',
              method.args = list(family = Gamma(link = 'inverse')),
              lwd = 1.5) +
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

GGally::ggpairs(data = env.vars, aes(color = array),
                columns = c('dpct', 'Average noise', 'Average temperature', 'Tilt angle',
                  'sst', 'dt', 'wdir', 'wspd', 'gst', 'wvht', 'dpd', 'apd', 'mwd',
                  'pres', 'atmp', 'wtmp'))

# Drop sst, wtmp in favor of Average temperature
# Drop apd and dpd in favor of Average noise
# Drop gst and wvht in favor of wspd
GGally::ggpairs(data = env.vars, aes(color = array),
                columns = c('dpct', 'Average noise', 'Average temperature', 'Tilt angle',
                            'dt', 'wdir', 'wspd', 'mwd', 'pres', 'atmp'))
# Clean up
rm(freq.cast, met.data, rec.cast, sst, i)

# d50 ~ environment modeling ----
mod_data <- env.vars[, c('date', 'array', 'dpct', 'Average noise',
                         'Average temperature', 'Tilt angle',
                         'dt', 'wdir', 'wspd', 'mwd', 'pres', 'atmp')]
mod <- glm(dpct ~ ., family = gaussian(link = 'inverse'),
           data = mod_data[mod_data$array == 'MD WEA',
                           !names(mod_data) %in% c('date', 'array')])

mod_trim <- step(mod)
mod2 <- glm(dpct ~ .^2, family = gaussian(link = 'inverse'),
           data = mod_trim$model)

mod_trim <- step(mod2)

# test <- predict(glm(dpct ~ wspd, family = gaussian(link = 'inverse'),
#                     data = mod_data[mod_data$array == 'MD WEA',
#                                     !names(mod_data) %in% c('date', 'array')]),
#                 list('wspd' = seq(2, 18, 0.2)), type = 'response')
#
# plot(mod_data[mod_data$array == 'MD WEA', 'wspd'], mod_data[mod_data$array == 'MD WEA', 'dpct'])
# lines(seq(2, 18, 0.2), test)
