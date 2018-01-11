library(TelemetryR); library(lubridate); library(dplyr)

## Import detections ----
dets <- vemsort('p:/obrien/biotelemetry/detections/offshore md/fish migration')
dets <- mutate(dets,
               array = ifelse(grepl('I', station), 'Inner',
                              ifelse(grepl('O', station), 'Outer',
                              'MD WEA')))

load('p:/obrien/randomr/ACTactive.rda')
species <- left_join(data.frame(dets), ACTactive,
                     by = c('transmitter' = 'Tag.ID.Code.Standard')) %>%
  mutate(Common.Name = ifelse(grepl('striped', Common.Name, ignore.case = T),
                              'Striped bass',
                       ifelse(grepl('c stur|^stur', Common.Name, ignore.case = T),
                              'Atlantic sturgeon',
                       ifelse(grepl('white shark',  Common.Name, ignore.case = T),
                              'White shark',
                              Common.Name))))

## Import temperature data ----
rec.data <- readRDS("data and imports/rec_events.rds")
rec.data <- rec.data %>%
  filter(Description == 'Average temperature') %>%
  mutate(Date.Time = with_tz(Date.Time, tz = 'America/New_York'),
         Date.Time = floor_date(Date.Time, '6hour'),
         array = ifelse(grepl('I', Site), 'Inner',
                        ifelse(grepl('O', Site), 'Outer',
                               'MD WEA')),
         Data = as.numeric(Data),
         year = lubridate::year(Date.Time) - 2000,
         year = ifelse(month(Date.Time) %in% 10:12,
                       year + 1, year),
         season = ifelse(month(Date.Time) %in% 4:9,
                         'SprSum', 'AutWin'),
         season = paste0(season, year)) %>%
  group_by(Date.Time, season, array) %>%
  summarize(avg.temp = mean(Data))

## Create ECDF plot of detection returns, per array ----
species_all <- species %>%
  group_by(array, Common.Name, transmitter) %>%
  distinct(array, transmitter, .keep_all = T) %>%
  summarize(min = min(date.utc))

det_ecdfplot <- function(spec.plot, array, ylab = NULL, ...){

  data <- filter(species_all, grepl(spec.plot, Common.Name, ignore.case = T))
  data <- split(data, data$array)
  data <- lapply(data, function(x){ecdf(x$min)})

  plot(x = as_datetime(knots(data[[array]])),
       y = data[[array]](knots(data[[array]])),
       ylim = c(0, 1),
       xlim = c(ymd_hms('20161111 19:00:00'),
                ymd_hms('20170328 21:00:00')),
       xlab = 'Date',
       ylab = ifelse(is.null(ylab),
                     eval(expression(paste(
                       'Fraction detected in', array, 'array'))),
                     ylab),
       pch = 16, ...)

  # This is just to stop lines() from throwing warnings when we specify
  # graphical parameters it doesn't like.
  LLines <- function(..., log, axes, frame.plot, panel.first, panel.last) {
    lines(...)
  }
  LLines(x = as_datetime(knots(data[[array]])),
        y = data[[array]](knots(data[[array]])),
        type = 's', ...)
}

det_ecdfplot(spec.plot = 'sturg', array = 'MD WEA')

## Create temperature v ECDF plot ----
temp_ecdfplot <- function(array, spec.data){
  temp.data <- rec.data[rec.data$array == array,]

  par(mar = c(4, 4, 1, 4) + 0.1)
  plot(x = temp.data$Date.Time,
       y = temp.data$avg.temp,
       xlim = c(ymd_hms('20161111 19:00:00'),
                ymd_hms('20170328 21:00:00')),
       ylim = c(4.5, 16.5),
       xaxt = 'n',
       xlab = 'Date',
       ylab = 'Temperature (C)',
       type = 'l',
       lwd = 2)
  axis.POSIXct(1, at = seq(ymd('2016-11-11'),
                        ymd('2017-03-28'),
                        by = 'month'), format = '%m-%Y')

  par(new = T)
  det_ecdfplot(spec.plot = spec.data, array = array, ylab = '', axes = F,
               col = 'blue')
  axis(4, las = 0.5, col.axis = 'blue')
  mtext(side = 4, line = 2, text = eval(expression(paste(
    'Fraction detected in', array, 'array'))), col = 'blue')
}

temp_ecdfplot(array = 'MD WEA', spec.data = 'sturg')


## Same thing, but for season/array combinations ----
# Edit detection data
sp_season <- species %>%
  mutate(year = lubridate::year(date.local) - 2000,
         year = ifelse(lubridate::month(date.local) %in% 10:12,
                       year + 1, year),
         season = ifelse(lubridate::month(date.local) %in% 4:9,
                         'SprSum', 'AutWin'),
         season = paste0(season, year)) %>%
  group_by(array, Common.Name, transmitter, season) %>%
  distinct(array, transmitter, .keep_all = T) %>%
  summarize(min = min(date.local))

AS_season <- filter(sp_season, grepl('c stur', Common.Name, ignore.case = T),
                    season == 'AutWin17')

# Split according to array/season combos
AS_season <- split(AS_season, list(AS_season$array, AS_season$season))

# Select non-empty parts of list
AS_season <- AS_season[lapply(AS_season, function(x){dim(x)[1]}) > 0]

# Apply ecdf to list
AS_season <- lapply(AS_season, function(x){ecdf(x$min)})
plot(AS_season[[1]], verticals = T)

SB <- filter(sp_season, grepl('striped', Common.Name, ignore.case = T))


det_ecdfplot <- function(data, spec.plot, array,
                         season.plot = NULL, ylab = NULL, ...){

  data <- filter(data,
                 grepl(spec.plot, Common.Name, ignore.case = T),
                 season == season.plot)
  if(is.null(season.plot)){
    data <- split(data, data$array)
  } else{
    data <- split(data, list(data$array, data$season))
  }
  data <- lapply(data, function(x){ecdf(x$min)})

  date_lims <- if(grepl('Aut', season.plot)){
    ymd_hms(c(paste0('20', as.numeric(substr(season.plot, 7, 8)) - 1,
                     '1101 00:00:00'),
              paste0('20', substr(season.plot, 7, 8), '0331 23:59:59')),
            tz = 'America/New_York')
  } else {
    ymd_hms(c(paste0('20', substr(season.plot, 7, 8), '0401 00:00:00'),
              paste0('20', substr(season.plot, 7, 8), '0930 23:59:59')),
            tz = 'America/New_York')
  }

  subset <- grep(array, names(data))

  plot(x = as_datetime(knots(data[[subset]])),
       y = data[[subset]](knots(data[[subset]])),
       ylim = c(0, 1),
       xlim = date_lims,
       xlab = 'Date',
       ylab = ifelse(is.null(ylab),
                     eval(expression(paste(
                       'Fraction detected in', array, 'array,', season.plot))),
                     ylab),
       pch = 16, ...)

  # This is just to stop lines() from throwing warnings when we specify
  # graphical parameters it doesn't like.
  LLines <- function(..., log, axes, frame.plot, panel.first, panel.last) {
    lines(...)
  }
  LLines(x = as_datetime(knots(data[[subset]])),
         y = data[[subset]](knots(data[[subset]])),
         type = 's', ...)
}
det_ecdfplot(data = sp_season, 'str', 'MD WEA', season.plot = 'SprSum17')

temp_ecdfplot <- function(spec.data, array, spec.plot, season.plot){
  temp.data <- rec.data[rec.data$array == array &
                          rec.data$season == season.plot,]

  date_lims <- if(grepl('Aut', season.plot)){
    ymd_hms(c(paste0('20', as.numeric(substr(season.plot, 7, 8)) - 1,
                     '1001 00:00:00'),
              paste0('20', substr(season.plot, 7, 8), '0331 23:59:59')),
            tz = 'America/New_York')
  } else {
    ymd_hms(c(paste0('20', substr(season.plot, 7, 8), '0401 00:00:00'),
              paste0('20', substr(season.plot, 7, 8), '0930 23:59:59')),
            tz = 'America/New_York')
  }

  par(mar = c(4, 4, 1, 4) + 0.1)
  plot(x = temp.data$Date.Time,
       y = temp.data$avg.temp,
       xlim = date_lims,
       # ylim = c(4.5, 16.5),
       xaxt = 'n',
       xlab = 'Date',
       ylab = 'Temperature (C)',
       type = 'l',
       lwd = 2)
  axis.POSIXct(1, at = seq(date_lims[1], date_lims[2], by = 'month'),
               format = '%m-%Y')

  par(new = T)
  det_ecdfplot(data = spec.data, spec.plot = spec.plot, array = array,
               season.plot = season.plot, ylab = '', axes = F, col = 'blue')
  axis(4, las = 0.5, col.axis = 'blue')
  mtext(side = 4, line = 2, text = eval(expression(paste(
    'Fraction detected in', array, 'array,', season.plot))), col = 'blue')
}
temp_ecdfplot(spec.data = sp_season, array = 'MD WEA',
              spec.plot = 'str', season.plot = 'AutWin17')
temp_ecdfplot(spec.data = sp_season, array = 'MD WEA',
              spec.plot = 'str', season.plot = 'SprSum17')
temp_ecdfplot(spec.data = sp_season, array = 'MD WEA',
              spec.plot = 'str', season.plot = 'AutWin18')

