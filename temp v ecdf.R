library(TelemetryR); library(lubridate); library(dplyr)
dets <- vemsort('p:/obrien/biotelemetry/detections/offshore md/fish migration')
dets <- mutate(dets,
               array = ifelse(grepl('I', station), 'Inner',
                              ifelse(grepl('O', station), 'Outer',
                              'Array')))

load('p:/obrien/randomr/ACTactive.rda')
species <- left_join(data.frame(dets), ACTactive,
                     by = c('transmitter' = 'Tag.ID.Code.Standard')) %>%
  mutate(Common.Name = ifelse(grepl('striped', Common.Name, ignore.case = T),
                              'Striped bass', Common.Name)) %>%
  group_by(array, Common.Name, transmitter) %>%
  distinct(array, transmitter, .keep_all = T) %>%
  summarize(min = min(date.utc))

## Create ECDF plot of detection returns, per array ----
det_ecdfplot <- function(spec.plot, array, ...){

  data <- filter(species, grepl(spec.plot, Common.Name, ignore.case = T))
  data <- split(data, data$array)
  data <- lapply(data, function(x){ecdf(x$min)})

  plot(x = as_datetime(knots(data[[array]])),
       y = data[[array]](knots(data[[array]])),
       ylim = c(0, 1),
       xlim = c(ymd_hms('20161111 19:00:00'),
                ymd_hms('20170328 21:00:00')),
       xlab = 'Date',
       ylab = array,
       pch = 16, ...)
  lines(x = as_datetime(knots(data[[array]])),
        y = data[[array]](knots(data[[array]])),
        type = 's', ...)
}

det_ecdfplot(spec.plot = 'sturg', array = 'Array')
det_ecdfplot(spec.plot = 'sturg', array = 'Outer')
det_ecdfplot(spec.plot = 'sturg', array = 'Inner', col = 'blue')

## Bring in temperature data ----
rec.data <- readRDS("rec_events.rds")
rec.data <- rec.data %>%
  filter(Description == 'Average temperature',
         Date.Time > ymd_hms('20161111 19:00:00')) %>%
  mutate(array = ifelse(grepl('I', Site), 'Inner',
                        ifelse(grepl('O', Site), 'Outer',
                               'Array')),
         Data = as.numeric(Data)) %>%
  group_by(Date.Time, array) %>%
  summarize(avg.temp = mean(Data))

temp_ecdfplot <- function(array, spec.data, ...){
  temp.data <- rec.data[rec.data$array == array,]

  plot(x = temp.data$Date.Time,
       y = temp.data$avg.temp,
       # ylim = c(0, 1),
       xlim = c(ymd_hms('20161111 19:00:00'),
                ymd_hms('20170328 21:00:00')),
       xlab = 'Date',
       ylab = 'Temperature (C)',
       type = 'l')
  par(new = T)
  det_ecdf <- det_ecdfplot(spec.plot = 'sturg', array = 'Array', ...)
}



## Same thing, but for month/array combinations ----
species_mo <- species %>%
  mutate(month = lubridate::month(date.utc, label = T, abbr = T),
         month = factor(month,
                        levels = c('Nov', 'Dec', 'Jan', 'Feb', 'Mar'))) %>%
  group_by(array, Common.Name, transmitter, month) %>%
  distinct(array, transmitter, .keep_all = T) %>%
  summarize(min = min(date.utc))

AS_mo <- filter(species_mo, grepl('stur', Common.Name, ignore.case = T))

# Split according to array/month combos
AS_mo <- split(AS_mo, list(AS_mo$array, AS_mo$month))

# Select non-empty parts of list
AS_mo <- AS_mo[lapply(AS_mo, function(x){dim(x)[1]}) > 0]

# Apply ecdf to list
AS_mo <- lapply(AS_mo, function(x){ecdf(x$min)})
# plot(AS_mo[[2]], verticals = T)

SB <- filter(species_mo, grepl('striped', Common.Name, ignore.case = T))


# Dates are returned in Unix (seconds since Jan 1, 1970), so we need to convert
# back to a usable time. Will be something like:
# data$date <- as_datetime(data$date, tz = 'America/New_York')
