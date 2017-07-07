library(TelemetryR); library(dplyr)
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

AS <- filter(species, grepl('stur', Common.Name, ignore.case = T))


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
