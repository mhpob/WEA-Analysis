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
                              'Striped bass', Common.Name),
         month = lubridate::month(date.utc, label = T, abbr = T),
         month = factor(month,
                        levels = c('Nov', 'Dec', 'Jan', 'Feb', 'Mar'))) %>%
  group_by(array, Common.Name, transmitter, month) %>%
  distinct(array, transmitter, .keep_all = T) %>%
  summarize(min = min(date.utc))

AS <- filter(species, grepl('stur', Common.Name, ignore.case = T))
SB <- filter(species, grepl('striped', Common.Name, ignore.case = T))

AS <- split(AS, list(AS$array, AS$month))
AS <- AS[lapply(AS, function(x){dim(x)[1]}) > 0]
AS <- lapply(AS, function(x){ecdf(x$min)})
CF_data$date.local <- as_datetime(CF_data$x, tz = 'America/New_York')


final_detect_cf <- function(data, pct.stop){
  new.data <- hi_lo(data)

  #   Note: no built-in function, so I'm going to cheat and use a built-in
  #   function from the ggplot2 package
  CF_plot <- ggplot() + stat_ecdf(data = new.data, aes(x = high)) +
    scale_x_datetime(date_breaks = '3 months')
  print(CF_plot)

  #   Pull data from the plot. Dates are returned in Unix (seconds since
  #   Jan 1, 1970), so we need to convert back to a usable time.
  CF_data <- data.frame(ggplot_build(CF_plot)$data)
  CF_data$date.local <- as_datetime(CF_data$x, tz = 'America/New_York')

  #   Find date where X% of tags stopped transmitting
  stopped <- min(CF_data[CF_data$date.local >= pct.stop, 'date.local'])
  cat((as.numeric(pct.stop) * 100), '% of last detections occurred before',
      as.character(stopped))
}
