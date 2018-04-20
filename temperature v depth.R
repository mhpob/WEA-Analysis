rec.data <- readRDS("data and imports/rec_events.rds")
library(dplyr)
data <- rec.data %>%
  filter(grepl('45[567]$|463$', Receiver),
         Date.Time >= '2017-05-01',
         Date.Time < '2017-12-02',
         Description == 'Average temperature') %>%
  mutate(Data = as.numeric(Data),
         Depth = case_when(
           grepl('456', Receiver) ~ '40 m',
           grepl('463', Receiver) ~ '35 m',
           grepl('455', Receiver) ~ '30 m',
           grepl('457', Receiver) ~ '20 m'
         ),
         Depth = factor(Depth, levels = c('20 m', '30 m', '35 m', '40 m')))

library(ggplot2)
ggplot() + geom_line(data = data, aes(x = Date.Time, y = Data, color = Depth)) +
  labs(x = 'Date', y = 'Hourly Temperature (c)') +
  scale_x_datetime(date_labels = '%b', date_breaks = 'month') +
  theme_bw()

