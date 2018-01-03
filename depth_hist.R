library(dplyr); library(ggplot2)

sb <- readRDS('data and imports/boem_sb.rds')

ggplot() + geom_histogram(data = sb, aes(depth), binwidth = 1, na.rm = T) +
  labs(x = 'Depth (m)', y = 'Count') +
  theme_bw()

sb.wea <- filter(sb, grepl('[IAO][NMS]', station))

ggplot() + geom_histogram(data = sb.wea, aes(depth), binwidth = 1, na.rm = T) +
  labs(x = 'Depth (m)', y = 'Count') +
  lims(x = c(0, 25)) +
  theme_bw()
