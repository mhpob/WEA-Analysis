library(readxl); library(TelemetryR); library(dplyr)

detections <- vemsort('p:/obrien/biotelemetry/detections')

boem_sb <- detections %>%
  filter(transmitter %in% paste('A69-9002',
                                seq(6757, 6796, 1),
                                sep = '-')) %>%
  mutate(depth = -1.2129 + sensor.value * 0.3032)

tag_data <- read_excel('p:/obrien/biotelemetry/md wea habitat/data/wea tagging data.xlsx')

boem_sb <- left_join(boem_sb, tag_data, by = c('transmitter' = 'Tag ID'))

saveRDS(boem_sb, file = 'data and imports/boem_sb.rds')
