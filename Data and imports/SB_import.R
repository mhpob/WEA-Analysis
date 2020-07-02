library(readxl); library(TelemetryR); library(dplyr)

cl <- parallel::makeCluster(parallel::detectCores() - 1)
detections <- vemsort('p:/obrien/biotelemetry/detections',
                      clust = cl, prog_bar = T,
                      creation_date = '2016-10-31')
parallel::stopCluster(cl)

# Clean up
detections <- subset(detections,select=-c(12))

boem_sb <- detections %>%
  filter(transmitter %in% paste('A69-9002',
                                seq(6757, 6796, 1),
                                sep = '-')) %>%
  mutate(depth = -1.2129 + sensor.value * 0.3032)

tag_data <- read_excel('p:/obrien/biotelemetry/All Transmitters.xlsx')
colnames(tag_data)[9] <- "Common.Name"
tag_data[tag_data$Common.Name %in% 'Striped Bass','Common.Name'] <- 'Striped bass'
sb_data <- filter(tag_data,Common.Name=="Striped bass")

all_sb <- left_join(detections, sb_data,
                     by = c('transmitter' = 'Tag ID Code Standard'))

all_sb <- all_sb[complete.cases(all_sb$Common.Name),]

saveRDS(all_sb, file = 'data and imports/all_sb.rds')
