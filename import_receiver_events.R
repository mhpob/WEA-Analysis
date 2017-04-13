# Load location information (merging tables is easier than nested ifelse)
rec_deployments <- readxl::read_excel(
  'p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx')

library(dplyr)
rec_deployments <- filter(rec_deployments, Date == '20161111') %>%
  mutate(`Dep VR2AR` = paste('VR2AR', `Dep VR2AR`, sep = '-')) %>%
  select(`Site ID`, `Dep VR2AR`, `Dep Lat_DD`, `Dep Long_DD`)

# Load receiver events, merge with location information
rec_events <- read.csv(
  'p:/obrien/biotelemetry/md wea habitat/data/wea_recevents_201611_201703.csv',
  stringsAsFactors = F)

rec_events$Date.Time <- lubridate::ymd_hms(rec_events$Date.Time)
rec_events <- left_join(rec_events, rec_deployments,
                        by = c('Receiver' = 'Dep VR2AR'))

names(rec_events) <- c('Date.Time', 'Receiver', 'Description', 'Data', 'Units',
                       'Site', 'Lat', 'Long')


saveRDS(rec_events, file = 'rec_events.rds')
