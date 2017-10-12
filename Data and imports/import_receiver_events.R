# Load location information (merging tables is easier than nested ifelse)
rec_deployments <- readxl::read_excel(
  'p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx')

library(dplyr)
rec_deployments <- rec_deployments %>%
  filter(`Cruise ID` %in% c(201611, 201703),
         !is.na(`Dep VR2AR`)) %>%
  mutate(`Dep VR2AR` = paste('VR2AR', `Dep VR2AR`, sep = '-'),
         month = lubridate::month(Date)) %>%
  select(`Cruise ID`, Date, `Site ID`, `Dep VR2AR`, `Dep Lat_DD`, `Dep Long_DD`)

# Load receiver events, merge with location information
rec_events_2016 <- read.csv(
  'p:/obrien/biotelemetry/md wea habitat/data/wea_recevents_201611_201703.csv',
  stringsAsFactors = F)
rec_events_2017 <- read.csv(
  'p:/obrien/biotelemetry/md wea habitat/data/wea_recevents_201703_201708.csv',
  stringsAsFactors = F)

rec_events <- rbind(rec_events_2016, rec_events_2017)

rec_events <- rec_events %>%
  distinct(Date.Time, Receiver, Description, .keep_all = T) %>%
  mutate(Date.Time = lubridate::ymd_hms(Date.Time),
         Cruise = as.numeric(ifelse(Date.Time <= '2017-03-30', 201611,
                       ifelse(Date.Time > '2017-03-30' &
                                Date.Time <= '2017-08-27', 201703, '!!!'))))

rec_events <- left_join(rec_events, rec_deployments,
                        by = c('Receiver' = 'Dep VR2AR',
                               'Cruise' = 'Cruise ID'))

names(rec_events) <- c('Date.Time', 'Receiver', 'Description', 'Data', 'Units',
                       'Tend.Cruise', 'C.Date', 'Site', 'Lat', 'Long')

rec_events <- rec_events[, names(rec_events) != 'C.Date']

saveRDS(rec_events, file = 'data and imports/rec_events.rds')
