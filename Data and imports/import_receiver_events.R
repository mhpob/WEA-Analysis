# Load location information (merging tables is easier than nested ifelse)
rec_deployments <- readxl::read_excel(
  'p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx')

library(lubridate); library(dplyr)
rec_deployments <- rec_deployments %>%
  filter(`Cruise ID` %in% c(201611, 201703, 201708, 201712,
                            201804, 201808, 201812),
         !is.na(`Dep VR2AR`)) %>%
  mutate(`Dep VR2AR` = paste('VR2AR', `Dep VR2AR`, sep = '-')) %>%
  select(`Cruise ID`, Date, `Site ID`, `Dep VR2AR`, `Dep Lat_DD`, `Dep Long_DD`)

# Load receiver events, merge with location information
rec_events_201703 <- read.csv(
  'p:/obrien/biotelemetry/md wea habitat/data/wea_recevents_201611_201703.csv',
  stringsAsFactors = F)
rec_events_201708 <- read.csv(
  'p:/obrien/biotelemetry/md wea habitat/data/wea_recevents_201703_201708.csv',
  stringsAsFactors = F)
rec_events_201712 <- read.csv(
  'p:/obrien/biotelemetry/md wea habitat/data/wea_recevents_201708_201712.csv',
  stringsAsFactors = F)
rec_events_201804 <- read.csv(
  'p:/obrien/biotelemetry/md wea habitat/data/wea_recevents_201712_201804.csv',
  stringsAsFactors = F)
rec_events_201808 <- read.csv(
  'p:/obrien/biotelemetry/md wea habitat/data/wea_recevents_201804_201808.csv',
  stringsAsFactors = F)
rec_events_201812 <- read.csv(
  'p:/obrien/biotelemetry/md wea habitat/data/wea_recevents_201808_201812.csv',
  stringsAsFactors = F)

rec_events <- rbind(rec_events_201703, rec_events_201708, rec_events_201712,
                    rec_events_201804, rec_events_201808, rec_events_201812)

rec_events <- rec_events %>%
  distinct(Date.Time, Receiver, Description, .keep_all = T) %>%
  mutate(Date.Time = lubridate::ymd_hms(Date.Time))

cruise_id_key <- filter(rec_events, Description == 'Data Upload') %>%
  mutate(Cruise = substr(Data, 14, 19),
         Cruise = case_when(Cruise == 201703 ~ 201611,
                            Cruise == 201708 ~ 201703,
                            Cruise == 201712 ~ 201708,
                            Cruise == 201804 ~ 201712,
                            Cruise == 201808 ~ 201804,
                            Cruise == 201812 ~ 201808)) %>%
  group_by(Receiver, Cruise) %>%
  summarize(min = min(Date.Time),
            max = max(Date.Time)) %>%
  arrange(Receiver, min) %>%
  data.frame()

for(i in 2:dim(cruise_id_key)[1]){
  if(cruise_id_key$min[i] == cruise_id_key$max[i]){
    cruise_id_key$min[i] <- cruise_id_key$max[i - 1] + 1
  }

}

rec_events <- left_join(rec_events, cruise_id_key, by = 'Receiver')
rec_events$test <- rec_events$Date.Time %within%
  interval(rec_events$min, rec_events$max)
rec_events <- rec_events[rec_events$test == T, ]
rec_events <- rec_events[, 1:6]


rec_events <- left_join(rec_events, rec_deployments,
                        by = c('Receiver' = 'Dep VR2AR',
                               'Cruise' = 'Cruise ID'))

names(rec_events) <- c('Date.Time', 'Receiver', 'Description', 'Data', 'Units',
                       'Tend.Cruise', 'C.Date', 'Site', 'Lat', 'Long')

rec_events <- rec_events[, names(rec_events) != 'C.Date']

saveRDS(rec_events, file = 'data and imports/rec_events.rds')
