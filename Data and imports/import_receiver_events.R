# Load location information (merging tables is easier than nested ifelse)
rec_deployments <- readxl::read_excel(
  'p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx')

library(lubridate); library(data.table)
rec_deployments <- setDT(rec_deployments)[!is.na(`Cruise ID`) &
                                            !is.na(`Dep VR2AR`)]
rec_deployments <- rec_deployments[, `Dep VR2AR` := paste('VR2AR', `Dep VR2AR`, sep = '-')]
rec_deployments <- rec_deployments[, .(`Cruise ID`, Date, `Site ID`, `Dep VR2AR`,
                                      `Dep Lat_DD`, `Dep Long_DD`)]

# Special case of VR2AR-546462, which was considered lost on 201712 cruise, but
#   found on 201804 cruise
rec_deployments <- rbind(rec_deployments,
                         rec_deployments[grepl('546462', `Dep VR2AR`) &
                                           `Cruise ID` == 201708,])
rec_deployments[133, 1:2] <- list(201712, 20171220)


# Load receiver events, merge with location information
rec_events <- lapply(list.files('p:/obrien/biotelemetry/md wea habitat/data',
                                pattern = 'RecEvents.*.csv',
                                full.names = T),
                     fread, col.names = function(.) tolower(gsub('[) (/]', '', .)))

rec_events <- rbindlist(rec_events, use.names = F)

rec_events <- unique(rec_events, by = c('datetime', 'receiver', 'description'))
rec_events <- rec_events[, datetime := lubridate::ymd_hms(datetime)]
rec_events <- rec_events[, datetime_fover := datetime]
setkey(rec_events, receiver, datetime, datetime_fover)


cruise_id_key <- rec_events[description == 'Data Upload']
cruise_id_key <- cruise_id_key[, cruise := substr(data, 14, 19)]
cruise_id_key <- cruise_id_key[, cruise := fcase(cruise == 201703, 201611,
                                                 cruise == 201708, 201703,
                                                 cruise == 201712, 201708,
                                                 cruise == 201804, 201712,
                                                 cruise == 201808, 201804,
                                                 cruise == 201812, 201808)]
cruise_id_key <- cruise_id_key[, .(min = min(datetime),
                                   max = max(datetime)),
                               by = c('receiver', 'cruise')]

# Special case of VR2AR-546462, which was considered lost on 201712 cruise, but
#   found on 201804 cruise
cruise_id_key[grepl('462$', receiver) &
                cruise == 201712, 'min'] <- lubridate::ymd_hms('2017-12-20 08:38:00')
cruise_id_key <- rbind(cruise_id_key,
                       list('VR2AR-546462', 201708,
                            lubridate::ymd_hms('2017-08-23 12:13:21'),
                            lubridate::ymd_hms('2017-12-20 08:37:59')))

cruise_id_key <- cruise_id_key[order(receiver, min)]


for(i in 2:nrow(cruise_id_key)){
  if(cruise_id_key$min[i] == cruise_id_key$max[i]){
    cruise_id_key$min[i] <- cruise_id_key$max[i - 1] + 1
  }

}
setkey(cruise_id_key, receiver, min, max)

rec_events <- foverlaps(rec_events, cruise_id_key,
                         by.x = c('receiver', 'datetime', 'datetime_fover'),
                         by.y = c('receiver', 'min', 'max'),
                         mult = 'last')
rec_events <- rec_events[, .(receiver, cruise, datetime, description, data, units)]


rec_events <- rec_events[rec_deployments, on = c(receiver = 'Dep VR2AR',
                                                 cruise = 'Cruise ID')]

names(rec_events) <- c('receiver', 'tend.cruise', 'datetime', 'description',
                       'data', 'units', 'c.date', 'site', 'lat', 'long')

# There are issues with the tilt angle of receivers due to deployment differences
# which create different baselines. Calculate the median per receiver per deployment
# period, and subtract it from every value. Then take the absolute value.
tilt <- rec_events[description == 'Tilt angle']
tilt <- tilt[, .(med.tilt = median(as.numeric(data), na.rm = T)),
             by = c('receiver', 'tend.cruise')]
tilt <- rec_events[description == 'Tilt angle'][tilt, on = c('receiver', 'tend.cruise')]
tilt <- tilt[, data := abs(as.numeric(data) - med.tilt)]
tilt <- tilt[, -'med.tilt']

rec_events <- rec_events[description != 'Tilt angle']
rec_events <- rbind(rec_events, tilt)

rec_events <- rec_events[, -'c.date']

saveRDS(rec_events, file = 'data and imports/rec_events.rds')
