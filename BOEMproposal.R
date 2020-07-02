km2lat <- function(km) km / 110.574
km2long <- function(km, lat) km / (111.32 * cos(lat))

# Points from http://www.boem.gov/uploadedFiles/BOEM/Renewable_Energy_Program/State_Activities/Maryland%20Call%20for%20Information.pdf
pt13lat <- 38.379842
pt22lat <- 38.314765
pt23lat <- 38.293135
pt24lat <- 38.293093
pt32lat <- 38.293562

pt10long <- -74.725097
pt11long <- -74.725179
pt18long <- -74.670412
pt19long <- -74.670510
pt30long <- -74.780567

Nlat <- pt13lat - km2lat(0.75 * 1.609344)
Nlong0 <- pt11long + (pt11long - pt10long)
Slat <- pt23lat + 1.5*(pt23lat-pt32lat)/3.5 + km2lat(.75*1.609344)
Slong0 <- pt11long + (pt11long - pt10long)*4

Mlat <- pt13lat - km2lat(3 * 1.609344)
Mlong0 <- pt11long + 2.5*(pt11long - pt10long)

WEA <- data.frame(ID = c('N4', 'N3', 'N2', 'N1',
                         'S4', 'S3', 'S2', 'S1',
                         'M4', 'M3', 'M2', 'M1'),
                  lat = c(rep(Nlat, 4),
                          rep(Slat, 4),
                          rep(Mlat, 4)),
                  long = rep(c(Nlong0, Nlong0 - km2long(3.2, Nlat),
                               Nlong0 - km2long(6.4, Nlat),
                               Nlong0 - km2long(9.6, Nlat)), 3))

IWEA <- data.frame(ID = c('NI2', 'NI1',
                          'SI2', 'SI1'),
                   lat = c(rep(Nlat, 2),
                           rep(Slat, 2)),
                   long =rep(c(WEA[WEA$ID == 'N1', 'long'] - km2long(8, Nlat),
                               WEA[WEA$ID == 'N1', 'long'] - km2long(16, Nlat)),
                             2))

OWEA <- data.frame(ID = c('NO1', 'NO2',
                          'SO1', 'SO2'),
                   lat = c(rep(Nlat, 2),
                           rep(Slat, 2)),
                   long =rep(c(WEA[WEA$ID == 'N4', 'long'] + km2long(8, Nlat),
                               WEA[WEA$ID == 'N4', 'long'] + km2long(16, Nlat)),
                             2))

sites <- rbind(WEA, IWEA, OWEA)

write.csv(sites, 'p:/obrien/biotelemetry/boem proposal/sites2.csv', row.names = F)

TelemetryR::GEcircle(sites$lat, sites$long, 800, 'red', west = F)

## Fish Stuff ----
library(TelemetryR); library(dplyr)
dets <- vemsort('p:/obrien/biotelemetry/detections/offshore')

dets <- filter(dets, grepl('v-|t-|a-', dets$station, ignore.case = T)) 
species <- left_join(data.frame(dets), ACTtrans,
                     by = c('transmitter' = 'Tag.ID.Code.Standard'))
species <- filter(species, Common.Name %in%
                    c("not identified", "Atlantic sturgeon",
                      "Striped Bass", "Striped bass"))
species$names <- ifelse(species$Common.Name %in%
                          c('Atlantic sturgeon', 'not identified'),
                        'Atl. Sturgeon', 'Striped bass')


# test <- species %>%
#   filter(duplicated(trans.num) | duplicated(trans.num, fromLast = TRUE)) %>% 
#   arrange(trans.num, date.local)
# t2 <- NULL
# 
# for(i in 1:(dim(test)[1]-1)){
#   if(test[i, 'trans.num'] != test[i+1, 'trans.num'] |
#      test[i, 'station'] != test[i+1, 'station']){
#     t2 <- rbind(t2, test[i,], test[i + 1,])
#   }
# }
# 
# t3 <- t2 %>% 
#   group_by(trans.num) %>% 
#   summarize(tdiff = difftime(max(date.local), min(date.local), units = 'days')) %>% 
#   filter(tdiff >= 1)
# t4 <- t2 %>% filter(trans.num %in% t3$trans.num)

species$month <- lubridate::month(species$date.local)
sp.ranges <- species %>% group_by(names, lubridate::month(species$date.local),
                                  lubridate::year(species$date.local)) %>%
  summarize(total = n())

## Bubble Figure ----
library(lubridate); library(ggplot2); library(TelemetryR); library(dplyr)
# dets <- vemsort('p:/obrien/biotelemetry/detections/offshore')
# 
# dets <- filter(dets, grepl('v-|t-|a-', dets$station, ignore.case = T)) 
# species <- left_join(data.frame(dets), ACTtrans,
#                      by = c('transmitter' = 'Tag.ID.Code.Standard'))
# species <- filter(species, Common.Name %in%
#                     c("not identified", "Atlantic sturgeon",
#                       "Striped Bass", "Striped bass"))
# species$names <- ifelse(species$Common.Name %in%
#                           c('Atlantic sturgeon', 'not identified'),
#                         'Atl. Sturgeon', 'Striped bass')                      

species$dist <- ifelse(species$station == 'V-2', 5,
                ifelse(species$station == 'V-3', 7,
                ifelse(species$station == 'T-1C', 12,
                ifelse(species$station == 'A-5C', 32,
                ifelse(species$station ==  'T-2C', 50, 65)))))

species$season <- ifelse(species$date.utc > ymd('20141101') &
                           species$date.utc <= ymd('20150228'),
                         'Winter 2015',
                  ifelse(species$date.utc >= ymd('2015-03-01') &
                           species$date.utc <= ymd('2015-05-31'),
                         'Spring 2015',
                  ifelse(species$date.utc >= ymd('2015-06-01') &
                           species$date.utc <= ymd('2015-08-31'),
                         'Summer 2015',
                  ifelse(species$date.utc >= ymd('2015-09-01') &
                           species$date.utc <= ymd('2015-12-31'),
                         'Fall 2015', '!!!'))))
species$season <- factor(species$season, levels = c('Winter 2015', 'Spring 2015',
                                                    'Summer 2015', 'Fall 2015'),
                         ordered = T)

test <- species %>% distinct(season, trans.num, dist) %>% 
  group_by(season, names, dist) %>% summarize(count = n())

proposed.dist <- data.frame(x = c(4.3, 12.3, 20.3, 23.5,
                                  26.7, 29.9, 37.9, 45.9))

ggplot() + geom_point(data = test, aes(x = dist, y = names,
                                       size = count, color = names)) +
  facet_grid(season~.) +
  scale_size(range = c(3, 17),
             breaks = c(1,5,10, 25, 55, 60)) +
  scale_color_manual(values = c('darkred', 'navy'),
                     guide = F) +
  # geom_rug(data = proposed.dist, aes(x = x), lwd = 1) +
  annotate('rect', xmin = 19, xmax = 36, ymin = 0.5, ymax = 2.5,
            color = 'black', fill = NA, lwd = 1) +
  labs(x = 'Distance offshore (km)', y = NULL,
       size = paste('Fish', 'Detected', sep = '\n')) +
  theme_bw()

#### BSB DNR/BOEM ---------------------------------
met <- c(-74 - (45/60) - (13/3600), 38 + (21/60) + (10/3600))
library(rgdal)
blocks <- readOGR('c:/users/secor lab/desktop/gis products/md mammals/wind_planning_areas',
                  'Wind_Planning_Areas_06_20_2014')
midatl <- readOGR('c:/users/secor lab/desktop/gis products/chesapeake/midatlantic',
                  'matl_states_land')
circ_500m <- TelemetryR::ptcirc(met, 500)
circ_1k <- TelemetryR::ptcirc(met, 1000)
circ_2k <- TelemetryR::ptcirc(met, 2000)
circ_4k <- TelemetryR::ptcirc(met, 4000)
circ_20k <- TelemetryR::ptcirc(met, 20000)

library(ggplot2)
blocks <- fortify(blocks)
midatl <- fortify(midatl)

ggplot() + 
  geom_path(data = circ_20k, aes(x = long, y = lat), color= 'green', lwd = 1.5) +
  geom_path(data = circ_4k, aes(x = long, y = lat), color = 'orange', lwd = 1.5) +
  geom_path(data = circ_2k, aes(x = long, y = lat), color = 'orange', lwd = 1.5) +
  geom_path(data = circ_1k, aes(x = long, y = lat), color = 'red', lwd = 1) +
  geom_path(data = circ_500m, aes(x = long, y = lat), color = 'red', lwd = 1) +
  geom_path(data = blocks, aes(long, lat, group = group)) +
  geom_polygon(data = midatl, aes(long, lat, group = group)) +
  coord_map(xlim = c(-75.22, -74.45), ylim = c(38.15, 38.56)) +
  geom_point(aes(x = met[1], y = met[2])) +
  theme_bw() +
  labs(x = 'Longitude', y = 'Latitude')

