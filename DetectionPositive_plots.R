
# distinct day/site/transmitter combos
# group by day/site/species
# summarize transmitter (count) by day/site/specie combo
# bubble plot scaled by DPDs
sites <- readxl::read_excel(
  'p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx')
sites <- sites[sites$`Cruise ID` == 201708,]

library(TelemetryR); library(dplyr)
dets <- vemsort('p:/obrien/biotelemetry/detections/offshore md/fish migration')
load('p:/obrien/randomr/ACTactive.rda')

species <- left_join(data.frame(dets), ACTactive,
                     by = c('transmitter' = 'Tag.ID.Code.Standard')) %>%
  mutate(month = lubridate::month(date.local),
    season = ifelse(month %in% 4:6, 'Spring',
                    ifelse(month %in% 7:9, 'Summer',
                    ifelse(month %in% 10:12, 'Autumn',
                           'Winter'))),
    season = ordered(season, levels = c('Autumn', 'Winter', 'Spring', 'Summer')),
    Common.Name = ifelse(grepl('striped', Common.Name, ignore.case = T),
                              'Striped bass',
                       ifelse(grepl('c stur|^stur', Common.Name, ignore.case = T),
                              'Atlantic sturgeon',
                       ifelse(grepl('white shark',  Common.Name, ignore.case = T),
                              'White shark',
                              Common.Name)))) %>%
  filter(Common.Name %in% c('Striped bass', 'Atlantic sturgeon', 'White shark'),
         date.local <= '2017-10-01')

n_unique <- species %>%
  group_by(station, Common.Name, season) %>%
  distinct(station, transmitter, .keep_all = T) %>%
  summarize(fish = n())


hold <- species %>%
  # round to a chosen time group (hour, day, etc.)
  mutate(time.group = lubridate::floor_date(date.local, 'day')) %>%
  # distinct time group/station combos
  distinct(time.group, station, Common.Name, .keep_all = T) %>%
  group_by(station, season, Common.Name) %>%
  summarize(DP = n()) %>%
  left_join(n_unique) %>%
  left_join(sites, by = c('station' = 'Site ID'))


library(raster)
midstates <- shapefile('p:/obrien/midatlantic/matl_states_land.shp')
wea <- shapefile('c:/users/secor/desktop/gis products/md mammals/wind_planning_areas/Wind_Planning_Areas_06_20_2014.shp')

library(ggplot2)
midstates <- fortify(midstates)
wea <- fortify(wea)
base_map <- ggplot() +
  geom_polygon(data = midstates, aes(x = long, y = lat, group = group),
               fill  = 'lightgrey', color = 'black') +
  geom_polygon(data = wea, aes(x = long, y = lat, group = group),
               fill = NA, color = 'black') +
  coord_map(xlim = c(-75.2, -74.45), ylim = c(38.2, 38.5)) +
  labs(x = 'Longitude', y = 'Latitude')


base_map +
  geom_point(data = filter(hold, grepl('stur', Common.Name)),
             aes(x = `Dep Long_DD`, y = `Dep Lat_DD`, size = fish, color = DP)) +
  facet_wrap(~ season, ncol = 2) +
  scale_color_gradient(low = 'black', high = 'red') +
  labs(title = 'Atlantic Sturgeon', color = 'DPD', size = '# Fish') +
  guides(size = guide_legend(order = 1)) +
  theme_bw()

base_map +
  geom_point(data = filter(hold, grepl('bass', Common.Name)),
             aes(x = `Dep Long_DD`, y = `Dep Lat_DD`, size = fish, color = DP)) +
  facet_wrap(~ season, ncol = 2) +
  scale_color_gradient(low = 'black', high = 'red') +
  labs(title = 'Striped Bass', color = 'DPD', size = '# Fish') +
  guides(size = guide_legend(order = 1)) +
  theme_bw()

base_map +
  geom_point(data = filter(hold, grepl('White', Common.Name)),
             aes(x = `Dep Long_DD`, y = `Dep Lat_DD`, size = fish, color = DP)) +
  facet_wrap(~ season, ncol = 2) +
  scale_color_gradient(low = 'black', high = 'red') +
  labs(title = 'White Sharks', color = 'DPD', size = '# Fish') +
  guides(size = guide_legend(order = 1)) +
  theme_bw()