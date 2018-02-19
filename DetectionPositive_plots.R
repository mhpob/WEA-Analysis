
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
         year = ifelse(month %in% 10:12 |
                         (month == 9 & lubridate::day(date.local) %in% 16:30),
                       lubridate::year(date.local) - 1999,
                       lubridate::year(date.local) - 2000),
         season = ifelse(month %in% 4:8 |
                           (month == 3 & lubridate::day(date.local) %in% 16:31) |
                           (month == 9 & lubridate::day(date.local) %in% 1:15),
                         'Spring/Summer', 'Autumn/Winter'),
         season = factor(paste0(season, year),
                         levels = c('Autumn/Winter17', 'Spring/Summer17',
                                    'Autumn/Winter18')),
         Common.Name = case_when(
           grepl('striped', Common.Name, ignore.case = T) ~ 'Striped bass',
           grepl('c stur|^stur', Common.Name, ignore.case = T) ~ 'Atlantic sturgeon',
           grepl('white shark',  Common.Name, ignore.case = T) ~ 'White shark',
           grepl('dusky',  Common.Name, ignore.case = T) ~ 'Dusky shark',
           T ~ Common.Name)) %>%
  filter(Common.Name %in%
           c('Striped bass', 'Atlantic sturgeon', 'White shark', 'Dusky shark'))

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
  labs(x = NULL, y = NULL)


base_map +
  geom_point(data = filter(hold, grepl('Dusk', Common.Name)),
             aes(x = `Dep Long_DD`, y = `Dep Lat_DD`, size = fish, color = DP)) +
  facet_wrap(~ season, ncol = 2) +
  scale_color_gradient(low = 'black', high = 'red') +
  labs(title = 'Dusky sharks', color = 'DPD', size = '# Fish') +
  guides(size = guide_legend(order = 1)) +
  theme_bw()
ggsave('AS_DPD.jpeg', path = 'p:/rothermel/boemfinal/DPD', dpi = 300,
       width = 6.5, height = 4, units = 'in', scale = 1.4)

base_map +
  geom_point(data = filter(hold, grepl('bass', Common.Name)),
             aes(x = `Dep Long_DD`, y = `Dep Lat_DD`, size = fish, color = DP)) +
  facet_wrap(~ season, ncol = 2) +
  scale_color_gradient(low = 'black', high = 'red') +
  labs(title = 'Striped bass', color = 'DPD', size = '# Fish') +
  guides(size = guide_legend(order = 1)) +
  theme_bw()
ggsave('SB_DPD.jpeg', path = 'p:/rothermel/boemfinal/DPD', dpi = 300,
       width = 6.5, height = 4, units = 'in', scale = 1.4)

base_map +
  geom_point(data = filter(hold, grepl('White', Common.Name)),
             aes(x = `Dep Long_DD`, y = `Dep Lat_DD`, size = fish, color = DP)) +
  facet_wrap(~ season, ncol = 2) +
  scale_color_gradient(low = 'black', high = 'red') +
  labs(title = 'White sharks', color = 'DPD', size = '# Fish') +
  guides(size = guide_legend(order = 1)) +
  theme_bw()
ggsave('WS_DPD.jpeg', path = 'p:/rothermel/boemfinal/DPD', dpi = 300,
       width = 6.5, height = 4, units = 'in', scale = 1.4)
