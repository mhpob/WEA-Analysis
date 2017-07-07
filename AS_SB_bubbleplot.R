library(readxl)
sites <- read_excel('p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx')
sites <- sites[sites$Date == 20161111,]

library(TelemetryR); library(dplyr)
dets <- vemsort('p:/obrien/biotelemetry/detections/offshore md/fish migration')
load('p:/obrien/randomr/ACTactive.rda')

species <- left_join(data.frame(dets), ACTactive,
                     by = c('transmitter' = 'Tag.ID.Code.Standard')) %>%
  mutate(Common.Name = ifelse(grepl('striped', Common.Name, ignore.case = T),
                              'Striped bass', Common.Name))

n_spec_all <- species %>% group_by(station, Common.Name) %>%
  distinct(station, transmitter, .keep_all = T) %>%
  summarize(n = n())

n_spec_month <- species %>%
  mutate(month = lubridate::month(date.utc, label = T, abbr = T),
         month = factor(month, levels = c('Nov', 'Dec', 'Jan', 'Feb', 'Mar'))) %>%
  group_by(station, Common.Name, month) %>%
  distinct(station, transmitter, .keep_all = T) %>%
  summarize(n = n())

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


plot_data_all <- filter(n_spec_all,
               grepl('sturg|striped', Common.Name, ignore.case = T)) %>%
  left_join(sites, by = c('station' = 'Site ID'))

base_map +
  geom_point(data = filter(plot_data_all, grepl('Atl', Common.Name)),
              aes(x = `Dep Long_DD`, y = `Dep Lat_DD`, size = n, color = n)) +
  scale_color_gradient(low = 'black', high = 'red') +
  labs(title = 'Unique Atlantic Sturgeon', color = '', size = '')



plot_data_month <- filter(n_spec_month,
                          grepl('sturg|striped', Common.Name, ignore.case = T)) %>%
  left_join(sites, by = c('station' = 'Site ID'))

base_map +
  geom_point(data = filter(plot_data_month, grepl('Atl', Common.Name)),
             aes(x = `Dep Long_DD`, y = `Dep Lat_DD`, size = n),
             color = 'red') +
  labs(title = 'Unique Atlantic Sturgeon', size = 'Number') +
  facet_wrap(~ month) +
  theme(legend.position = c(0.85, 0.25))
