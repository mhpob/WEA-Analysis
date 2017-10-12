library(TelemetryR); library(lubridate);library(dplyr)
dets <- vemsort('p:/obrien/biotelemetry/detections/offshore md/fish migration')
load('p:/obrien/randomr/ACTactive.rda')

elasmo <- left_join(data.frame(dets), ACTactive,
                    by = c('transmitter' = 'Tag.ID.Code.Standard')) %>%
  filter(grepl('shark|ray|skate|blacktip|bonnethead|tiger|hammerhead|dogfish',
               Common.Name, ignore.case = T)) %>%
  mutate(Common.Name = tolower(Common.Name))

plot_data <- elasmo %>%
  filter(Common.Name == 'white shark') %>%
  mutate(mnth = month(date.local, label = T, abbr = T),
         mnth = factor(mnth, levels = c('Nov', 'Dec', 'Jan', 'Feb', 'Mar',
                                         'Apr', 'May', 'Jun', 'Jul', 'Aug')))%>%
  distinct(mnth, transmitter, station, lat, long) %>%
  group_by(mnth, station, lat, long) %>%
  summarize(n = n())

# No detections in January, using a dummy variable to plot empty map
plot_data <- rbind(data.frame(plot_data), data.frame(mnth = 'Jan', station = NA,
                                                     lat = NA, long = NA, n =NA))


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
  geom_point(data = plot_data,
             aes(x = long, y = lat, size = n),
             color = 'red') +
  labs(title = 'Unique Great White Sharks', size = 'Number') +
  facet_wrap(~ mnth, ncol = 4) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.15))
