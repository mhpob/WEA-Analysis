library(readxl)
sites <- read_excel('p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx')
sites <- sites[grepl('2017', sites$Date),]

library(raster)
midstates <- shapefile('p:/obrien/midatlantic/matl_states_land.shp')
wea <- shapefile('c:/users/secor/desktop/gis products/md mammals/wind_planning_areas/Wind_Planning_Areas_06_20_2014.shp')

library(ggplot2)
midstates <- fortify(midstates)
wea <- fortify(wea)

ggplot() +
  geom_point(data = sites, aes(x = `Dep Long_DD`, y = `Dep Lat_DD`)) +
  geom_text(data = sites,
            aes(x = `Dep Long_DD`, y = `Dep Lat_DD`,
                    label = `Site ID`, vjust = 1.2)) +
  geom_polygon(data = midstates, aes(x = long, y = lat, group = group),
                        fill  = 'lightgrey', color = 'black') +
  geom_polygon(data = wea, aes(x = long, y = lat, group = group),
               fill = NA, color = 'black') +
  coord_map(xlim = c(-75.2, -74.45), ylim = c(38.2, 38.5)) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme_bw()

