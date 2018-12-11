library(dplyr); library(sf)
coast <- read_sf('c:/users/secor/desktop/gis products/geopackage files/atlcoast.gpkg') %>%
  st_transform(4326)

library(readxl)
sites <- read_excel('p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx') %>%
  select(cruise = `Cruise ID`, site = `Site ID`,
         lat = `Dep Lat_DD`, long = `Dep Long_DD`) %>%
  filter(cruise == 201712,
         grepl('AN3|IS2', site)) %>%
  rbind(data.frame(cruise = c(201712, 201712),
                   site = c(44009, 'OCIM2'),
                   lat = c(38.461, 38.328),
                   long = c(-74.703, -75.091))) %>%
  cbind(type = c(rep('range', 6), rep('env', 2))) %>%
  st_as_sf(coords = c('long', 'lat'), crs = 4326)

library(ggplot2)

ggplot() + geom_sf(data = coast) +
  geom_sf(data = sites, aes(color = type), size = 3) +
  coord_sf(xlim = c(-75.15, -74.7),
           ylim = c(38.3 ,38.5)) +
  guides(color = F) +
  theme_bw()
