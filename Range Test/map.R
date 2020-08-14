library(dplyr); library(sf)
md_coast <- read_sf('p:/obrien/midatlantic/matl_states_land.shp') %>%
  filter(!is.na(STATE_NAME)) %>%
  st_transform(4326)

atl_coast <- read_sf('data and imports/mapping/natural earth/ne_10m_coastline.shp')%>%
  st_transform(4326)

bathy <- read_sf('data and imports/mapping/bathymetry/test.gpkg',
                 query = 'select * from test where elev_m > -50') %>%
  st_transform(4326) %>%
  filter(as.numeric(st_length(.))> 30000) %>%
  st_crop(xmin = -75.15, xmax = -74.6, ymin = 38.2, ymax = 38.5)

library(readxl)
sites <- read_excel('p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx') %>%
  select(cruise = `Cruise ID`, site = `Site ID`,
         lat = `Dep Lat_DD`, long = `Dep Long_DD`) %>%
  filter(cruise == 201712,
         grepl('AN3|IS|AN2', site)) %>%
  st_as_sf(coords = c('long', 'lat'), crs = 4326)


library(png)
img <- readPNG('range test/manuscript/cartoon_right.png')
arr <- grid::rasterGrob(img, interpolate = TRUE)



library(ggplot2); library(ggspatial)

## Make a one-column map
main <- ggplot() +
  geom_sf(data = bathy, color = 'gray', lwd = 0) +
  geom_label(data = data.frame(
    lab = c('10', '20', '30'),
    x = c(-75.1, -74.838, -74.695),
    y = c(38.27, 38.425, 38.425)
  ),
  aes(x = x, y = y, label = lab),
  label.size = 0, label.padding = unit(0.1, "lines"), size = 1) +
  geom_sf(data = md_coast, size = 0) +
  geom_sf(data = sites, aes(color = ifelse(grepl('A', site), T, F),
                            shape = ifelse(site %in% c('IS1', 'AN2'), T, F)),
          show.legend = F) +
  scale_color_manual(values = c('#D55E00', '#0072B2')) +
  scale_shape_manual(values = c(1, 4)) +
  coord_sf(xlim = c(-75.15, -74.675),
           ylim = c(38.25, 38.45)) +
  annotation_scale(text_cex = 0.5, height = unit(1, 'mm')) +
  annotate('segment', arrow = arrow(length = unit(2, 'mm')),
           x = -74.85, xend = -74.775,
           y = 38.394, yend = 38.375) +
  annotate('segment', arrow = arrow(length = unit(2, 'mm')),
           x = -74.93, xend = -74.93,
           y = 38.35, yend = 38.315) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text = element_text(size = 12 / .pt),
        axis.text.y = element_text(angle = 45, vjust = 0),
        panel.grid = element_line(size = 0),
        axis.ticks = element_line(size = 0),
        plot.margin = unit(c(0, 0, 0, 0), 'mm'))

inset <- ggplotGrob(
  ggplot() +
    geom_sf(data = atl_coast, lwd = 0) +
    coord_sf(xlim = c(-77.5, -69.5),
             ylim = c(35.5, 43)) +
    annotate('rect', xmin = -75.15, xmax = -74.6, ymin = 38.2, ymax = 38.5,
             fill = NA, color = 'red', size  = 0)+
    theme_void() +
    theme(panel.background = element_rect(fill = 'white'))
)



tiff("range test/manuscript/figures/Figure1.tif",
     width = 85, height = 46.7, units = 'mm', compression = 'lzw', res = 600,
     pointsize = 10)


main +
  annotation_custom(inset, xmin = -74.8, xmax = -74.62,
                    ymin = 38.245, ymax = 38.35) +
  annotation_custom(arr, xmin = -75, xmax = -74.85, ymin = 38.35, ymax = 38.45)


dev.off()

###

## Probably need to make a bigger one...