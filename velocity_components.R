# Load other scripts and then clean up after them.
source('site-to-site bearings.R')
rm(earth.bear)

source('site_distances.R')
rm(geo16, midstates, ras.back, ras.water, sites, trans16, distances, lc.dist)


library(TelemetryR); library(dplyr)
dets <- vemsort('p:/obrien/biotelemetry/detections/offshore md/fish migration')
load('p:/obrien/randomr/ACTactive.rda')

species <- left_join(data.frame(dets), ACTactive,
                     by = c('transmitter' = 'Tag.ID.Code.Standard')) %>%
  mutate(Common.Name = ifelse(grepl('striped', Common.Name, ignore.case = T),
                              'Striped bass',
                         ifelse(grepl('c stur|^stur', Common.Name, ignore.case = T),
                                     'Atlantic sturgeon',
                         ifelse(grepl('white shark',  Common.Name, ignore.case = T),
                                     'White shark',
                                     Common.Name)))) %>%
  filter(Common.Name %in% c('Striped bass', 'Atlantic sturgeon', 'White shark')) %>%
  arrange(date.local)

species_split <- split(species, species$transmitter)
# Make sure we always have at least two detections
species_split <- species_split[sapply(species_split, function(x) dim(x)[1] > 1)]

tdiff_func <- function(x){
  tdiff <- NULL
  for(i in seq(1, nrow(x) - 1, 1)){
    tdiff[i] <- as.double(
      difftime(x$date.local[i + 1], x$date.local[i], units = 'hour')
    )
  }

  from <- x$station[seq(1, nrow(x) - 1, 1)]
  to <- x$station[seq(2, nrow(x), 1)]
  species <- rep(x$Common.Name[1], times = nrow(x) - 1)

  output <- data.frame(species, from, to, tdiff, stringsAsFactors = F)
  output
}

tdiff <- lapply(species_split, tdiff_func)
tdiff <- do.call(rbind, tdiff)

# Select the instances when the fish moves from site to site
tdiff <- tdiff[tdiff$to != tdiff$from,]

site2site <- tdiff %>%
  left_join(combos) %>%
  left_join(tidy_dists) %>%
  mutate(n_s = distance / tdiff * cos(bearing),
         w_e = distance / tdiff * sin(bearing))


library(ggplot2)
ggplot() + geom_histogram(data = site2site,
                     aes(x = n_s)) +
  facet_wrap(~species, ncol = 1, scales = 'free_y') +
  labs(x = 'North-to-South Movement (km/hr)', y = 'Count') +
  theme_bw()

ggplot() + geom_histogram(data = site2site,
                          aes(x = w_e)) +
  facet_wrap(~species, ncol = 1, scales = 'free_y') +
  labs(x = 'West-to-East Movement (km/hr)', y = 'Count') +
  theme_bw()

