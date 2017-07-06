library(dplyr)
library(raster)

WEAinterp <- function(data, coordinates, res = 1){
  if(is.data.frame(data) == F) data <- as.data.frame(data)

  water_qual <- sp::SpatialPointsDataFrame(coords = coordinates,
                                           data = data,
                                           proj4string = sp::CRS('+proj=longlat'))
  water_qual <- sp::spTransform(water_qual,
                                sp::CRS('+proj=utm +zone=18 +datum=NAD83 +units=km'))
  water_qual@data <- water_qual@data %>%
    dplyr::mutate(reasting = res * round(water_qual@coords[, 1] / res),
                  rnorthing = res * round(water_qual@coords[, 2] / res)) %>%
    dplyr::group_by(reasting, rnorthing) %>%
    dplyr::summarize(median = median(data)) %>%
    as.data.frame()
  water_qual@coords <- as.matrix(water_qual@data[, c('reasting','rnorthing')])


  # Make grid to interpolate over, make into spatial object with same
  # projection as shapefile. Use SpatialPixels to more-easily convert to
  # raster later.
  grid <- expand.grid(
    seq(min(water_qual@data$reasting), max(water_qual@data$reasting), res),
    seq(min(water_qual@data$rnorthing), max(water_qual@data$rnorthing), res))
  grid <- sp::SpatialPixelsDataFrame(points = grid[, 1:2], data = grid,
                                     tolerance = 0.99,
                                     proj4string =
                                       sp::CRS('+proj=utm +zone=18 +datum=NAD83 +units=km'))

  wea.ras <- raster::raster(grid, layer = 1)

  ## Create transition matrix that represets a pairwise product of
  ## cells' conductance
  wea.trans <- gdistance::transition(wea.ras, function(x) x[1] * x[2], 8)
  wea.trans <- gdistance::geoCorrection(wea.trans)

  dist <- gdistance::costDistance(wea.trans, water_qual)
  dist <- as.matrix(as.dist(dist, diag = T, upper = T))

  ## Interpolation steps
  # Use transition matrix to calculate corrected distances (as the fish swim)
  pred.dist <- gdistance::costDistance(wea.trans, grid)
  pred.dist <- as.matrix(as.dist(pred.dist, diag = TRUE, upper = TRUE))

  # Calculate distances between observed and predicted values
  op.dist <- gdistance::costDistance(wea.trans, water_qual, grid)
  op.dist <- as.matrix(op.dist, diag = TRUE, upper = TRUE)

  # Fit geostatistical model for the data
  vg <- geoR::variog(coords = water_qual@data[, c('reasting','rnorthing')],
                     data = water_qual@data[, c('reasting','rnorthing', 'median')],
                     max.dis = 600, dists.mat = dist, messages = F)
  # ML fit

    vpar <- SpatialTools::maxlik.cov.sp(
      as.matrix(cbind(1, water_qual@data[, c('reasting','rnorthing')])),
      water_qual@data[, 'median'],
      coords = as.matrix(water_qual@data[, c('reasting','rnorthing')]),
      sp.type = "matern", range.par = 600, error.ratio = 0.5,
      D = dist, reml = T, control = list(trace = 0))

  # Define spatial structure of prediction matrix from fitted spatial model
  V0 <- SpatialTools::cov.sp(coords = as.matrix(
    water_qual@data[, c('reasting','rnorthing')]),
    sp.type = "matern", sp.par = vpar$sp.par,
    error.var = vpar$error.var, smoothness = vpar$smoothness,
    finescale.var = 0,
    pcoords = as.matrix(grid@data[,1:2]),
    D = dist, Dp = pred.dist, Dop = op.dist)

  # Apply spatial structure nd model to predict values
  krige <- SpatialTools::krige.uk(water_qual@data[, 'median'],
                                  V = V0$V, Vop = V0$Vop, Vp = V0$Vp,
                                  X = as.matrix(cbind(1,
                                                      water_qual@data[, c('reasting','rnorthing')])),
                                  Xp = as.matrix(cbind(1, grid@data[,1:2])),
                                  nsim = 0)

  grid[['value']] <- krige$pred
  grid[['se']] <- krige$mspe
  grid
}

midstates <- shapefile('p:/obrien/midatlantic/matl_states_land.shp')
wea <- shapefile('c:/users/secor/desktop/gis products/md mammals/wind_planning_areas/Wind_Planning_Areas_06_20_2014.shp')
midstates <- sp::spTransform(midstates,
                             sp::CRS('+proj=utm +zone=18 +datum=NAD83 +units=km'))
wea <- sp::spTransform(wea,
                       sp::CRS('+proj=utm +zone=18 +datum=NAD83 +units=km'))
library(ggplot2)
midstates <- fortify(midstates)
wea <- fortify(wea)

rec_events <- readRDS('rec_events.rds')



noise <- rec_events %>%
  filter(Description == 'Average noise',
         Date.Time > '2016-11-11') %>%
  mutate(Data = as.numeric(Data),
         Data = ifelse(Data >= 300, 300, Data),
         day = lubridate::floor_date(Date.Time, unit = 'day')) %>%
  group_by(day, Lat, Long) %>%
  summarize(avg = mean(Data)) %>%
  as.data.frame()

dates <- unique(noise$day)

temperature <- rec_events %>%
  filter(Description == 'Average temperature',
         Date.Time > '2016-11-11') %>%
  mutate(Data = as.numeric(Data),
         day = lubridate::floor_date(Date.Time, unit = 'day')) %>%
  group_by(day, Lat, Long) %>%
  summarize(avg = mean(Data)) %>%
  as.data.frame()


# Temperature, noise, tilt interpolation going forward?
# Animations/plots v detections?
library(animation)


saveHTML({
for (i in 1:length(dates)){
  anim.dat <- temperature[temperature$day == dates[i],]
  data <- anim.dat$avg
  coordinates <- matrix(c(anim.dat$Long, anim.dat$Lat), ncol = 2)

  interp.results <- WEAinterp(data, coordinates, res = 1)
  interp.data <- interp.results@data

  test <- ggplot() +
    geom_polygon(data = midstates, aes(x = long, y = lat, group = group),
                 fill  = 'lightgrey', color = 'black') +
    coord_cartesian(xlim = c(490, 545), ylim = c(4230, 4257)) +
    geom_raster(data = interp.data , aes(x = Var1, y = Var2, fill = value)) +
    scale_fill_continuous(limits = c(4, 17)) +
    annotate('text', label = as.Date(dates[i]), x = 540, y = 4255) +
    geom_polygon(data = wea, aes(x = long, y = lat, group = group),
                 fill = NA, color = 'black')

  print(test)
}}, interval = 0.2, ani.height = 720, ani.width = 1280)
