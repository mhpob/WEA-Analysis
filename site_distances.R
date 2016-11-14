library(gdistance)

sites <- read.csv('p:/obrien/biotelemetry/md wea habitat/proposal/sites.csv')
row.names(sites) <- sites$ID
sites <- sites[, c(3, 2)]

midstates <- shapefile('p:/obrien/midatlantic/matl_states_land.shp')

ras.back <- raster(extent(-75.2, -74.4, 38.25, 38.45),
                  resolution = 1/1200, #3 arc-second grid = 1200, 5 = 720, 10 = 360
                  vals = 1,
                  crs = proj4string(midstates))
ras.water <- mask(ras.back, midstates, inverse = T)

plot(ras.water)
points(sites$long, sites$lat)

trans16 <- transition(ras.water, transitionFunction = function(x){1}, 16)
geo16 <- geoCorrection(trans16, type = 'c')

lc.dist <- function (trans, loc, res = c("dist", "path")){
  # Code directly stolen then slightly edited from marmap package
  if (res == "dist") {
    cost <- costDistance(trans, as.matrix(loc))/1000
    return(round(cost, digits = 2))
  }
  if (res == "path") {
    nb.loc <- nrow(loc)
    path <- list()
    comb <- combn(1:nb.loc, 2)
    pb <- txtProgressBar(min = 0, max = ncol(comb), style = 3)
    for (i in 1:ncol(comb)) {
      origin <- sp::SpatialPoints(loc[comb[1, i], ])
      goal <- sp::SpatialPoints(loc[comb[2, i], ])
      temp <- gdistance::shortestPath(trans, origin, goal,
                                      output = "SpatialLines")
      path[[i]] <- temp@lines[[1]]@Lines[[1]]@coords
      setTxtProgressBar(pb, i)
    }
    close(pb)
    return(path)
  }
}

distances <- lc.dist(geo16, sites, res = 'dist')
# Convert km to nautical miles
distances <- 0.5399568 * distances
# Transit times (min) @ 19 kt
transit <- round(distances / 19 * 60)