library(gdistance)

sites <- readxl::read_excel('p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx')
sites <- as.data.frame(
  sites[sites$`Cruise ID` == '201804' & !is.na(sites$`Dep Long_DD`),])
row.names(sites) <- sites$`Site ID`
sites <- sites[, c('Dep Long_DD', 'Dep Lat_DD')]
names(sites) <- c('long', 'lat')

bailey_sites <- read.csv('p:/obrien/biotelemetry/md wea habitat/data/bailey_sites.csv',
                         stringsAsFactors = F)
row.names(bailey_sites) <- bailey_sites$Site
bailey_sites <- bailey_sites[, c('Longitude', 'Latitude')]
names(bailey_sites) <- c('long', 'lat')

bsb_sites <- read.csv('p:/obrien/biotelemetry/ocmd-bsb/receiver.csv',
                         stringsAsFactors = F)
row.names(bsb_sites) <- bsb_sites$Site
bsb_sites <- bsb_sites[, c('Long_DD', 'Lat_DD')]
names(bsb_sites) <- c('lat', 'long')

sites <- rbind(data.frame('long' = -75.1033, 'lat' = 38.3274, row.names = 'OCMD'),
               sites, bsb_sites)


ras.back <- raster(extent(-75.2, -74.3, 38.25, 38.45),
                  resolution = 1/1200, #3 arc-second grid = 1200, 5 = 720, 10 = 360
                  vals = 1,
                  crs = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))

plot(ras.back)
points(sites$long, sites$lat)

trans16 <- transition(ras.back, transitionFunction = function(x){1}, 16)
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

# Transit times (min) @ 19 kt (Carson)
# transit <- round(distances / 19 * 60)
# Transit times (min) @ 7 kt (Sea born)
transit <- round(distances / 19 * 60)

write.csv(as.matrix(transit), 'timesAUG2.csv')

# Convert matrix to 3-column data frame
tidy_dists <- reshape2::melt(data = as.matrix(distances),
                             varnames = c('from', 'to'),
                             value.name = 'distance',
                             as.is = T)
