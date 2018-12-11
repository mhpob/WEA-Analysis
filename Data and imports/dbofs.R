library(ncdf4) # need to use dev build: https://github.com/mdsumner/ncdf4 (20181102)


  nc.data <- nc_open(
    file.path('http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/DBOFS/MODELS/',
              '201806',
              '/nos.dbofs.regulargrid.n006.',
              '20180620',
              '.t12z.nc', fsep = ''))

  long <- ncvar_get(nc.data, 'Longitude',
                    start = c(224, 116),
                    count = c(4, 2))
  lat <- ncvar_get(nc.data, 'Latitude',
                   start = c(224, 116),
                   count = c(4, 2))


  salinity <- ncvar_get(nc.data, 'salt',
                        start = c(224, 116, 1, 1),
                        count = c(4, 2, 10, 1))

#thermocline depth-----

temperature <- ncvar_get(nc.data, 'temp',
                         start = c(224, 116, 1, 1),
                         count = c(4, 2, 10, 1))

library(raster)
ras <- lapply(1:10, function(x) raster(temperature[,,x]))
ras <- stack(ras)
ras <- t(ras)

extent(ras) <- extent(c(range(long), range(lat)))
dep <- ncvar_get(nc.data, 'Depth')

hold <- sapply(1:9, function(x){
  mean(as.matrix(
    (ras[[x + 1]]-ras[[x]])/(dep[x + 1]-dep[x])
  ))
})

dep[which.min(hold)]
----

nc.data <- nc_open(
  file.path('http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/DBOFS/MODELS/',
            '201811',
            '/nos.dbofs.regulargrid.n006.',
            '20181101',
            '.t06z.nc', fsep = ''))

j <- ncvar_get(nc.data, 'Longitude')
k <- ncvar_get(nc.data, 'Latitude')
sa <- ncvar_get(nc.data, 'salt', start = c(224, 116, 1, 1), count = c(4, 2, 10, 1))



unique(which(j <= -74.760 & j >= -74.775, arr.ind = T)[,1])
unique(which(k <= 38.370 & k >= 38.365, arr.ind = T)[,2])



j <- data.frame(
  date = as.POSIXct(ncvar_get(nc.data, 'ocean_time'), origin = '2016-01-01 00:00:00', tz = 'GMT'),

)

nc_close(nc.data)




names(nc.data$var)

dim(ncvar_get(nc.data, 'zeta'))

j <- ncvar_get(nc.data, 'ocean_time')
as.POSIXct(j, origin = '2016-01-01 00:00:00', tz = 'GMT')


j <- data.frame(
  lat = as.vector(ncvar_get(nc.data, 'Latitude')),
  lon = as.vector(ncvar_get(nc.data, 'Longitude')),
  bathy = as.vector(ncvar_get(nc.data, 'h')),
  salt_15 = as.vector(ncvar_get(nc.data, 'salt')[,,8]),
  temp_15 = as.vector(ncvar_get(nc.data, 'temp')[,,8]),
  salt_20 = as.vector(ncvar_get(nc.data, 'salt')[,,10]),
  temp_20 = as.vector(ncvar_get(nc.data, 'temp')[,,10]),
  ssh = as.vector(ncvar_get(nc.data, 'zetatomllw'))
  )

plot(j$lon, j$lat, col = j$bathy, pch = '.')
plot(j$lon, j$lat, col = j$salt_15, pch = '.')
plot(j$lon, j$lat, col = j$salt_20, pch = '.')
plot(j$lon, j$lat, col = j$temp_15, pch = '.')
plot(j$lon, j$lat, col = j$temp_20, pch = '.')
plot(j$lon, j$lat, col = j$ssh, pch = '.')

library(raster)
library(rasterVis)

j <- stack('c:/users/secor/downloads/nos.dbofs.regulargrid.n006.20181031.t18z.nc',
           varname = 'salt',
           bands = 1)


k <- do.call(stack, lapply(1:21,
                           function(x) raster('c:/users/secor/downloads/nos.dbofs.regulargrid.n006.20181031.t18z.nc',
                                              varname = 'salt',
                                              level = x)))

levelplot(k, layers = 7:10)
levelplot(k)


q <- raster(nc.data, varname = 'salt', ndcf = T)
