library(readxl)
sites <- read_excel('p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx')
sites <- sites[grepl('2017', sites$Date) & !is.na(sites$`Dep Long_DD`),]

j <- expand.grid(1:nrow(sites), 1:nrow(sites))

s.combn <- data.frame('Site' = t(combn(sites$'Site ID', 2)),
                      'Long' = t(combn(sites$'Dep Long_DD', 2)),
                      'Lat' = t(combn(sites$'Dep Lat_DD', 2)),
                      bearing = earth.bear(s.combn$Long.1, s.combn$Lat.1,
                                           s.combn$Long.2, s.combn$Lat.2))



earth.bear <- function (long1, lat1, long2, lat2){
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  bear <- atan2(sin(dlon) * cos(b1), cos(a1) * sin(b1) - sin(a1) *
                  cos(b1) * cos(dlon))
  deg <- (bear%%(2 * pi)) * (180/pi)
  return(deg)
}


