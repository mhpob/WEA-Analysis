sites <- readxl::read_excel('p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx')
sites <- sites[sites$`Cruise ID` == '201708' & !is.na(sites$`Dep Long_DD`),]

earth.bear <- function (long1, lat1, long2, lat2, radians = F){
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  bear <- atan2(sin(dlon) * cos(b1), cos(a1) * sin(b1) - sin(a1) *
                  cos(b1) * cos(dlon))
  output <- (bear%%(2 * pi))
  if(radians == F){
    output <- output * (180/pi)
  }
  return(output)
}


# s.combn <- data.frame('Site' = t(combn(sites$'Site ID', 2)),
#                       'Long' = t(combn(sites$'Dep Long_DD', 2)),
#                       'Lat' = t(combn(sites$'Dep Lat_DD', 2)))
# s.combn$bearing  <-  earth.bear(s.combn$Long.1, s.combn$Lat.1,
#                                 s.combn$Long.2, s.combn$Lat.2)


combos <- data.frame('Site' = expand.grid(sites$'Site ID', sites$'Site ID',
                                          stringsAsFactors = F),
                      'Long' = expand.grid(sites$'Dep Long_DD', sites$'Dep Long_DD'),
                      'Lat' = expand.grid(sites$'Dep Lat_DD', sites$'Dep Lat_DD'))
combos$bearing  <-  earth.bear(combos$Long.Var1, combos$Lat.Var1,
                               combos$Long.Var2, combos$Lat.Var2, radians = T)
names(combos) <- c('from', 'to', 'long.from', 'long.to', 'lat.from', 'lat.to',
                   'bearing')
