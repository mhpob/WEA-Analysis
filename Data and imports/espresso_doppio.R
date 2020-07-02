library('ncdf4')
roms <- data.frame()
for(i in seq_along(sites$`Site ID`)){
  path <- file.path('p:/obrien/biotelemetry/md wea habitat/data/SST',
                    paste(sites$`Site ID`[i], 'nc', sep = '.'))
  nc.data <- nc_open(path)
  sst.temp <- data.frame(
    site = strsplit(basename(nc.data$filename), '\\.')[[1]][1],
    date = as.Date.POSIXct(ncvar_get(nc.data, 'time')),
    sst = as.numeric(ncvar_get(nc.data, 'analysed_sst')),
    lat = as.numeric(ncvar_get(nc.data, 'latitude')),
    long = as.numeric(ncvar_get(nc.data, 'longitude')),
    stringsAsFactors = F
  )
  sst <- rbind(sst, sst.temp)
}


nc.data <- nc_open('c:/users/secor/desktop/espresso.nc')
names(nc.data$var)

time <- as.POSIXct(ncvar_get(nc.data, 'time_run')*3600, origin = '2013-05-18 00:00:00', tz = 'UTC')

data <- data.frame(
  date = as.Date.POSIXct(ncvar_get(nc.data, 'time')),
  sst = as.numeric(ncvar_get(nc.data, 'analysed_sst')),
  lat = as.numeric(ncvar_get(nc.data, 'latitude')),
  long = as.numeric(ncvar_get(nc.data, 'longitude')),
  stringsAsFactors = F
)


plot(ncvar_get(nc.data, 'lon_rho'), ncvar_get(nc.data, 'lat_rho'))
points(`lat` ~ `long`, data = sites, col ='red')


'http://tds.marine.rutgers.edu/thredds/ncss/roms/espresso/2013_da/avg/ESPRESSO_Real-Time_v2_Averages_Best?var=h&var=rho&var=salt&var=temp&var=u&var=v&latitude=38.304&longitude=-74.634&disableProjSubset=on&horizStride=1&time_start=2016-11-10T12%3A00%3A00Z&time_end=2018-09-19T12%3A00%3A00Z&timeStride=1&vertCoord=&addLatLon=true&accept=netcdf'



http://tds.marine.rutgers.edu/erddap/tabledap/DOPPIOobs.nc?time%2Cdepth%2Clatitude%2Clongitude%2Cvalue%2Ctype%2Cprovenance&time%3E=2016-11-10&time%3C=2018-12-06&depth%3E=-60&latitude%3E=38.1&latitude%3C=38.5&longitude%3E=-75.1&longitude%3C=-74.3


http://tds.marine.rutgers.edu/thredds/ncss/roms/espresso/2013_da/avg/ESPRESSO_Real-Time_v2_Averages_Best?var=h&var=rho&var=salt&var=temp&var=u&var=v&north=38.5&west=-75.1&east=-74.3&south=38.1&disableProjSubset=on&horizStride=1&time_start=2016-11-10T12%3A00%3A00Z&time_end=2018-09-19T12%3A00%3A00Z&timeStride=1&vertCoord=&addLatLon=true&accept=netcdf

