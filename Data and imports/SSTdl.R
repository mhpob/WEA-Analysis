library(readxl)
sites <- read_excel('p:/obrien/biotelemetry/md wea habitat/data/vr2ar deployment_recovery log.xlsx')
sites <- sites[grepl('201712', sites$Date),]

from <- '2016-10-01'
to <- '2018-12-31'

lat <- round(sites$`Dep Lat_DD`, 3)
long <- round(sites$`Dep Long_DD`, 3)


netCDFurls <- paste0(
  'https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.nc?analysed_sst[(',
  from, 'T00:00:00Z):1:(', to, 'T00:00:00Z)][(',
  lat, '):1:(', lat, ')][(', long, '):1:(', long, ')]')

# The next bit can take a while...
# Needed to use download.file(method = 'libcurl') to make this work. The default
# (method = 'wininet') timed out. Note that libcurl is optional on windows; if
# there are issues, use capabilities("libcurl") to make sure libcurl is supported.
for(i in seq_len(nrow(sites))){
  path <- file.path('p:/obrien/biotelemetry/md wea habitat/data/SST',
                    paste(sites$`Site ID`[i], 'nc', sep = '.'))

  download.file(netCDFurls[i],
                destfile = path,
                method = 'libcurl',
                mode = 'wb')
}


library('ncdf4')
sst <- data.frame()
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

saveRDS(sst,
        'p:/obrien/biotelemetry/md wea habitat/wea-analysis/data and imports/sst.rds')
