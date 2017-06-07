library(TelemetryR); library(dplyr)
dec <- vemsort('p:/obrien/biotelemetry/detections/offshore MD/fish migration')

#join ACT data, tidy data frame
load('p:/obrien/randomr/ACTall.rda')
species <- left_join(data.frame(dec), ACTall,
                     by = c('transmitter' = 'Tag.ID.Code.Standard'))

dets <- subset(species,select=-c(4:7,13:16,28:35))

dets[dets$Common.Name %in% 'Striped Bass','Common.Name'] <- 'Striped bass'

detects <- dets[!(dets$Common.Name == "" | is.na(dets$Common.Name)), ]

cal <- function(part){grepl(part, detects[, 'station'])}

detects$Array <- ifelse(cal('A'), 'Array',
                 ifelse(cal('O'), 'Outer','Inner'))


#new explorations
ggplot(data=detects, aes(x=Common.Name)) +
    geom_bar(stat="count")

#V-Track experiments, just striped bass
#Data into V-Track format
library(VTrack)

bass <- filter(detects,Common.Name=="Striped bass")
vbass <- subset(bass,select=c(1,2,3,4))
vbass <- vbass[c(1,3,2,4)]
names(vbass) <- c("DATETIME","TRANSMITTERID","RECEIVERID","STATIONNAME")
vbass$TRANSMITTERID <- do.call(rbind,strsplit(as.character(vbass$TRANSMITTERID), split="-"))[,3]

#V-Track Points and Distance Matrix

weaArray <- unique(bass[,c("station","lat","long")])
names(weaArray) <- c("LOCATION","LATITUDE","LONGITUDE")
weaArray$RADIUS <- 0
row.names(weaArray) <- NULL
plot(weaArray$LONGITUDE, weaArray$LATITUDE)

#SB WEA filtering, movement analysis