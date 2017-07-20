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

detects$Array <- ifelse(cal('A'), 'MD WEA',
                 ifelse(cal('O'), 'Outer','Inner'))

#V-Track experiments, just striped bass
#Data into V-Track format
library(VTrack)

bass <- filter(detects,Common.Name=="Striped bass")
vbassn <- subset(bass,select=c(1,2,3,4,5,6,8,11,12,14,16,17))
vbass <- vbassn[c(1,4,5,6,2,7,12,11,3,10,8,9)]

vsbass <- ReadInputData(infile=vbass,
                         fAATAMS=TRUE,
                         fVemcoDualSensor=FALSE,
                         dateformat = "%Y-%m-%d %H:%M:%S")

(TransmitterList <- ExtractUniqueValues(vsbass,2))
(ReceiverList <- ExtractUniqueValues(vsbass,5))

#V-Track Points and Distance Matrix

weaArray <- unique(bass[,c("station","lat","long")])
names(weaArray) <- c("LOCATION","LATITUDE","LONGITUDE")
weaArray$RADIUS <- 0
row.names(weaArray) <- NULL
plot(weaArray$LONGITUDE, weaArray$LATITUDE)
weaArrayDM <- GenerateDirectDistance(weaArray)

#Plotting one trasmitter
par(mfrow=c(1,1),las=1,bty="l")
TID <- ExtractData(vsbass,sQueryTransmitterList = c("53488"))
TID_S <- ExtractUniqueValues(TID,6)# extract station names
plot(as.Date(TID$DATETIME),
     as.numeric(as.factor(as.character(TID$STATIONNAME))),
     ylab="",xlab="DATETIME", type = 's',
     yaxt="n",pch=16,cex.axis=0.9,cex=0.7,
     main=unique(TID[,2]))
axis(side=2,las=1, at=seq(1,length(TID_S),1),cex.axis=0.7,
     labels = TID_S)

TID <- ExtractData(vsbass,sQueryTransmitterList = c("11423":"57598"))

TIDRes<- RunResidenceExtraction(TID,
                                "STATIONNAME",
                                2,
                                60*60*12,
                                sDistanceMatrix=weaArrayDM)

TIDReslog <- TIDRes$residenceslog
TIDResresid <- TIDRes$residences
TIDResmov <- TIDRes$nonresidences

pchDURATION <- ifelse(TIDResresid$DURATION<60,0.1,
                      ifelse(TIDResresid$DURATION<(60*60),0.5,
                             ifelse(TIDResresid$DURATION<(60*60*24),1,3)))
par(mfrow=c(1,1),las=1,bty="l",oma=c(2,1,1,1))
plot(as.Date(TIDResresid$STARTTIME),
 as.numeric(as.factor(as.character(TIDResresid$STATIONNAME))),
 ylab="RECEIVERID",xlab="DATETIME",
 yaxt="n",pch=21,cex=pchDURATION,col="steel blue",
 cex.axis=0.9,main=unique(TIDResresid$TRANSMITTERID))
axis(side=2,las=1, at=seq(1,length(TID_S),1),cex.axis=0.7,
 labels = TID_S)

totalDur <-   tapply(TIDResresid$DURATION,
  TIDResresid$STATIONNAME,sum)
totalDurT <- data.frame(LOCATION =names(totalDur),
  DURATION = as.vector(totalDur))
XYDuration <- merge(weaArray,totalDurT)
pointcex <- XYDuration$DURATION/sum(XYDuration$DURATION)*4

points(weaArray$LONGITUDE,weaArray$LATITUDE,pch=1,cex=0.5,
  col="yellow")
points(XYDuration$LONGITUDE,XYDuration$LATITUDE,
       cex=pointcex, pch=16,col='steel blue')

Vmove_D <- RunTimeProfile(TIDResmov,"STARTTIME","day")
week_mov <- tapply(Vmove_D$DISTANCE,Vmove_D$DATETIME,sum)
plot(as.vector(week_mov)~as.numeric(names(week_mov)),pch=16,
  type='b',xlab="Week number",ylab="Distance moved (km)",
  main="")

#V-Track Points and Distance Matrix

weaArray <- unique(bass[,c("station","lat","long")])
names(weaArray) <- c("LOCATION","LATITUDE","LONGITUDE")
weaArray$RADIUS <- 0
row.names(weaArray) <- NULL
plot(weaArray$LONGITUDE, weaArray$LATITUDE)
weaArrayDM <- GenerateDirectDistance(weaArray)


# To correct, explore further
#
#SB WEA filtering, movement analysis
par(mfrow=c(1,1),bty="l",oma=c(1,1,1,1))
plot(as.Date(vbass$DATETIME), as.numeric(as.factor(as.numeric(as.character(
  vbass$TRANSMITTERID)))),
  ylab="TRANSMITTERID",xlab="DATETIME",
  yaxt="n",pch=16,cex=0.2)
axis(side=2, at=seq(1,length(TransmitterList),1),
     labels=TransmitterList[order(as.numeric(
       TransmitterList))],cex.axis=0.7)

ggplot(data=vbass,aes(x=DATETIME,y=TRANSMITTERID,group=RECEIVERID,colour=RECEIVERID))+geom_point()

#Sturgeon WEA filtering, movement analysis
sturg <- filter(detects,Common.Name=="Atlantic sturgeon")
vsturg <- subset(sturg,select=c(1,2,3,4))
vsturg <- vsturg[c(1,3,2,4)]
names(vsturg) <- c("DATETIME","TRANSMITTERID","RECEIVERID","STATIONNAME")
vsturg$TRANSMITTERID <- do.call(rbind,strsplit(as.character(vsturg$TRANSMITTERID), split="-"))[,3]

ggplot(data=v,aes(x=DATETIME,y=TRANSMITTERID,group=RECEIVERID,colour=RECEIVERID))+geom_point()


#New Things
GenerateAnimationKMLFile_Track(vsbass, # VTrack archive file
TransmitterList[53488], # Transmitter id
weaArray, # spoints file
"Track1.kml", # file name
"cc69deb3") # colour of the track

#
#
#
#
#
#Sturgeon-only analysis

sturg <- filter(detects,Common.Name=="Atlantic sturgeon")
vsturgn <- subset(sturg,select=c(1,2,3,4,5,6,8,11,12,14,16,17))
vsturg <- vsturgn[c(1,4,5,6,2,7,12,11,3,10,8,9)]

vssturg <- ReadInputData(infile=vsturg,
                         fAATAMS=TRUE,
                         fVemcoDualSensor=FALSE,
                         dateformat = "%Y-%m-%d %H:%M:%S")

(TransmitterList <- ExtractUniqueValues(vssturg,2))
(ReceiverList <- ExtractUniqueValues(vssturg,5))

#V-Track Points and Distance Matrix

weaArray <- unique(sturg[,c("station","lat","long")])
names(weaArray) <- c("LOCATION","LATITUDE","LONGITUDE")
weaArray$RADIUS <- 0
row.names(weaArray) <- NULL
plot(weaArray$LONGITUDE, weaArray$LATITUDE)
weaArrayDM <- GenerateDirectDistance(weaArray)

TID <- ExtractData(vssturg,sQueryTransmitterList = c("11182":"57019"))

TIDRes<- RunResidenceExtraction(TID,
                                "STATIONNAME",
                                2,
                                60*60*12,
                                sDistanceMatrix=weaArrayDM)

TIDReslog <- TIDRes$residenceslog
TIDResresid <- TIDRes$residences
TIDResmov <- TIDRes$nonresidences

