setwd("p:/obrien/biotelemetry/md wea habitat/wea-analysis")
rec <- readRDS("rec_events.rds")
library(dplyr)
glimpse(rec)
rec <- subset(rec,select=-c(9:11))

#Temperature
temp<- rec %>%
  filter(Description=="Average temperature",
         Date.Time > '2016-11-13')%>%
  mutate(Data=as.numeric(Data))


#Plot Temp by Site/Location
cal <- function(part){grepl(part, temp[, 'Site'])}

temp$Array <- ifelse(cal('A'), 'MD WEA',
              ifelse(cal('O'), 'Outer','Inner'))

ggplot(data=temp,
       aes(x=Date.Time,y=Data,group=Array,colour=Array))+
  geom_line()+
  geom_point()+
  labs(x="Date", y="Average Temperature")+
  theme_bw()


temp$Transect <- ifelse(cal('N'), 'Northern',
                 ifelse(cal('S'), 'Southern','Middle'))

ggplot(data=temp,
       aes(x=Date.Time,y=Data,group=Transect,colour=Transect))+
  geom_line()+
  geom_point()+
  labs(x="Date", y="Average Temperature")+
  theme_bw()

tempagg <- aggregate(Data ~ lubridate::floor_date(Date.Time, unit = 'day') + Transect, data = temp, FUN = mean)

names(tempagg)[1]<-paste("Date")

tempagg$Transect <- factor(tempagg$Transect,ordered = T,levels=c('Northern','Middle','Southern'))

ggplot(data = tempagg,aes(x = `Date`,y = Data,color = Transect)) +
  geom_point() +
  geom_line()+
  theme_bw()+
  labs(x="Date", y="Average Temperature")


#Noise
k <- filter(rec,Description=="Average noise",Date.Time > '2016-11-13')
k$Data<-as.numeric(k$Data)
ggplot(data=k,aes(x=Date.Time,y=Data,group=Site,colour=Site))+
  geom_line()+geom_point()

#Aggregate by date

noise <- rec %>%
  filter(Description == 'Average noise',Date.Time > '2016-11-13') %>%
  mutate(Data = as.numeric(Data),
         day = lubridate::floor_date(Date.Time, unit = 'day')) %>%
  group_by(day, Site, Lat, Long) %>%
  summarize(avg = mean(Data)) %>%
  as.data.frame()


gal <- function(part){grepl(part, noise[, 'Site'])}
noise$Array <- ifelse(gal('A'), 'MD WEA',
               ifelse(gal('O'), 'Outer','Inner'))

#Aggregate by array

noiseagg <- aggregate(avg ~ day + Array, data = noise, FUN = mean)

ggplot(data=noiseagg,
       aes(x=day,y=avg,group=Array,colour=Array))+
  geom_line()+
  geom_point()+
  labs(x="Date", y="Average Noise")+
  facet_wrap(~Array,nrow=3)+
  theme_bw()


sup <- function(part){grepl(part, noise[, 'Site'])}
noise$Location <- ifelse(sup('N'), 'Northern',
                  ifelse(sup('S'), 'Southern','Middle'))

noise$Transect <- factor(noise$Location,ordered = T,levels=c('Northern','Middle','Southern'))

noiseagg2 <- aggregate(avg ~ day + Transect, data = noise, FUN = mean)

ggplot(data=noiseagg2,
       aes(x=day,y=avg,group=Transect,colour=Transect))+
  geom_line()+
  geom_point()+
  labs(x="Date", y="Average Noise")+
  facet_wrap(~Transect,nrow=3)+
  theme_bw()

#Tilt Analysis
tilt <- rec %>%
  filter(Description == 'Tilt angle',Date.Time > '2016-11-13') %>%
  mutate(Data = as.numeric(Data),
         day = lubridate::floor_date(Date.Time, unit = 'day')) %>%
  group_by(day, Site, Lat, Long) %>%
  summarize(avg = mean(Data)) %>%
  as.data.frame()

gal <- function(part){grepl(part, tilt[, 'Site'])}
tilt$Array <- ifelse(gal('A'), 'MD WEA',
               ifelse(gal('O'), 'Outer','Inner'))

tiltagg <- aggregate(avg ~ day + Array, data = tilt, FUN = mean)

ggplot(data=tiltagg,
       aes(x=day,y=avg,group=Array,colour=Array))+
  geom_line()+
  geom_point()+
  labs(x="Date", y="Tilt Angle")+
  facet_wrap(~Array,nrow=3)+
  theme_bw()

sup <- function(part){grepl(part, tilt[, 'Site'])}
tilt$Location <- ifelse(sup('N'), 'Northern',
                  ifelse(sup('S'), 'Southern','Middle'))

tilt$Transect <- factor(tilt$Location,ordered = T,levels=c('Northern','Middle','Southern'))

tiltagg2 <- aggregate(avg ~ day + Transect, data = tilt, FUN = mean)

ggplot(data=tiltagg2,
       aes(x=day,y=avg,group=Transect,colour=Transect))+
  geom_line()+
  geom_point()+
  labs(x="Date", y="Average Tilt")+
  facet_wrap(~Transect,nrow=3)+
  theme_bw()

#detections

dets <- rec %>%
  filter(grepl("Detections",Description),Date.Time > '2016-11-13') %>%
  mutate(Data = as.numeric(Data),
         day = lubridate::floor_date(Date.Time, unit = 'day')) %>%
  group_by(day, Site, Lat, Long) %>%
  summarize(avg = mean(Data)) %>%
  as.data.frame()

pss <- function(part){grepl(part, dets[, 'Site'])}
dets$Array <- ifelse(pss('A'), 'Array',
               ifelse(pss('O'), 'Outer','Inner'))

ggplot(data=dets,aes(x=day,y=avg,group=Array,colour=Array))+
  geom_line()+geom_point()+labs(x="Date", y="Detections")+facet_wrap(~Array,nrow=3)
