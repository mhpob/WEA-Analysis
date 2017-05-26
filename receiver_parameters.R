setwd("p:/obrien/biotelemetry/md wea habitat/wea-analysis")
rec <- readRDS("rec_events.rds")
library(dplyr)
glimpse(rec)
View(rec)
rec <- subset(rec,select=-c(9:11))

#Temperature
j <- filter(rec,Description=="Average temperature")
class(j$Data)
j$Data<-as.numeric(j$Data)
hist(j$Data)
ggplot(data=j,aes(x=Date.Time,y=Data,group=Site,colour=Site))+geom_line()+
  geom_point()

#Plot Temp by Site/Location
cal <- function(part){grepl(part, j[, 'Site'])}

j$Array <- ifelse(cal('A'), 'Array',
           ifelse(cal('O'), 'Outer','Inner'))

ggplot(data=j,aes(x=Date.Time,y=Data,group=Array,colour=Array))+
  geom_line()+geom_point()

j$Location <- ifelse(cal('N'), 'North',
              ifelse(cal('S'), 'South','Middle'))

ggplot(data=j,aes(x=Date.Time,y=Data,group=Location,colour=Location))+
  geom_line()+geom_point()