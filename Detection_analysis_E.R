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


#Species Histogram
ggplot(data=detects, aes(x=Common.Name)) +
    geom_bar(stat="count")

#BOEm Quarterly Updates
fsh <- subset(detects,select=c(2,3,17))
fish<-distinct(fsh)
bass <- filter(fish,Common.Name=="Striped bass")
sbass<-group_by(bass,transmitter)
cbass<-as.tbl(sbass)
tibbass<-cbass %>% summarise(receiver = n())

ggplot(tibbass, aes(receiver))+geom_histogram(binwidth=1,color = 'black')+theme_bw()+labs(x="Number of Receivers", y="Count")

sturg <- filter(fish,Common.Name=="Atlantic sturgeon")
ssturg<-group_by(sturg,transmitter)
csturg<-as.tbl(ssturg)
tibsturg<-csturg %>% summarise(receiver = n())

ggplot(tibsturg, aes(receiver))+geom_histogram(binwidth=1,color = 'black')+theme_bw()+labs(x="Number of Receivers", y="Count")

#by Array
fsh <- subset(detects,select=c(2,3,17,20))
fish<-distinct(fsh)
bass <- filter(fish,Common.Name=="Striped bass")
Ibass <- filter(bass,Array=="Inner")
Isbass<-group_by(Ibass,transmitter)
Icbass<-as.tbl(Isbass)
Itibbass<-Icbass %>% summarise(receiver = n())

ggplot(Itibbass, aes(receiver))+geom_histogram(breaks=seq(1, 13, by=.5))+theme_bw()+labs(x="# Transmitter Detections", y="Count")+ggtitle("Inner Array")

fsh <- subset(detects,select=c(2,3,17,20))
fish<-distinct(fsh)
bass <- filter(fish,Common.Name=="Striped bass")
MDbass <- filter(bass,Array=="MD WEA")
MDsbass<-group_by(MDbass,transmitter)
MDcbass<-as.tbl(MDsbass)
MDtibbass<-MDcbass %>% summarise(receiver = n())

ggplot(MDtibbass, aes(receiver))+geom_histogram(breaks=seq(1, 13, by=.5))+theme_bw()+labs(x="# Transmitter Detections", y="Count")+ggtitle("MD WEA Array")

fsh <- subset(detects,select=c(2,3,17,20))
fish<-distinct(fsh)
bass <- filter(fish,Common.Name=="Striped bass")
Obass <- filter(bass,Array=="Outer")
Osbass<-group_by(Obass,transmitter)
Ocbass<-as.tbl(Osbass)
Otibbass<-Ocbass %>% summarise(receiver = n())

ggplot(Otibbass, aes(receiver))+geom_histogram(breaks=seq(1, 13, by=.5))+theme_bw()+labs(x="# Transmitter Detections", y="Count")+ggtitle("Outer Array")

#Sturgeon Analysis
fsh <- subset(detects,select=c(2,3,17,20))
fish<-distinct(fsh)
sturg <- filter(fish,Common.Name=="Atlantic sturgeon")
Isturg <- filter(sturg,Array=="Inner")
Issturg<-group_by(Isturg,transmitter)
Icsturg<-as.tbl(Issturg)
Itibsturg<-Icsturg %>% summarise(receiver = n())

ggplot(Itibsturg, aes(receiver))+geom_histogram(breaks=seq(1, 5, by=.5))+theme_bw()+labs(x="# Transmitter Detections", y="Count")+ggtitle("Inner Array")

MDsturg <- filter(sturg,Array=="MD WEA")
MDssturg<-group_by(MDsturg,transmitter)
MDcsturg<-as.tbl(MDssturg)
MDtibsturg<-MDcsturg %>% summarise(receiver = n())

ggplot(MDtibsturg, aes(receiver))+geom_histogram(breaks=seq(1, 5, by=.5))+theme_bw()+labs(x="# Transmitter Detections", y="Count")+ggtitle("MD WEA Array")

Osturg <- filter(sturg,Array=="Outer")
Ossturg<-group_by(Osturg,transmitter)
Ocsturg<-as.tbl(Ossturg)
Otibsturg<-Ocsturg %>% summarise(receiver = n())

ggplot(Otibsturg, aes(receiver))+geom_histogram(breaks=seq(1, 5, by=.5))+theme_bw()+labs(x="# Transmitter Detections", y="Count")+ggtitle("Outer Array")

#BOEM report, all species together
fsh <- subset(detects,select=c(2,3,17,20))
fish<-distinct(fsh)
sfish<-group_by(fish,transmitter)
cfish<-as.tbl(sfish)
tibfish<-cfish %>% summarise(receiver = n())

ggplot(tibfish, aes(receiver))+geom_histogram(binwidth=1,color="black")+theme_bw()+labs(x="Number of Receivers", y="Count")+ggtitle("All Species")

Ifish <- filter(fish,Array=="Inner")
Isfish<-group_by(Ifish,transmitter)
Icfish<-as.tbl(Isfish)
Itibfish<-Icfish %>% summarise(receiver = n())

ggplot(Itibfish, aes(receiver))+geom_histogram(breaks=seq(1, 13, by=.5))+theme_bw()+labs(x="Number of Receivers", y="Count")+ggtitle("Inner Array")

MDfish <- filter(fish,Array=="MD WEA")
MDsfish<-group_by(MDfish,transmitter)
MDcfish<-as.tbl(MDsfish)
MDtibfish<-MDcfish %>% summarise(receiver = n())

ggplot(MDtibfish, aes(receiver))+geom_histogram(breaks=seq(1, 13, by=.5))+theme_bw()+labs(x="Number of Receivers", y="Count")+ggtitle("MD WEA Array")

Ofish <- filter(fish,Array=="Outer")
Osfish<-group_by(Ofish,transmitter)
Ocfish<-as.tbl(Osfish)
Otibfish<-Ocfish %>% summarise(receiver = n())

ggplot(Otibfish, aes(receiver))+geom_histogram(breaks=seq(1, 13, by=.5))+theme_bw()+labs(x="Number of Receivers", y="Count")+ggtitle("Outer Array")

#Creating faceted figure
namevector <- c("Array")
Otibfish[,namevector] <- 'Outer'
MDtibfish[,namevector] <- 'MD WEA'
Itibfish[,namevector] <- 'Inner'

big.fish<- rbind(Itibfish,MDtibfish,Otibfish)

ggplot(big.fish, aes(receiver))+geom_histogram(bins = 10, color = 'black')+scale_x_continuous(breaks=c(1,3,5,7,9))+theme_bw()+labs(x="Number of Receivers", y="Count")+facet_wrap(~Array,nrow=3)



## Mike's code (sorry for being a jerk) ----
## Turns out that separating it by array doesn't make sense. Sorry for the goose
## chase. Leaving this here in case you need it for reference.
## All species
# all.spec.counts <- dec %>%
#   group_by(transmitter) %>%
#   summarize(num.detect = n()) %>%
#   filter(num.detect > 1) %>% #make sure that we've heard everything at least 2x
#   left_join(dec) %>%
#   mutate(Array = ifelse(grepl('A', station), 'MD WEA',
#                  ifelse(grepl('O', station), 'Outer','Inner'))) %>%
#   distinct(transmitter, receiver, Array) %>%
#   count(transmitter, Array)
#
# ggplot() +
#   geom_histogram(data = all.spec.counts, aes(x = n), bins = 10, color = 'black',
#                  fill = 'lightgray') +
#   labs(x = 'Number of Receivers', y = 'Number of Transmitters') +
#   scale_x_continuous(breaks = seq(1, 10, 2)) +
#   facet_wrap(~ Array, ncol = 1) +
#   theme_bw()

## Species-specific
# sb.counts <- detects %>%
#   group_by(transmitter) %>%
#   summarize(num.detect = n()) %>%
#   filter(num.detect > 1) %>% #make sure that we've heard everything at least 2x
#   left_join(detects) %>%
#   mutate(Array = ifelse(grepl('A', station), 'MD WEA',
#                  ifelse(grepl('O', station), 'Outer','Inner'))) %>%
#   filter(grepl('striped', Common.Name, ignore.case = T)) %>%
#   distinct(transmitter, receiver, Array) %>%
#   count(transmitter, Array)
#
# ggplot() +
#   geom_histogram(data = sb.counts, aes(x = n), bins = 10, color = 'black',
#                  fill = 'lightgray') +
#   labs(x = 'Number of Receivers', y = 'Number of Transmitters') +
#   scale_x_continuous(breaks = seq(1, 10, 2)) +
#   facet_wrap(~ Array, ncol = 1) +
#   theme_bw()
#
# as.counts <- detects %>%
#   group_by(transmitter) %>%
#   summarize(num.detect = n()) %>%
#   filter(num.detect > 1) %>% #make sure that we've heard everything at least 2x
#   left_join(detects) %>%
#   mutate(Array = ifelse(grepl('A', station), 'MD WEA',
#                         ifelse(grepl('O', station), 'Outer','Inner'))) %>%
#   filter(grepl('sturg', Common.Name, ignore.case = T)) %>%
#   distinct(transmitter, receiver, Array) %>%
#   count(transmitter, Array)
#
# ggplot() +
#   geom_histogram(data = as.counts, aes(x = n), bins = 5, color = 'black',
#                  fill = 'lightgray') +
#   labs(x = 'Number of Receivers', y = 'Number of Transmitters') +
#   scale_x_continuous(breaks = seq(1, 10, 2)) +
#   facet_wrap(~ Array, ncol = 1) +
#   theme_bw()

## The stuff that actually matters... ----
## All species
all.spec.counts <- dec %>%
  group_by(transmitter) %>%
  summarize(num.detect = n()) %>%
  filter(num.detect > 1) %>% #make sure that we've heard everything at least 2x
  left_join(dec) %>%
  distinct(transmitter, receiver) %>%
  count(transmitter)

# dim(all.spec.counts[all.spec.counts$n >= 2,])

ggplot() +
  geom_histogram(data = all.spec.counts, aes(x = n), binwidth = 1,
                 color = 'black',
                 fill = 'lightgray') +
  labs(x = 'Number of Receivers', y = 'Number of Transmitters') +
  scale_x_continuous(breaks = seq(1, max(all.spec.counts$n), 2)) +
  theme_bw()

## Species-specific
#Striped bass
# sb.counts <- detects %>%
#   group_by(transmitter) %>%
#   summarize(num.detect = n()) %>%
#   filter(num.detect > 1) %>% #make sure that we've heard everything at least 2x
#   left_join(detects) %>%
#   filter(grepl('striped', Common.Name, ignore.case = T)) %>%
#   distinct(transmitter, receiver) %>%
#   count(transmitter)
#
# ggplot() +
#   geom_histogram(data = sb.counts, aes(x = n), binwidth = 1, color = 'black',
#                  fill = 'lightgray') +
#   labs(x = 'Number of Receivers', y = 'Number of Transmitters') +
#   scale_x_continuous(breaks = seq(1, max(sb.counts$n), 2)) +
#   theme_bw()

#Sturgeon
# as.counts <- detects %>%
#   group_by(transmitter) %>%
#   summarize(num.detect = n()) %>%
#   filter(num.detect > 1) %>% #make sure that we've heard everything at least 2x
#   left_join(detects) %>%
#   filter(grepl('sturg', Common.Name, ignore.case = T)) %>%
#   distinct(transmitter, receiver) %>%
#   count(transmitter)
#
# ggplot() +
#   geom_histogram(data = as.counts, aes(x = n), binwidth = 1, color = 'black',
#                  fill = 'lightgray') +
#   labs(x = 'Number of Receivers', y = 'Number of Transmitters') +
#   scale_x_continuous(breaks = seq(1, max(as.counts$n), 2)) +
#   theme_bw()
