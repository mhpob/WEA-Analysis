library(TelemetryR); library(dplyr)
dets <- vemsort('p:/obrien/biotelemetry/detections/offshore MD/fish migration')
dets <- dets %>%
  group_by(transmitter) %>%
  summarize(num.detect = n()) %>%
  filter(num.detect > 1) %>% #make sure that we've heard everything at least 2x
  left_join(dets)


load('p:/obrien/randomr/ACTactive.rda')

species <- left_join(data.frame(dets), ACTactive,
                     by = c('transmitter' = 'Tag.ID.Code.Standard')) %>%
  mutate(Common.Name = ifelse(grepl('striped', Common.Name, ignore.case = T),
                              'Striped bass', Common.Name))

spl.dets <- split(species, factor(species$transmitter))
spl.tracks <- lapply(spl.dets, track, dates = 'date.utc', ids = 'station')
valid <- lapply(spl.tracks, function(x) dim(x)[1])
spl.tracks <- spl.tracks[valid > 1]

station.movement <- function(x){
  x[, 'end.station'] <- c(x[2 : dim(x)[1], 'station'],
                     NA)
  x[-dim(x)[1],]
}


l <- lapply(spl.tracks, station.movement)
m <- do.call(rbind, l)
n <- as.data.frame(xtabs(~ station + end.station, data = m))


library(circlize)
circos.par(gap.after = c(rep(1, 3), 5, rep(1, 10), 5, rep(1, 3), 5))
chordDiagram(n, directional = 1,
             order = c('IN1', 'IS1', 'IN2', 'IS2', 'AN1', 'AM1', 'AS1', 'AN2',
                       'AM2', 'AS2', 'AN3', 'AM3', 'AS3', 'AN4', 'AM4', 'AS4',
                       'ON1', 'OS1', 'ON2', 'OS2'),
             grid.col = c(IN1 = 'blue', IN2 = 'blue', IS1 = 'blue', IS2 = 'blue',
                          AN1 = 'green', AN2 = 'green', AN3 = 'green',
                          AN4 = 'green', AM1 = 'green', AM2 = 'green',
                          AM3 = 'green', AM4 = 'green', AS1 = 'green',
                          AS2 = 'green', AS3 = 'green', AS4 = 'green',
                          ON1 = 'purple', ON2 = 'purple', OS1 = 'purple',
                          OS2 = 'purple'),
             col = colorRamp2(range(n$Freq), c('yellow', 'red')))
circos.clear()

#North-South movement
NS_move <- function(part){grepl(part, n[, 'station'])}

NSn <- mutate(n, station = ifelse(NS_move('AN'), 'ArrayN',
                           ifelse(NS_move('AM'), 'ArrayM',
                           ifelse(NS_move('AS'), 'ArrayS',
                           ifelse(NS_move('ON'), 'OuterN',
                           ifelse(NS_move('OS'), 'OuterS',
                           ifelse(NS_move('IN'), 'InnerN','InnerS')))))))

NS_move2 <- function(part){grepl(part, NSn[, 'end.station'])}

NSn2 <- mutate(NSn, end.station = ifelse(NS_move2('AN'), 'ArrayN',
                           ifelse(NS_move2('AM'), 'ArrayM',
                           ifelse(NS_move2('AS'), 'ArrayS',
                           ifelse(NS_move2('ON'), 'OuterN',
                           ifelse(NS_move2('OS'), 'OuterS',
                           ifelse(NS_move2('IN'), 'InnerN','InnerS')))))))

circos.par(gap.after = c(rep(1, 1), 5, rep(1, 2), 5, rep(1, 1), 5))
chordDiagram(NSn2, directional = 1,
             order = c('InnerN', 'InnerS',
                       'ArrayN', 'ArrayM', 'ArrayS',
                       'OuterN', 'OuterS'),
             grid.col = c(InnerN = 'blue4', InnerS = 'blue',
                          ArrayN = 'red', ArrayM = 'orange', ArrayS = 'yellow',
                          OuterN = 'purple4', OuterS = 'purple'))

#Striped Bass analysis
bass <- filter(species,Common.Name=="Striped bass")

b.spl.dets <- split(bass, factor(bass$transmitter))
b.spl.tracks <- lapply(b.spl.dets, track, dates = 'date.utc', ids = 'station')
b.valid <- lapply(b.spl.tracks, function(x) dim(x)[1])
b.spl.tracks <- b.spl.tracks[b.valid > 1]

b.station.movement <- function(x){
  x[, 'end.station'] <- c(x[2 : dim(x)[1], 'station'],
                     NA)
  x[-dim(x)[1],]
}


b.l <- lapply(b.spl.tracks, b.station.movement)
b.m <- do.call(rbind, b.l)
b.n <- as.data.frame(xtabs(~ station + end.station, data = b.m))

circos.par(gap.after = c(rep(1, 2), 5, rep(1, 10), 5, rep(1, 3), 5))
chordDiagram(b.n, directional = 1,
             order = c('IN1', 'IS1', 'IN2', 'IS2', 'AN1', 'AM1', 'AS1', 'AN2',
                       'AM2', 'AS2', 'AN3', 'AM3', 'AS3', 'AN4', 'AM4', 'AS4',
                       'ON1', 'OS1', 'ON2', 'OS2'),
             grid.col = c(IN1 = 'blue', IN2 = 'blue', IS1 = 'blue', IS2 = 'blue',
                          AN1 = 'green', AN2 = 'green', AN3 = 'green',
                          AN4 = 'green', AM1 = 'green', AM2 = 'green',
                          AM3 = 'green', AM4 = 'green', AS1 = 'green',
                          AS2 = 'green', AS3 = 'green', AS4 = 'green',
                          ON1 = 'purple', ON2 = 'purple', OS1 = 'purple',
                          OS2 = 'purple'),
             col = colorRamp2(range(b.n$Freq), c('yellow', 'red')))
circos.clear()

#Striped Bass N/S Movements

b.NS_move <- function(part){grepl(part, b.n[, 'station'])}

b.NSn <- mutate(b.n, station = ifelse(b.NS_move('AN'), 'ArrayN',
                           ifelse(b.NS_move('AM'), 'ArrayM',
                           ifelse(b.NS_move('AS'), 'ArrayS',
                           ifelse(b.NS_move('ON'), 'OuterN',
                           ifelse(b.NS_move('OS'), 'OuterS',
                           ifelse(b.NS_move('IN'), 'InnerN','InnerS')))))))

b.NS_move2 <- function(part){grepl(part, b.NSn[, 'end.station'])}

b.NSn2 <- mutate(b.NSn, end.station = ifelse(b.NS_move2('AN'), 'ArrayN',
                           ifelse(b.NS_move2('AM'), 'ArrayM',
                           ifelse(b.NS_move2('AS'), 'ArrayS',
                           ifelse(b.NS_move2('ON'), 'OuterN',
                           ifelse(b.NS_move2('OS'), 'OuterS',
                           ifelse(b.NS_move2('IN'), 'InnerN','InnerS')))))))

circos.par(gap.after = c(rep(1), 5, rep(1, 2), 5, rep(1), 5))
chordDiagram(b.NSn2, directional = 1,
             order = c('InnerN', 'InnerS',
                       'ArrayN', 'ArrayM', 'ArrayS',
                       'OuterN', 'OuterS'),
             grid.col = c(InnerN = 'blue4', InnerS = 'blue',
                          ArrayN = 'red', ArrayM = 'orange', ArrayS = 'yellow',
                          OuterN = 'purple4', OuterS = 'purple'))

#Sturgeon Analysis

sturg <- filter(species,Common.Name=="Atlantic sturgeon")

s.spl.dets <- split(sturg, factor(sturg$transmitter))
s.spl.tracks <- lapply(s.spl.dets, track, dates = 'date.utc', ids = 'station')
s.valid <- lapply(s.spl.tracks, function(x) dim(x)[1])
s.spl.tracks <- s.spl.tracks[s.valid > 1]

s.station.movement <- function(x){
  x[, 'end.station'] <- c(x[2 : dim(x)[1], 'station'],
                     NA)
  x[-dim(x)[1],]
}


s.l <- lapply(s.spl.tracks, s.station.movement)
s.m <- do.call(rbind, s.l)
s.n <- as.data.frame(xtabs(~ station + end.station, data = s.m))

circos.par(gap.after = c(rep(1, 3), 5, rep(1, 10), 5, rep(1), 5))
chordDiagram(s.n, directional = 1,
             order = c('IN1', 'IS1', 'IN2', 'IS2', 'AN1', 'AM1', 'AS1', 'AN2',
                       'AM2', 'AS2', 'AN3', 'AM3', 'AS3', 'AN4', 'AM4', 'AS4',
                       'ON1', 'OS1', 'ON2', 'OS2'),
             grid.col = c(IN1 = 'blue', IN2 = 'blue', IS1 = 'blue', IS2 = 'blue',
                          AN1 = 'green', AN2 = 'green', AN3 = 'green',
                          AN4 = 'green', AM1 = 'green', AM2 = 'green',
                          AM3 = 'green', AM4 = 'green', AS1 = 'green',
                          AS2 = 'green', AS3 = 'green', AS4 = 'green',
                          ON1 = 'purple', ON2 = 'purple', OS1 = 'purple',
                          OS2 = 'purple'),
             col = colorRamp2(range(s.n$Freq), c('yellow', 'red')))
circos.clear()

#Sturgeon N/S Movements

s.NS_move <- function(part){grepl(part, s.n[, 'station'])}

s.NSn <- mutate(s.n, station = ifelse(s.NS_move('AN'), 'ArrayN',
                           ifelse(s.NS_move('AM'), 'ArrayM',
                           ifelse(s.NS_move('AS'), 'ArrayS',
                           ifelse(s.NS_move('ON'), 'OuterN',
                           ifelse(s.NS_move('OS'), 'OuterS',
                           ifelse(s.NS_move('IN'), 'InnerN','InnerS')))))))

s.NS_move2 <- function(part){grepl(part, s.NSn[, 'end.station'])}

s.NSn2 <- mutate(s.NSn, end.station = ifelse(s.NS_move2('AN'), 'ArrayN',
                           ifelse(s.NS_move2('AM'), 'ArrayM',
                           ifelse(s.NS_move2('AS'), 'ArrayS',
                           ifelse(s.NS_move2('ON'), 'OuterN',
                           ifelse(s.NS_move2('OS'), 'OuterS',
                           ifelse(s.NS_move2('IN'), 'InnerN','InnerS')))))))

circos.par(gap.after = c(rep(1), 5, rep(1, 2), 5, rep(1), 5))
chordDiagram(s.NSn2, directional = 1,
             order = c('InnerN', 'InnerS',
                       'ArrayN', 'ArrayM', 'ArrayS',
                       'OuterN', 'OuterS'),
             grid.col = c(InnerN = 'blue4', InnerS = 'blue',
                          ArrayN = 'red', ArrayM = 'orange', ArrayS = 'yellow',
                          OuterN = 'purple4', OuterS = 'purple'))