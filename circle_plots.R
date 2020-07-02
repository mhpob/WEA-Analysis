library(TelemetryR); library(dplyr)
dets <- vemsort('p:/obrien/biotelemetry/detections/offshore MD/fish migration')
dets <- dets %>%
  group_by(transmitter) %>%
  summarize(num.detect = n()) %>%
  filter(num.detect > 1) %>% #make sure that we've heard everything at least 2x
  left_join(dets)


load('p:/obrien/randomr/ACTall.rda')

species <- left_join(data.frame(dets), ACTall,
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
