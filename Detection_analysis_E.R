library(TelemetryR); library(dplyr)
dec <- vemsort('p:/obrien/biotelemetry/detections/offshore MD/fish migration')

#number unique fish detected, per station
n_dec_all <- dec %>%
  group_by(station) %>%
  summarize(a = n_distinct(transmitter))

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
