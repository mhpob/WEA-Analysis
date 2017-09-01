library(TelemetryR); library(dplyr)
dets <- vemsort('p:/obrien/biotelemetry/detections/offshore md/fish migration')

# levels(factor(dets$transmitter))

# Number of unique fish detected
n_dets_all <- dets %>%
  summarize(a = n_distinct(transmitter))

# n_dets_201605 <- dets %>%
#   filter(date.local > lubridate::ymd('20160228', tz = "America/New_York"),
#          date.local <= lubridate::ymd('20160517', tz = "America/New_York")) %>%
#   group_by(station) %>%
#   summarize(a = n_distinct(transmitter))

# Detections in wind energy area.
load('p:/obrien/randomr/ACTactive.rda')
species <- left_join(data.frame(dets), ACTactive,
                     by = c('transmitter' = 'Tag.ID.Code.Standard')) %>%
  mutate(Common.Name = ifelse(grepl('striped', Common.Name, ignore.case = T),
                              'Striped bass', Common.Name))

n_spec_all <- species %>% group_by(Common.Name) %>%
  distinct(transmitter, .keep_all = T) %>%
  summarize(n = n())

# Note that the following were double-tagged by K. Dunton:
# 27721/22456
# 27725/22461
# 27731/22465
# 27744/22468

n_spec_201605 <- species %>%
  filter(date.local > lubridate::ymd('20160228', tz = "America/New_York"),
         date.local <= lubridate::ymd('20160517', tz = "America/New_York")) %>%
  group_by(station, Common.Name) %>%
  distinct(station, transmitter) %>%
  summarize(n = n())

# Species, Institution, PI
PI_all <- species %>%
  group_by(Common.Name, Release.Location, Primary.Researcher,
           Primary.Tagging.Organization) %>%
  distinct(transmitter, .keep_all = T) %>%
  summarize(n = n())

PI_201605 <- species %>%
  filter(date.local > lubridate::ymd('20160228', tz = "America/New_York"),
         date.local <= lubridate::ymd('20160517', tz = "America/New_York")) %>%
  group_by(Common.Name, Release.Location, Primary.Researcher,
           Primary.Tagging.Organization) %>%
  distinct(transmitter) %>%
  summarize(n = n())

# Species in wind energy area.
WEA_species <- filter(species, grepl('I|O', station))
levels(factor(WEA_species$Common.Name))
