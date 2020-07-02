library(TelemetryR); library(dplyr)
dets <- vemsort('p:/obrien/biotelemetry/detections/offshore md/fish migration')

dets <- dets %>%
  mutate(cruise = case_when(date.local <= '2017-03-29' ~ '201703',
                            date.local > '2017-03-29' &
                              date.local <= '2017-08-23' ~ '201708',
                            date.local > '2017-08-23' &
                              date.local <= '2017-12-20' ~ '201712',
                            date.local > '2017-12-20' &
                              date.local <= '2018-04-11' ~ '201804',
                            date.local > '2018-04-11' &
                              date.local <= '2018-08-08' ~ '201808',
                            T ~ '201812'),
         array = case_when(grepl('I', station) ~ 'Inner',
                           grepl('A', station) ~ 'MD WEA',
                           T ~ 'Outer'))


# Detections in wind energy area.
load('p:/obrien/randomr/ACTall.rda')
species <- ACTall %>%
  mutate(Common.Name = tolower(Common.Name),
         Common.Name = case_when(
           grepl('^sturgeon', Common.Name) ~ 'atlantic sturgeon',
           grepl('black sea', Common.Name) ~ 'black sea bass',
           T ~ Common.Name)) %>%
  right_join(dets,
            by = c('Tag.ID.Code.Standard' = 'transmitter'))

# Total detections
sp_overall <- species %>%
  distinct(station, date.local, .keep_all = T) %>%
  group_by(cruise) %>%
  summarize(n = n())

# By array
sp_array <- species %>%
  distinct(station, date.local, .keep_all = T) %>%
  mutate(Common.Name = ifelse(Common.Name == 'not identified', 'atlantic sturgeon',
                              Common.Name)) %>% # assume the not identified are sturgeon
  group_by(cruise, array, Common.Name) %>%
  distinct(Tag.ID.Code.Standard, .keep_all = T) %>%
  summarize(n = n())

# Number of PIs
PIs <- species %>%
  group_by(cruise) %>%
  distinct(cruise, Primary.Researcher) %>%
  #take this out to check the actual number. Skomal is often repeated and NAs are given a category
  summarize(n())

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
