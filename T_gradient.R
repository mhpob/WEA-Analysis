sst <- readRDS('data and imports/sst.rds')
rec_events <- readRDS('data and imports/rec_events.rds')

library(lubridate); library(dplyr)
bwt <- filter(rec_events, Description == 'Average temperature') %>%
  mutate(date = date(Date.Time),
         Data = as.numeric(Data)) %>%
  group_by(date, Site, Lat, Long) %>%
  summarize(bwt = mean(Data)) %>%
  data.frame()
names(bwt) <- tolower(names(bwt))

all <- left_join(bwt, sst, by = c('date', 'site')) %>%
  mutate(dt = sst - bwt,
         array = case_when(grepl('I', site) ~ 'Inner',
                           grepl('A', site) ~ 'MD WEA',
                           grepl('O', site) ~ 'Outer'),
         transect = case_when(grepl('N', site) ~ 'Northern',
                              grepl('M', site) ~ 'Middle',
                              grepl('S', site) ~ 'Soutern'),
         transect = factor(transect, levels = c('Northern', 'Middle', 'Southern')))

ribbon <- all %>%
  filter(!is.na(dt),
         date >= '2016-11-12') %>%
  group_by(date, array) %>%
  summarize(mean = mean(dt),
            max = max(dt),
            min = min(dt))

library(ggplot2)
ggplot(data = ribbon) +
  geom_ribbon(aes(x = date, ymin = min, ymax = max,
                  group = array, fill = array), alpha = 0.8, color = 'gray') +
  labs(x = NULL, y = expression(Delta*T~'(°C)'), fill = NULL) +
  scale_x_date(date_breaks = '2 month', date_labels = '%b-%Y') +
  theme_bw()


