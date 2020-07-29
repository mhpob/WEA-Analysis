sst <- readRDS('data and imports/sst.rds')
rec_events <- readRDS('data and imports/rec_events.rds')

library(lubridate); library(dplyr)
bwt <- filter(rec_events, description %in% c('Average temperature', 'Temperature')) %>%
  mutate(date = date(datetime),
         data = as.numeric(data)) %>%
  group_by(date, site, lat, long) %>%
  dplyr::summarize(bwt = mean(data)) %>%
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
  labs(x = NULL, y = expression(Delta*T~'(?C)'), fill = NULL) +
  scale_x_date(date_breaks = '2 month', date_labels = '%b%y') +
  geom_line(data = filter(all, site %in% c('IS2', 'AN3')),
            aes(x = date, y = dt, group = site)) +
  theme_bw()+
  theme(text = element_text(size = 22))

