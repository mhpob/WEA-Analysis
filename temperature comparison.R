library(ggplot2)
library(tidyr)
library(dplyr)


data <- readRDS('data and imports/rangetest_logit_binary.RDS') %>%
  rename(bwt = `Average temperature`) %>%
  distinct(date, array, .keep_all = T) %>%
  pivot_longer(cols = c(dt, sst, bwt), names_to = 'type') %>%
  mutate(group = ifelse(type == 'dt', "Delta*T", 'Observed'),
         group = factor(group, levels = c('Observed', "Delta*T"), ordered = T),
         type = toupper(type),
         site = ifelse(array == 'Inner', 'Nearshore', 'Offshore'))




ggplot() +
  geom_line(data = data,
            aes(x = date, y = value, lty = site, color = type)) +
  scale_color_manual(breaks = c('SST', 'BWT'),
                     values = c('blue', 'black', 'red')) +
  facet_wrap(~ group, ncol = 1,
             labeller = label_parsed) +
  labs(x = NULL, y = 'Temperature (C)') +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.box = 'horizontal',
        legend.margin = margin(0, 2, 2,2),
        legend.text = element_text(size = 12))
