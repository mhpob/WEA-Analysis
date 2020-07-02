library(ggplot2)
library(tidyr)
library(dplyr)


data <- readRDS('data and imports/rangetest_logit_binary.RDS') %>%
  rename(bwt = `Average temperature`,
         noise = `Average noise`) %>%
  distinct(date, array, .keep_all = T) %>%
  pivot_longer(cols = c(dt, sst, bwt, noise), names_to = 'type') %>%
  mutate(group = case_when(type == 'dt' ~ "Delta*T~(degree*C)",
                           type == 'noise' ~ 'Average~noise~(mV)',
                           T ~ 'Observed~temperature~(degree*C)'),
         group = factor(group, levels = c('Observed~temperature~(degree*C)',
                                          "Delta*T~(degree*C)",
                                          'Average~noise~(mV)'), ordered = T),
         type = toupper(type),
         site = ifelse(array == 'Inner', 'Nearshore', 'Offshore'))




ggplot() +
  geom_line(data = data,
            aes(x = date, y = value, lty = site, color = type)) +
  scale_color_manual(breaks = c('SST', 'BWT'),
                     values = c('blue', 'black', 'black', 'red')) +
  facet_wrap(~ group, ncol = 1,
             labeller = label_parsed,
             strip.position = 'right', scales = 'free_y') +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.box = 'horizontal',
        legend.margin = margin(0, 2, 2,2),
        legend.text = element_text(size = 12))


names(data)
ggplot() +geom_line(data=data, aes(x = date, y = `Average noise`, lty = site))
