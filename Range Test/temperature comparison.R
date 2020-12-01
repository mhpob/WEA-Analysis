library(ggplot2)
library(tidyr)
library(dplyr)


data <- readRDS('data and imports/rangetest_logit_binary_pt0.RDS') %>%
  rename(bwt = average_temperature,
         noise = average_noise) %>%
  distinct(date, array, .keep_all = T) %>%
  pivot_longer(cols = c(dt, sst, bwt, noise), names_to = 'type') %>%
  mutate(group = case_when(type == 'dt' ~ "Delta*T~(degree*C)",
                           type == 'noise' ~ 'Noise~(mV)',
                           T ~ 'Temperature~(degree*C)'),
         group = factor(group, levels = c('Temperature~(degree*C)',
                                          "Delta*T~(degree*C)",
                                          'Noise~(mV)'), ordered = T),
         type = toupper(type),
         site = ifelse(array == 'Inner', 'Nearshore', 'Mid-shelf'))

lims <- data.frame(group = rep(levels(data$group), each = 2),
                   y = c(0, max(data[grepl('SST|BWT', data$type),]$value),
                         0, max(data[data$type == 'DT',]$value),
                         150, max(data[data$type == 'NOISE',]$value)),
                   hline = c(rep(NA, 4), 300, 300))
lims$group <- factor(lims$group, ordered = T, levels = unique(lims$group))


tiff("range test/manuscript/figures/Figure2_colortemps.tif",
     width = 170, height = 100, units = 'mm', compression = 'lzw', res = 600,
     pointsize = 6)


ggplot() +
  geom_hline(data = lims, aes(yintercept = hline),
             color = "#F0E442") +
  geom_line(data = data,
            aes(x = date, y = value, lty = site, color = type)) +
  scale_color_manual(breaks = c('SST', 'BWT'),
                     values = c(SST = '#0072B2', NOISE = 'black', DT = 'black', BWT = '#D55E00')) +
  scale_linetype_manual(values = c('solid', 'dashed')) +
  facet_wrap(~ group, ncol = 1,
             labeller = label_parsed,
             strip.position = 'right', scales = 'free_y') +

  geom_blank(data = lims, aes(y = y)) +
  scale_y_continuous(breaks = function(.){
    if(max(.) > 600) c(200, 400, 600, 800) else c(0, 5, 10, 15, 20, 25)
  }) +
  scale_x_date(date_breaks = 'month', date_labels = '%b', expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.box = 'horizontal',
        legend.margin = margin(0, 2, 2,2),
        legend.text = element_text(size = 12)) +
  guides(lty = guide_legend(order = 1))


dev.off()


tiff("range test/manuscript/figures/Figure2_colorsites.tif",
     width = 170, height = 100, units = 'mm', compression = 'lzw', res = 600,
     pointsize = 8)

ggplot() +
  geom_hline(data = lims, aes(yintercept = hline),
             color = "#F0E442") +
  geom_line(data = data,
            aes(x = date, y = value, color = site, lty = type)) +
  scale_color_manual(values = c('#0072B2','#D55E00')) +
  scale_linetype_manual(breaks = c('BWT', 'SST'),
                        values = c(DT = 'solid', NOISE = 'solid',
                                   BWT = 'solid', SST = 'dotted')) +
  facet_wrap(~ group, ncol = 1,
             labeller = label_parsed,
             strip.position = 'right', scales = 'free_y')+
  geom_blank(data = lims, aes(y = y)) +
  scale_y_continuous(breaks = function(.){
    if(max(.) > 600) c(200, 400, 600, 800) else c(0, 5, 10, 15, 20, 25)
  }) +
  scale_x_date(date_breaks = 'month', date_labels = '%b', expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.box = 'horizontal',
        legend.margin = margin(0, 2, 2,2),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(order = 1))


dev.off()
