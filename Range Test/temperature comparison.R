## Prep data for wind rug as requested by reviewer
library(data.table)


ndbc_dl <- function(url){
  buoy_header <- fread(url,
                       nrows = 2, header = F)
  buoy_header <- buoy_header[, lapply(.SD, function(.) gsub('#', '', .))]
  buoy_header <- buoy_header[, lapply(.SD, function(.) gsub('/', '_', .))]
  buoy_header <- buoy_header[, lapply(.SD, function(.) paste(., collapse = '_'))]

  buoy_data <- fread(url,
                     skip = 2,
                     col.names = as.character(buoy_header))
  buoy_data <- buoy_data[WSPD_m_s != 99,]


  buoy_data <- buoy_data[, date := as.Date(paste(YY_yr, MM_mo, DD_dy, sep = '-'))]
}


coast_2017 <- ndbc_dl('https://www.ndbc.noaa.gov/view_text_file.php?filename=44009h2017.txt.gz&dir=data/historical/stdmet/')
coast_2018 <- ndbc_dl('https://www.ndbc.noaa.gov/view_text_file.php?filename=44009h2018.txt.gz&dir=data/historical/stdmet/')

coast <- rbind(coast_2017, coast_2018)
coast <- coast[date >= '2017-12-21' & date <= '2018-12-04']


coast <- coast[,
               .(mean = mean(WSPD_m_s),
                 max = max(WSPD_m_s)),
               by = date]

coast[, gale := ifelse(mean >= 17.2, T, F)]



# No records after 2018-11-15; supplement with Ocean City Inlet records
ocmd <- ndbc_dl('https://www.ndbc.noaa.gov/view_text_file.php?filename=ocim2h2018.txt.gz&dir=data/historical/stdmet/')

ocmd <- ocmd[date > '2018-11-15' & date <= '2018-12-04']

# Coastal buoy are hourly averages, so do the same here
ocmd <- ocmd[, .(wspd = mean(WSPD_m_s),
                 date = as.Date(paste(YY_yr, MM_mo, DD_dy, sep = '-'))),
             by = c('YY_yr', 'MM_mo', 'DD_dy', 'hh_hr')]

ocmd <- ocmd[,
             .(mean = mean(wspd),
               max = max(wspd)), by = date]
ocmd[, gale := ifelse(mean >= 17.2, T, F)]

coast <- rbind(coast, ocmd)

coast <- coast[mean >= 10.8]

coast[, group := 'Noise~(mV)']
coast[, group := factor(group, levels = c('Temperature~(degree*C)',
                                          "Delta*T~(degree*C)",
                                          'Noise~(mV)'), ordered = T)]



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

  ggnewscale::new_scale_color() +
  geom_rug(data = coast,
           aes(x = date, color = gale), length = unit(0.07, "npc"),
           inherit.aes = T,
           show.legend = F) +
  scale_color_manual(values = c('black', 'red')) +

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



# tiff("range test/manuscript/figures/Figure2_colortemps.tif",
#      width = 170, height = 100, units = 'mm', compression = 'lzw', res = 600,
#      pointsize = 6)
#
#
# ggplot() +
#   geom_hline(data = lims, aes(yintercept = hline),
#              color = "#F0E442") +
#   geom_line(data = data,
#             aes(x = date, y = value, lty = site, color = type)) +
#   scale_color_manual(breaks = c('SST', 'BWT'),
#                      values = c(SST = '#0072B2', NOISE = 'black', DT = 'black', BWT = '#D55E00')) +
#   scale_linetype_manual(values = c('solid', 'dashed')) +
#   facet_wrap(~ group, ncol = 1,
#              labeller = label_parsed,
#              strip.position = 'right', scales = 'free_y') +
#
#   geom_blank(data = lims, aes(y = y)) +
#   scale_y_continuous(breaks = function(.){
#     if(max(.) > 600) c(200, 400, 600, 800) else c(0, 5, 10, 15, 20, 25)
#   }) +
#   scale_x_date(date_breaks = 'month', date_labels = '%b', expand = c(0, 0)) +
#   labs(x = NULL, y = NULL) +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill = NA),
#         strip.text = element_text(size = 10),
#         axis.text = element_text(size = 11),
#         axis.title = element_text(size = 12),
#         legend.title = element_blank(),
#         legend.position = c(0.2, 0.9),
#         legend.box = 'horizontal',
#         legend.margin = margin(0, 2, 2,2),
#         legend.text = element_text(size = 12)) +
#   guides(lty = guide_legend(order = 1))
#
#
# dev.off()