# Packages ----
library(ggplot2); library(tidyr); library(data.table); library(mgcv)



# Import data ----
data <- readRDS('data and imports/rangetest_median_sitespec.RDS')
data <- data[data$distance > 0 & data$distance < 2400,]
names(data) <- gsub(' ', '_', tolower(names(data)))
data$array <- as.factor(gsub(' ', '', data$array))
data$station <- as.factor(data$station)
data$freq <- data$success / (data$success + data$fail)
data$wgt <- (data$success + data$fail) / mean((data$success + data$fail))
data$date <- as.factor(data$date)

data <- data[order(data$station, data$distance, data$date),]
data$ts.start <- ifelse(data$date == '2017-12-21' |
                          (data$date == '2018-08-08' & data$station == 'AN3' & data$distance == 800) |
                          (data$date == '2018-08-08' & data$station == 'AN3_250' & data$distance == 550), T, F)



# Import winning model ----
##  see "mgcv model selection.R" for selection

m_ndt <- bam(freq ~ distance + s(station, array, bs = 're') +
               s(average_noise, k = 40) +
               s(dt, k = 40),
             family = quasibinomial(),
             data = data,
             weights = data$wgt,
             discrete = T,
             rho = 0.466396825956338,
             AR.start = data$ts.start)

med_noise <- unique(data.table(data), by = c('date', 'station'))[, median(average_noise)]


newdata <- data.frame(
  array = 'MDWEA',
  station = 'AN3',
  dt = seq(min(data$dt), max(data$dt), length.out = 100),
  average_noise = med_noise,
  distance = 800
)

preds <- predict(m_ndt,
                 newdata = newdata,
                 type = 'link',
                 se.fit = T,
                 # terms = 's(dt)')
                 exclude = 's(array,station)')

newdata$pred <- m_ndt$family$linkinv(preds$fit)
newdata$lci <- m_ndt$family$linkinv(preds$fit - 1.96 * preds$se.fit)
newdata$uci <- m_ndt$family$linkinv(preds$fit + 1.96 * preds$se.fit)




ggplot(data = newdata) +
  geom_ribbon(aes(x = dt, ymin = lci, ymax = uci), fill = 'gray') +
  geom_line(aes(x=dt, y = pred)) +
  geom_rug(data = unique(data.table(data), by = 'dt'),
           aes(x = dt)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = 'ΔT (°C)', y = 'Predicted detectability') +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.margin = unit(c(0.05, 0.8, 0.05, 0.1), "cm"))




newdata <- data.frame(
  array = 'MDWEA',
  station = 'AN3',
  dt = 2.5,
  average_noise = seq(min(data$average_noise), max(data$average_noise), length.out = 100),
  distance = 800
)

preds <- predict(m_ndt,
                 newdata = newdata,
                 type = 'link',
                 se.fit = T,
                 # terms = 's(dt)')
                 exclude = 's(array,station)')

newdata$pred <- m_ndt$family$linkinv(preds$fit)
newdata$lci <- m_ndt$family$linkinv(preds$fit - 1.96 * preds$se.fit)
newdata$uci <- m_ndt$family$linkinv(preds$fit + 1.96 * preds$se.fit)




ggplot(data = newdata) +
  geom_ribbon(aes(x = average_noise, ymin = lci, ymax = uci), fill = 'gray') +
  geom_line(aes(x=average_noise, y = pred)) +
  geom_rug(data = unique(data.table(data), by = 'average_noise'),
           aes(x = average_noise)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = 'ΔT (°C)', y = 'Predicted detectability') +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.margin = unit(c(0.05, 0.8, 0.05, 0.1), "cm"))









# Tensor product prediction ----
# new_data_te <- expand.grid(
#   average_noise = seq(min(data$average_noise) -
#                         0.05 * (max(data$average_noise) - min(data$average_noise)),
#                       max(data$average_noise) +
#                         0.05 * (max(data$average_noise) - min(data$average_noise)),
#                       length.out = 100),
#   dt = seq(min(data$dt) -
#              0.01 * (max(data$dt) - min(data$dt)),
#            max(data$dt) +
#              0.01 * (max(data$dt) - min(data$dt)),
#            length.out = 100)
# )
# new_data_te$array <- 'Inner' # dummy input, as this will be excluded
# new_data_te$distance <- 800
#
#
# preds <- predict(mod_int, new_data_te, type = 'link', se  = T,
#               exclude = 's(array)')
# linkinv <- mod_int$family$linkinv
#
# new_data_te$pred <- linkinv(preds$fit)
# new_data_te$lci <- linkinv(preds$fit - 1.96 * preds$se.fit)
# new_data_te$uci <- linkinv(preds$fit + 1.96 * preds$se.fit)
#
# # Flag points that are too far from observed values
# too.far <- exclude.too.far(new_data_te$dt, new_data_te$average_noise,
#                            data[data$distance == 800,]$dt, data[data$distance == 800,]$average_noise,
#                            dist = 0.2)
#
# new_data_te <- new_data_te[!too.far,]
#
#
#
#
# preds <- predict(mod_int, type = 'link', se = T,
#                  exclude = 's(array)')
# data$pred_int <- linkinv(preds$fit)
# data$lci_int <- linkinv(preds$fit - 1.96 * preds$se.fit)
# data$uci_int <- linkinv(preds$fit + 1.96 * preds$se.fit)
# data$signif <- ifelse(data$uci_int < 0.5, 'less',
#                       ifelse(data$lci_int > 0.5, 'more', 'non'))
# data$signif <- factor(data$signif,
#                       levels = c('more', 'non', 'less'), ordered = T)
#
# # Observed @ 800m plus or minus mean model RMSE
# data <- data[data$distance == 800,]
# data$rmse_signif <- ifelse(data$freq < (0.5 - 0.187), 'less',
#                            ifelse(data$freq > (0.5 + 0.187), 'more', 'non'))
# data$rmse_signif <- factor(data$rmse_signif,
#                            levels = c('more', 'non', 'less'), ordered = T)
#
#
# # Tensor product visualization ----
# te_vis <- ggplot(data = new_data_te, aes(x = dt, y = average_noise, z = pred)) +
#   # Draw contours
#   geom_contour_filled() +
#   # Highlight 50%
#   geom_contour(breaks = 0.5, color = 'red', size = 1) +
#   # Draw points, colored by significance
#   geom_point(inherit.aes = F,
#              data = data,
#              aes(x = dt, y = average_noise, color = rmse_signif),
#              show.legend = F) +
#   geom_hline(yintercept = 231, linetype = 'dashed') +
#   scale_color_manual(values = c('turquoise', 'gray85', 'yellow4')) +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(x = NULL, y = 'Ambient noise at 69 kHz (mV)', fill = 'Pred. detectability') +
#   theme_bw() +
#   theme(axis.text.x = element_blank(),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.justification = c(1, 1),
#         legend.position = c(0.99, 0.99),
#         legend.background = element_rect(fill = 'gray96'),
#         panel.background = element_rect(fill = 'gray96'),
#         panel.grid = element_blank(),
#         plot.margin = unit(c(0.05, 0.8, 0.05, 0.1), "cm"))
#
# te_vis
#
# # Cross section visualization ----
# ##  hold noise constant at its median
# med_noise <- median(data[data$distance == 800,]$average_noise)
#
# new_data_slice <- data.frame(array = 'Inner',
#                        distance = 800,
#                        average_noise = med_noise,
#                        dt = seq(min(data$dt) -
#                                   0.01 * (max(data$dt) - min(data$dt)),
#                                 max(data$dt) +
#                                   0.01 * (max(data$dt) - min(data$dt)),
#                                 length.out = 100)
# )
#
# preds <- predict(mod_int, new_data_slice, type = 'link', se  = T,
#                  exclude = 's(array)')
#
# new_data_slice$pred <- linkinv(preds$fit)
# new_data_slice$lci <- linkinv(preds$fit - 1.96 * preds$se.fit)
# new_data_slice$uci <- linkinv(preds$fit + 1.96 * preds$se.fit)
#
# slice_subset <- data[data$distance == 800 &
#                        data$average_noise %between% c(floor(med_noise - 5),
#                                                       ceiling(med_noise + 5)),]
#
# slice_vis <- ggplot() +
#   geom_ribbon(data = new_data_slice, aes(x = dt, ymin = lci, ymax = uci), fill = 'gray') +
#   geom_line(data = new_data_slice, aes(x = dt, y = pred)) +
#
#   ### Do we want a point overlay?
#   # geom_point(inherit.aes = F,
#   #            data = data[data$average_noise %between% c(floor(med_noise - 5),
#   #                                                       ceiling(med_noise + 5)),],
#   #            aes(x = dt, y = freq, color = rmse_signif),
#   #            show.legend = F) +
#   scale_color_manual(values = c('turquoise', 'gray85', 'yellow4')) +
#
#   geom_rug(data = slice_subset,
#            aes(x = dt)) +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
#   labs(x = 'ΔT (°C)', y = 'Predicted detectability', fill = NULL) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12),
#         plot.margin = unit(c(0.05, 0.8, 0.05, 0.1), "cm"))
#
#
# library(patchwork)
# tiff('range test/manuscript/figures/Figure3.tif', compression = 'lzw', res = 300,
#      width = 170, height = 160, units = 'mm', pointsize = 6)
#
# te_vis / slice_vis + plot_layout(heights = c(2, 1)) & theme(axis.text = element_text(size = 10),
#                                                             axis.title = element_text(size = 12),
#                                                             legend.text = element_text(size = 10))
#
# dev.off()

