# Packages ----
library(ggplot2); library(tidyr); library(data.table); library(mgcv)



# Import data ----
data <- readRDS('data and imports/rangetest_logit_binary_pt0.RDS')
names(data) <- gsub(' ', '_', tolower(names(data)))
data$array <- as.factor(gsub(' ', '', data$array))
data <- data[data$distance > 0,]
data$freq <- data$success / (data$success + data$fail)



# Import winning model ----
##  see "mgcv model selection.R" for selection

mod_int <- gam(cbind(success, fail) ~
                      distance +
                      te(average_noise, dt) +
                      s(array, bs = 're'),
                    family = binomial(),
                    data = data,
                    method = 'REML',
                    control = gam.control(nthreads = 4))



# Tensor product prediction ----
new_data <- expand.grid(
  average_noise = seq(min(data$average_noise) -
                        0.05 * (max(data$average_noise) - min(data$average_noise)),
                      max(data$average_noise) +
                        0.05 * (max(data$average_noise) - min(data$average_noise)),
                      length.out = 100),
  dt = seq(min(data$dt) -
             0.01 * (max(data$dt) - min(data$dt)),
           max(data$dt) +
             0.01 * (max(data$dt) - min(data$dt)),
           length.out = 100)
)
new_data$array <- 'Inner' # dummy input, as this will be excluded
new_data$distance <- 800


preds <- predict(mod_int, new_data, type = 'link', se  = T,
              exclude = 's(array)')
linkinv <- mod_int$family$linkinv

new_data$pred <- linkinv(preds$fit)
new_data$lci <- linkinv(preds$fit - 1.96 * preds$se.fit)
new_data$uci <- linkinv(preds$fit + 1.96 * preds$se.fit)

# Flag points that are too far from observed values
too.far <- exclude.too.far(new_data$dt, new_data$average_noise,
                           data[data$distance == 800,]$dt, data[data$distance == 800,]$average_noise,
                           dist = 0.05)

new_data <- new_data[!too.far,]




preds <- predict(mod_int, type = 'link', se = T,
                 exclude = 's(array)')
data$pred_int <- linkinv(preds$fit)
data$lci_int <- linkinv(preds$fit - 1.96 * preds$se.fit)
data$uci_int <- linkinv(preds$fit + 1.96 * preds$se.fit)
data$signif <- ifelse(data$uci_int < 0.5, 'less',
                      ifelse(data$lci_int > 0.5, 'more', 'non'))
data$signif <- factor(data$signif,
                      levels = c('more', 'non', 'less'), ordered = T)

# Observed @ 800m plus or minus mean model RMSE
data <- data[data$distance == 800,]
data$rmse_signif <- ifelse(data$freq < (0.5 - 0.187), 'less',
                           ifelse(data$freq > (0.5 + 0.187), 'more', 'non'))
data$rmse_signif <- factor(data$rmse_signif,
                           levels = c('more', 'non', 'less'), ordered = T)


# Tensor product visualization ----
ggplot(data = new_data, aes(x = dt, y = average_noise, z = pred)) +
  # Draw contours
  geom_contour_filled() +
  # Highlight 50%
  geom_contour(breaks = 0.5, color = 'red', size = 1) +
  # Draw points, colored by significance
  geom_point(inherit.aes = F,
             data = data,
             aes(x = dt, y = average_noise, color = rmse_signif),
             show.legend = F) +
  scale_color_manual(values = c('turquoise', 'gray85', 'yellow4')) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'ΔT (°C)', y = 'Ambient noise at 69 kHz (mV)', fill = NULL) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = 'gray96'),
        panel.grid = element_blank())



# Cross section visualization ----
##  hold noise constant at its median
med_noise <- median(data[data$distance == 800,]$average_noise)

new_data <- data.frame(array = 'Inner',
                       distance = 800,
                       average_noise = med_noise,
                       dt = seq(min(data$dt) -
                                  0.01 * (max(data$dt) - min(data$dt)),
                                max(data$dt) +
                                  0.01 * (max(data$dt) - min(data$dt)),
                                length.out = 100)
)

preds <- predict(mod_int, new_data, type = 'link', se  = T,
                 exclude = 's(array)')

new_data$pred <- linkinv(preds$fit)
new_data$lci <- linkinv(preds$fit - 1.96 * preds$se.fit)
new_data$uci <- linkinv(preds$fit + 1.96 * preds$se.fit)

slice_subset <- data[data$distance == 800 &
                       data$average_noise %between% c(floor(med_noise - 5),
                                                      ceiling(med_noise + 5)),]

ggplot() +
  geom_ribbon(data = new_data, aes(x = dt, ymin = lci, ymax = uci), fill = 'gray') +
  geom_line(data = new_data, aes(x = dt, y = pred)) +

  ### Do we want a point overlay?
  geom_point(inherit.aes = F,
             data = data[data$average_noise %between% c(floor(med_noise - 5),
                                                        ceiling(med_noise + 5)),],
             aes(x = dt, y = freq, color = rmse_signif),
             show.legend = F) +
  scale_color_manual(values = c('turquoise', 'gray85', 'yellow4')) +

  geom_rug(data = slice_subset,
           aes(x = dt)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = 'ΔT (°C)', y = 'Predicted probability of detection', fill = NULL) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12))
