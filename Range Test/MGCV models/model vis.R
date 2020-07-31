# Packages ----
library(ggplot2); library(data.table); library(mgcv)



# Import data ----
data <- readRDS('data and imports/rangetest_logit_binary_pt0.RDS')
names(data) <- gsub(' ', '_', tolower(names(data)))
data$array <- as.factor(gsub(' ', '', data$array))
data <- data[data$distance > 0,]




# Import winning models ----
##  see "mgcv model selection.R" for selection
mod_dndt <- gam(cbind(success, fail) ~
                  distance +
                  s(average_noise) +
                  s(dt) +
                  s(array, bs = 're'),
                family = binomial(),
                data = data,
                method = 'REML',
                control = gam.control(nthreads = 4))

mod_int <- gam(cbind(success, fail) ~
                      distance +
                      te(average_noise, dt) +
                      s(array, bs = 're'),
                    family = binomial(),
                    data = data,
                    method = 'REML',
                    control = gam.control(nthreads = 4))



# ~ distance + te(noise, dt) ----
## Prediction
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
                           dist = 0.2)

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

## Vis
ggplot(data = new_data, aes(x = dt, y = average_noise, z = pred)) +
  # Draw contours
  geom_contour_filled() +
  # Highlight 50%
  geom_contour(breaks = 0.5, color = 'red', size = 1) +
  # Draw points, colored by significance
  geom_point(inherit.aes = F,
             data = data[data$distance == 800,],
             aes(x = dt, y = average_noise, color = signif),
             show.legend = F) +
  scale_color_manual(values = c('turquoise', 'gray85', 'yellow4')) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'ΔT', y = 'Noise at 69 kHz', fill = NULL) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12))


# ~ distance + noise + dt ----
## Prediction



## Vis