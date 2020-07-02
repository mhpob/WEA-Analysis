# Packages

library(dplyr)
library(mgcv)

# Data to fit the model (range test stuff)

data <- readRDS('p:/obrien/biotelemetry/md wea habitat/wea-analysis/data and imports/rangetest_logit_binary.RDS') %>%
  mutate(array = as.factor(array),
         date_f = as.factor(date),
         cp_phase = ifelse(lubridate::month(date) %in% 5:8, 'p', 'a')) %>%
  rename(noise = `Average noise`,
         tilt = `Tilt angle`,
         bwt = `Average temperature`)



#### Fit the model
#If you decide to move array from a fixed to a random effect,
# just change to something like s(array, bs = 're'). I chose 5 knots to keep
# things from overfitting, and cubic splines over thin plate because thin plate
# really, really, ridiculously overfit the interaction term in areas where we
# had no data. You can change the main effect basis to whatever you want and it
# doesn't change much; cubic spline just performed better in the interaction.

model <- gam(cbind(success, fail) ~ distance + array +
               ti(dt, bs = 'cs', k = 5) +
               ti(noise, bs = 'cs', k = 5)+
               ti(dt, noise, bs = 'cs', k = 5),
             data = data,
             family = 'binomial',
             method = 'REML')

### Get data to predict
# Using the data we used to fit the model as the predicted data here. To use your
# data, you'll need a data frame with distance (dummy), array (dummy if array
# is switched to random effect), noise, and dt.

pred_data <- data %>%
  distinct(date_f, array, .keep_all = T) %>%
  filter(array == 'Inner') %>%
  mutate(distance = 0)

### Calculate the linear predictor matrix
# I don't actually know what I'm talking about and this is probably an incorrect
# way to look at it, but I view this as the number that is multiplied by the
# fitted model coefficients.
#
# The coefs of  the linear terms are always multiplied by 1 since that slope
# (the response) is constant no matter where the point is that we're trying to
# predict. The coefs of the curvy terms are multiplied by a whole bunch of funky
# numbers since the response is different depending on the value of the point we're
# trying to fit.

lpm <- predict(model, pred_data, type = 'lpmatrix')

# Note that since we specified 5 knots, there are 4 terms for the main effects
# (5 - 1) and 16 ((5 - 1) ^ 2) for the interaction.

colnames(lpm)

### Calculate D50
# log10(50% / (100% - 50%)) = intercept +
#                             dist_coef * D50 +
#                             array_coef * array +
#                             squiggly_coefs * squiggly_bits
# **ALGEBRAAAAaaaaa**
# 0 = intercept + dist_coef * D50 + array_coef * array + squiggly_coefs * squiggly_bits
# D50 = -(intercept + array_coef * array + squiggly_coefs * squiggly_bits) / dist_coef
#
# Note that the non-coefficient terms come from the linear predictor matrix

d50 <- (log10(0.5 / (1 - 0.5)) - (lpm[, -2] %*% coef(model)[-2])) / coef(model)[2]

# Generic function, assuming that "distance" is always the second term:

dXX <- function(model, lpmatrix, desired_efficiency){
  (log10(desired_efficiency / (1 - desired_efficiency)) -
     (lpmatrix[, -2] %*% coef(model)[-2])) / coef(model)[2]
}

# In english, that's:
# "log of desired percent divided by one minus desired percent...
#   ...minus the linear predictor matrix (without the distance term) times
#   the model coefficients (without the distance coefficient)...
#   ...divided by the distance coefficient

### Bind back to the data
# This will just spit out numbers in the same order as the data we wanted to
# predict. Should bind it back on to make sense of anything.

final_data <- cbind(pred_data, d50)


library(itsadug)

fvisgam(model, view = c('dt', 'noise'), transform = model$family$linkinv,
        cond = list(distance = 550, array = 'Inner'),
        ylab = 'Ambient noise (mV)', xlab = 'dT (C)',
        main = 'Noise - DeltaT Interaction')
points(data$dt, data$noise)

plot_smooth(model, 'noise', transform = model$family$linkinv,
            cond = list(distance = 550, dt = 0, array = 'Inner'),
            rug = T, xlab = 'Ambient noise (mV)', ylab = 'Detection Efficiency',
            ylim = c(0, 1))
plot_smooth(model, 'noise', transform = model$family$linkinv,
            cond = list(distance = 550, dt = 10, array = 'Inner'),
            rug = T, xlab = 'Ambient noise (mV)', ylab = 'Detection Efficiency',
            ylim = c(0, 1))

plot_smooth(model, 'dt', transform = model$family$linkinv,
            cond = list(distance = 550, noise = 200, array = 'Inner'),
            rug = T, xlab = 'dT (C)', ylab = 'Detection Efficiency',
            ylim = c(0, 1))
plot_smooth(model, 'dt', transform = model$family$linkinv,
            cond = list(distance = 550, noise = 300, array = 'Inner'),
            rug = T, xlab = 'dT (C)', ylab = 'Detection Efficiency',
            ylim = c(0, 1))
