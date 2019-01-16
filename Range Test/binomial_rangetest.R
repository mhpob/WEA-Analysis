library(dplyr)
library(mgcv)

data <- readRDS('data and imports/rangetest_logit_binary.RDS') %>%
  mutate(array = as.factor(array),
         date_f = as.factor(date))

mod <- gamm(cbind(success, fail) ~ distance + s(dt, bs = "ts"),
            random = list(array = ~1, date_f = ~1),
            data = data, family = 'binomial', method = 'REML')
summary(mod)
anova(mod)


# To calculate D50 of fitted model solve:
# logit(p) ~ intercept + smoother coefficients * variable + distance coefficient * distance
# log10(p / [1 - p]) ~ intercept + sm.coef * variable + dist.coef * distance
# 0 ~ intercept + sm.coef * variable + dist.coef * D50
# D50 ~ -(intercept + sm.coef * variable) / dist.coef
#
# In words: negative sum of linear predictor matrix without the distance term,
# multiplied by their coefficients,
# then all divided by the distance coefficient.
#
# You can get the sum of the coefs * variable by multiplying the predicted
# linear prediction matrix with the corresponding coefficients.

# Create the linear prediction matrix for dt = 5; other terms don't matter, but
# remove the effect of random terms.
lpm <- predict(mod$gam, data.frame(dt = 10,
                               distance = 1,
                               date_f = '2018-12-02',
                               array = 'MD WEA'),
               type = 'lpmatrix', exclude = c('s(date_f)', 's(array)'))
d50 <- -(lpm[, -2] %*% coef(mod$gam)[-2]) / coef(mod$gam)[2]



# SE prediction ----
# Using the posterior distribution of the distance parameter, simulate large number
# of parameter values. This assumes that the parameter estimates are multivariate
# normally distributed. Modified from help documentation of predict.gam().
#
# "Population prediction intervals": https://stackoverflow.com/a/35589393/7496818

ci.pred <- function(n, obj, lpreds){
  # n = number of draws
  # obj = fitted model
  # lpreds = linear predictor matrix of new data

  param_reps <- MASS::mvrnorm(n, coef(obj), vcov(obj))
  estimate_reps <- rep(0, n)

  for (i in 1:n){
    estimate_reps[i] <- -(lpreds[, -2] %*% param_reps[i, -2]) / param_reps[i, 2]
  }

  mod_d50 <- -(lpreds[, -2] %*% coef(obj)[-2]) / coef(obj)[2]
  sd_estimate <- sqrt(var(estimate_reps))

  out <- c(d50 = mod_d50,
           uci = mod_d50 + 1.96 * sd_estimate,
           lci = mod_d50 - 1.96 * sd_estimate)
  out
}

ci.pred(10000, mod$gam, lpm)

