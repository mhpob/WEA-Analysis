library(dplyr)
library(mgcv)

data <- readRDS('data and imports/rangetest_logit_binary.RDS') %>%
  mutate(array = as.factor(array),
         date_f = as.factor(date))

mod <- gam(cbind(success, fail) ~ distance + s(dt, bs = "ts") +
            s(array, bs = 're') + s(date_f, bs = 're'),
         data = data, family = 'binomial', method = 'REML')
summary(mod)
anova(mod)


# sum of linear predictor matrix without the distance term multiplied by their coefficients,
# divided by the distance coef
lpm <- predict(mod, data.frame(distance = 500,
                             dt = 5,
                             date_f = '2018-12-02',
                             array = 'MD WEA'),
               type = 'lpmatrix', exclude = c('s(date_f)', 's(array)'))
d50 <- -sum(lpm[1, -2] * coef(mod)[-2]) / coef(mod)[2]



mod2 <- gamm(cbind(success, fail) ~ distance + s(dt, bs = "ts"),
             random = list(date_f = ~1),
             data = data, family = 'binomial', method = 'REML')

lpm <- predict(mod2$gam, data.frame(distance = 1,
                               dt = 10),
               type = 'lpmatrix', exclude = c('s(date_f)'))
d50 <- -(lpm[, -2] %*% coef(mod2$gam)[-2]) / coef(mod2$gam)[2]



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

  param_reps <- MASS::mvrnorm(10000, coef(mod2$gam), vcov(mod2$gam))
  estimate_reps <- rep(0, 10000)

  for (i in 1:10000){
    estimate_reps[i] <- -(lpreds[, -2] %*% param_reps[i, -2]) / param_reps[i, 2]
  }

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

  mod_d50 <- -(lpreds[, -2] %*% coef(mod2$gam)[-2]) / coef(mod2$gam)[2]

  sd_estimate <- sqrt(var(estimate_reps))

  out <- c(d50 = mod_d50,
           uci = mod_d50 + 1.96 * sd_estimate,
           lci = mod_d50 - 1.96 * sd_estimate)
  out
}

ci.pred(10000, mod2$gam, lpm)

lpm <- predict(mod2$gam, data.frame(distance = 1,
                                    dt = -1),
               type = 'lpmatrix', exclude = c('s(date_f)'))

ci.pred(10000, mod2$gam, lpm)

lpm <- predict(mod2$gam, data.frame(distance = 1,
                                    dt = -0.5),
               type = 'lpmatrix', exclude = c('s(date_f)'))
ci.pred(10000, mod2$gam, lpm)
